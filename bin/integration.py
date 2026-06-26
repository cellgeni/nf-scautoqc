#!/usr/bin/env python
"""Run scVI integration with lower peak memory.

This trains scVI on non-doublet, non-cell-cycle, HVG-subset data, then re-reads
the input object and attaches cell-level integration outputs to the full-gene
object.
"""

import argparse
from contextlib import contextmanager
import copy
from datetime import datetime
import gc
import time

import h5py
import joblib
import matplotlib.pyplot as plt
import anndata
import numpy as np
import scanpy as sc
import scvi
from scipy import sparse


VERBOSE = False


def log(message, verbose=False):
    if verbose and not VERBOSE:
        return
    print(f"[{datetime.now().isoformat(timespec='seconds')}] {message}", flush=True)


@contextmanager
def log_step(message, verbose=False):
    start = time.monotonic()
    log(f"START {message}", verbose=verbose)
    try:
        yield
    finally:
        elapsed = time.monotonic() - start
        log(f"END {message} ({elapsed:.1f}s)", verbose=verbose)


def str_to_bool(value):
    """Parse Nextflow-style string booleans without treating "False" as true."""
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def close_backed(adata):
    if getattr(adata, "isbacked", False):
        adata.file.close()


def drop_unused_slots(adata):
    adata.raw = None
    adata.layers.clear()


def normalize_positions(indexer, size):
    if isinstance(indexer, slice):
        return np.arange(size, dtype=np.int64)[indexer]

    indexer = np.asarray(indexer)
    if indexer.dtype == bool:
        if indexer.size != size:
            raise ValueError(f"Boolean index has size {indexer.size}; expected {size}.")
        return np.flatnonzero(indexer)

    positions = indexer.astype(np.int64, copy=False)
    if positions.size and np.any(positions[1:] < positions[:-1]):
        raise ValueError("Fast CSR reader requires sorted row/column positions.")
    return positions


def read_csr_subset_from_h5ad(
    obj_path,
    obs_idx,
    var_idx,
    obs,
    var,
    block_size=10000,
    data_dtype=np.int32,
):
    """Read X/obs/var for a CSR-backed H5AD subset without AnnData view slicing."""
    with h5py.File(obj_path, "r") as handle:
        x_group = handle["X"]
        if x_group.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError("Fast reader only supports CSR-backed X.")

        n_obs, n_vars = map(int, x_group.attrs["shape"])
        obs_positions = normalize_positions(obs_idx, n_obs)
        var_positions = normalize_positions(var_idx, n_vars)
        log(
            "Fast CSR read requested: "
            f"{obs_positions.size} obs x {var_positions.size} vars "
            f"from source {n_obs} obs x {n_vars} vars",
            verbose=True,
        )

        column_map = np.full(n_vars, -1, dtype=np.int64)
        column_map[var_positions] = np.arange(var_positions.size, dtype=np.int64)

        data_ds = x_group["data"]
        indices_ds = x_group["indices"]
        indptr_ds = x_group["indptr"]

        new_indptr = np.empty(obs_positions.size + 1, dtype=np.int64)
        new_indptr[0] = 0
        cursor = 0
        out_row = 0

        while cursor < obs_positions.size:
            block_start = int(obs_positions[cursor])
            block_stop = min(block_start + block_size, n_obs)
            next_cursor = cursor
            while (
                next_cursor < obs_positions.size
                and obs_positions[next_cursor] < block_stop
            ):
                next_cursor += 1

            block_rows = obs_positions[cursor:next_cursor]
            block_indptr = indptr_ds[block_start : block_stop + 1]
            data_start = int(block_indptr[0])
            data_stop = int(block_indptr[-1])
            block_indices = indices_ds[data_start:data_stop]
            block_indptr = block_indptr - data_start
            column_keep = column_map[block_indices] >= 0

            for row in block_rows:
                rel_row = int(row) - block_start
                row_start = int(block_indptr[rel_row])
                row_stop = int(block_indptr[rel_row + 1])
                row_count = int(column_keep[row_start:row_stop].sum())
                new_indptr[out_row + 1] = new_indptr[out_row] + row_count
                out_row += 1

            cursor = next_cursor

        data = np.empty(int(new_indptr[-1]), dtype=data_dtype)
        indices = np.empty(int(new_indptr[-1]), dtype=np.int32)
        log(f"Fast CSR read output nnz: {int(new_indptr[-1])}", verbose=True)
        cursor = 0
        out_row = 0

        while cursor < obs_positions.size:
            block_start = int(obs_positions[cursor])
            block_stop = min(block_start + block_size, n_obs)
            next_cursor = cursor
            while (
                next_cursor < obs_positions.size
                and obs_positions[next_cursor] < block_stop
            ):
                next_cursor += 1

            block_rows = obs_positions[cursor:next_cursor]
            block_indptr = indptr_ds[block_start : block_stop + 1]
            data_start = int(block_indptr[0])
            data_stop = int(block_indptr[-1])
            block_indices = indices_ds[data_start:data_stop]
            block_data = data_ds[data_start:data_stop]
            block_indptr = block_indptr - data_start
            column_keep = column_map[block_indices] >= 0
            entry_keep = np.zeros(block_indices.size, dtype=bool)
            out_start = int(new_indptr[out_row])

            for row in block_rows:
                rel_row = int(row) - block_start
                row_start = int(block_indptr[rel_row])
                row_stop = int(block_indptr[rel_row + 1])
                row_keep = column_keep[row_start:row_stop]
                if row_keep.any():
                    entry_keep[row_start:row_stop] = row_keep
                out_row += 1

            out_stop = int(new_indptr[out_row])
            if out_stop > out_start:
                data[out_start:out_stop] = block_data[entry_keep].astype(
                    data_dtype, copy=False
                )
                indices[out_start:out_stop] = column_map[
                    block_indices[entry_keep]
                ].astype(np.int32, copy=False)

            cursor = next_cursor

    subset_x = sparse.csr_matrix(
        (data, indices, new_indptr),
        shape=(obs_positions.size, var_positions.size),
    )
    return anndata.AnnData(
        X=subset_x,
        obs=obs.iloc[obs_positions].copy(),
        var=var.iloc[var_positions].copy(),
    )


def read_minimal_view(source_ad, obs_idx, var_idx):
    """Read only X/obs/var from a backed AnnData view into memory."""
    if getattr(source_ad, "isbacked", False) and source_ad.filename is not None:
        try:
            return read_csr_subset_from_h5ad(
                source_ad.filename,
                obs_idx,
                var_idx,
                source_ad.obs,
                source_ad.var,
            )
        except ValueError as error:
            log(f"Falling back to AnnData backed slicing: {error}")

    with log_step("AnnData backed slicing fallback"):
        view = source_ad[obs_idx, var_idx]
        return anndata.AnnData(
            X=view.X.copy(),
            obs=view.obs.copy(),
            var=view.var.copy(),
        )


def get_training_masks(source_ad, from_scautoqc):
    if from_scautoqc:
        keep_cells = ~source_ad.obs["stringent_doublet"].to_numpy()
        keep_genes = ~source_ad.var["cc"].to_numpy()
    else:
        keep_cells = np.ones(source_ad.n_obs, dtype=bool)
        keep_genes = np.ones(source_ad.n_vars, dtype=bool)

    log(
        "Training masks: "
        f"{int(keep_cells.sum())}/{keep_cells.size} cells kept, "
        f"{int(keep_genes.sum())}/{keep_genes.size} genes kept"
    )
    return keep_cells, keep_genes


def select_cell_nuclei_hvg_union(source_ad, keep_cells, keep_genes, n_top_genes):
    """Select union HVGs and retain per-modality HVG annotations."""
    group_key = "cell_or_nuclei"
    if group_key not in source_ad.obs:
        raise KeyError(f"Expected '{group_key}' in adata.obs for HVG selection.")

    labels = source_ad.obs[group_key].astype(str)
    retained_labels = labels[keep_cells]
    unexpected_labels = sorted(set(retained_labels) - {"cell", "nuclei"})
    if unexpected_labels:
        raise ValueError(
            f"Expected '{group_key}' to contain only 'cell' or 'nuclei'; "
            f"found: {unexpected_labels}"
        )

    hvg_by_label = {
        "cell": np.zeros(source_ad.n_vars, dtype=bool),
        "nuclei": np.zeros(source_ad.n_vars, dtype=bool),
    }
    for label in ("cell", "nuclei"):
        mask = keep_cells & (labels == label).to_numpy()
        n_obs = int(mask.sum())
        if n_obs == 0:
            log(f"No {label} observations found for HVG selection; skipping.")
            continue

        log(f"{label} HVG selection will use {n_obs} cells", verbose=True)
        with log_step(f"Read {label} matrix for HVG selection"):
            modality_ad = read_minimal_view(source_ad, mask, keep_genes)
        drop_unused_slots(modality_ad)
        log(
            f"{label} HVG AnnData loaded: "
            f"{modality_ad.n_obs} obs x {modality_ad.n_vars} vars, "
            f"nnz={modality_ad.X.nnz if sparse.issparse(modality_ad.X) else 'dense'}",
            verbose=True,
        )
        with log_step(f"Compute {label} highly variable genes"):
            hvg_df = sc.pp.highly_variable_genes(
                modality_ad,
                flavor="seurat_v3",
                n_top_genes=min(n_top_genes, modality_ad.n_vars),
                subset=False,
                inplace=False,
            )
        selected = hvg_df["highly_variable"].to_numpy(dtype=bool)
        selected_names = modality_ad.var_names[selected]
        hvg_by_label[label] = source_ad.var_names.isin(selected_names)
        log(f"Selected {int(selected.sum())} HVGs from {n_obs} {label} observations.")
        del modality_ad, hvg_df, selected, selected_names, mask
        gc.collect()

    hvg_genes = hvg_by_label["cell"] | hvg_by_label["nuclei"]
    if not hvg_genes.any():
        raise ValueError("No cell or nuclei observations found for HVG selection.")

    hvg_group = np.full(source_ad.n_vars, "not_highly_variable", dtype=object)
    cell_only = hvg_by_label["cell"] & ~hvg_by_label["nuclei"]
    nuclei_only = hvg_by_label["nuclei"] & ~hvg_by_label["cell"]
    both = hvg_by_label["cell"] & hvg_by_label["nuclei"]
    hvg_group[cell_only] = "cell_only"
    hvg_group[nuclei_only] = "nuclei_only"
    hvg_group[both] = "both"

    hvg_info = source_ad.var.iloc[:, 0:0].copy()
    hvg_info["highly_variable"] = hvg_genes
    hvg_info["highly_variable_group"] = hvg_group

    log(
        "Using union HVGs for scVI: "
        f"{int(hvg_genes.sum())} total, "
        f"{int(cell_only.sum())} cell_only, "
        f"{int(nuclei_only.sum())} nuclei_only, "
        f"{int(both.sum())} both."
    )
    return hvg_genes, hvg_info


def select_global_hvgs(source_ad, keep_cells, keep_genes, n_top_genes):
    """Select HVGs from all retained cells regardless of cell/nuclei status."""
    log("Global HVG selection will use all retained cells.")
    with log_step("Read matrix for global HVG selection"):
        hvg_ad = read_minimal_view(source_ad, keep_cells, keep_genes)
    drop_unused_slots(hvg_ad)
    log(
        "Global HVG AnnData loaded: "
        f"{hvg_ad.n_obs} obs x {hvg_ad.n_vars} vars, "
        f"nnz={hvg_ad.X.nnz if sparse.issparse(hvg_ad.X) else 'dense'}",
        verbose=True,
    )
    with log_step("Compute global highly variable genes"):
        hvg_df = sc.pp.highly_variable_genes(
            hvg_ad,
            flavor="seurat_v3",
            n_top_genes=min(n_top_genes, hvg_ad.n_vars),
            subset=False,
            inplace=False,
        )

    selected = hvg_df["highly_variable"].to_numpy(dtype=bool)
    selected_names = hvg_ad.var_names[selected]
    hvg_genes = source_ad.var_names.isin(selected_names)

    hvg_group = np.full(source_ad.n_vars, "not_highly_variable", dtype=object)
    hvg_group[hvg_genes] = "global"
    hvg_info = source_ad.var.iloc[:, 0:0].copy()
    hvg_info["highly_variable"] = hvg_genes
    hvg_info["highly_variable_group"] = hvg_group

    log(f"Using global HVGs for scVI: {int(hvg_genes.sum())} total.")
    del hvg_ad, hvg_df, selected, selected_names
    gc.collect()
    return hvg_genes, hvg_info


def select_hvgs(source_ad, keep_cells, keep_genes, n_top_genes, hvg_strategy):
    if hvg_strategy == "cell_nuclei_union":
        return select_cell_nuclei_hvg_union(source_ad, keep_cells, keep_genes, n_top_genes)
    if hvg_strategy == "global":
        return select_global_hvgs(source_ad, keep_cells, keep_genes, n_top_genes)
    raise ValueError(f"Unknown HVG strategy: {hvg_strategy}")


def read_training_object(obj_path, from_scautoqc, n_top_genes, hvg_strategy):
    """Load only retained cells and selected HVGs for scVI training."""
    with log_step(f"Open input object backed: {obj_path}"):
        source_ad = sc.read_h5ad(obj_path, backed="r")
    try:
        keep_cells, keep_genes = get_training_masks(source_ad, from_scautoqc)
        with log_step(f"Select HVGs using strategy: {hvg_strategy}"):
            hvg_genes, hvg_info = select_hvgs(
                source_ad, keep_cells, keep_genes, n_top_genes, hvg_strategy
            )
        with log_step("Read final SCVI training matrix"):
            train_ad = read_minimal_view(source_ad, keep_cells, hvg_genes)
    finally:
        close_backed(source_ad)
        log("Closed backed input object", verbose=True)

    drop_unused_slots(train_ad)
    train_hvg_info = hvg_info.reindex(train_ad.var_names)
    for column in hvg_info.columns:
        train_ad.var[column] = train_hvg_info[column].to_numpy()
    train_ad.uns["_integration_obs_positions"] = np.flatnonzero(keep_cells)
    log(
        "SCVI training AnnData ready: "
        f"{train_ad.n_obs} obs x {train_ad.n_vars} vars, "
        f"nnz={train_ad.X.nnz if sparse.issparse(train_ad.X) else 'dense'}"
    )
    return train_ad, hvg_info


def stash_cell_level_outputs(adata, hvg_info):
    """Keep integration outputs before freeing the train object."""
    return {
        "obs_names": adata.obs_names.copy(),
        "obs_positions": adata.uns["_integration_obs_positions"].copy(),
        "var_hvg": hvg_info.copy(),
        "obs": adata.obs.copy(),
        "obsm": dict(adata.obsm.items()),
        "obsp": dict(adata.obsp.items()),
        "uns": {
            key: copy.deepcopy(adata.uns[key])
            for key in ("neighbors", "umap")
            if key in adata.uns
        },
    }


def attach_outputs_to_full_object(obj_path, outputs):
    """Re-read the original object, keep all genes, and attach integration outputs."""
    with log_step("Open full object backed for final output"):
        full_ad = sc.read_h5ad(obj_path, backed="r")
    try:
        with log_step("Read full-gene X for retained cells"):
            full_view = full_ad[outputs["obs_positions"], :]
            out_ad = anndata.AnnData(
                X=full_view.X.copy(),
                obs=outputs["obs"],
                var=full_view.var.copy(),
                uns=copy.deepcopy(outputs["uns"]),
            )
    finally:
        close_backed(full_ad)
        log("Closed backed full object")

    out_ad.obsm.update(outputs["obsm"])
    out_ad.obsp.update(outputs["obsp"])
    hvg_info = outputs["var_hvg"].reindex(out_ad.var_names)
    for column in hvg_info.columns:
        if hvg_info[column].dtype == bool:
            out_ad.var[column] = hvg_info[column].fillna(False).to_numpy(dtype=bool)
        else:
            out_ad.var[column] = hvg_info[column].fillna("not_highly_variable").to_numpy()
    log(
        "Final AnnData assembled: "
        f"{out_ad.n_obs} obs x {out_ad.n_vars} vars, "
        f"nnz={out_ad.X.nnz if sparse.issparse(out_ad.X) else 'dense'}"
    )

    return out_ad


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--obj", required=True, help="path of the pooled h5ad")
    parser.add_argument("--batch", required=True, help="batch key")
    parser.add_argument(
        "--n_top_genes",
        default=5000,
        type=int,
        help="number of top genes to use for HVG selection",
    )
    parser.add_argument(
        "--hvg_strategy",
        choices=("cell_nuclei_union", "global"),
        default="cell_nuclei_union",
        help=(
            "HVG selection strategy. cell_nuclei_union calculates HVGs separately "
            "for cell and nuclei observations and uses the union; global ignores "
            "cell_or_nuclei and calculates one HVG set from all retained cells."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print detailed matrix-reading and cleanup logs",
    )
    parser.add_argument(
        "--from_scautoqc",
        default=None,
        help="whether the input object is from scAutoQC",
    )
    args = parser.parse_args()
    global VERBOSE
    VERBOSE = args.verbose
    log(
        "integration.py started with "
        f"obj={args.obj}, batch={args.batch}, n_top_genes={args.n_top_genes}, "
        f"hvg_strategy={args.hvg_strategy}, from_scautoqc={args.from_scautoqc}"
    )

    arches_params = dict(
        n_layers=2,
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
    )

    from_scautoqc = str_to_bool(args.from_scautoqc)

    with log_step("Prepare training object"):
        train_ad, hvg_info = read_training_object(
            args.obj, from_scautoqc, args.n_top_genes, args.hvg_strategy
        )
    gc.collect()

    with log_step("Set up SCVI AnnData registry"):
        scvi.model.SCVI.setup_anndata(
            train_ad,
            batch_key=args.batch,
            continuous_covariate_keys=["log1p_n_counts", "percent_mito"],
        )

    with log_step("Initialise SCVI model"):
        vae = scvi.model.SCVI(train_ad, n_latent=20, **arches_params)

    with log_step("Train SCVI model"):
        vae.train(
            train_size=0.9,
            early_stopping_patience=45,
            max_epochs=100,
            batch_size=256,
            limit_train_batches=100,
        )

    with log_step("Write scVI model pickle"):
        joblib.dump(
            vae,
            "scvi_model.pkl",
        )

    with log_step("Write ELBO training plot"):
        plt.plot(vae.history["elbo_train"])[0].figure.savefig("elbo_training.png")
        plt.close("all")

    with log_step("Calculate latent representation"):
        train_ad.obsm["X_scvi"] = vae.get_latent_representation()

    with log_step("Write X_scVI embedding pickle"):
        joblib.dump(
            train_ad.obsm["X_scvi"],
            "Xscvi_embed.pkl",
        )

    with log_step("Compute neighbors from X_scvi"):
        sc.pp.neighbors(train_ad, use_rep="X_scvi")
    with log_step("Compute UMAP"):
        sc.tl.umap(train_ad, min_dist=0.3)

    with log_step("Stash cell-level integration outputs"):
        outputs = stash_cell_level_outputs(train_ad, hvg_info)

    del train_ad, hvg_info, vae
    gc.collect()

    with log_step("Attach outputs to full-gene object"):
        out_ad = attach_outputs_to_full_object(args.obj, outputs)
    del outputs
    gc.collect()

    with log_step("Write scautoqc_integrated.h5ad"):
        out_ad.write(
            "scautoqc_integrated.h5ad",
            compression="gzip",
        )
    log("integration.py finished")


if __name__ == "__main__":
    main()
