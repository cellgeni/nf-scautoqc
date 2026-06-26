#!/usr/bin/env python
"""
Gather STARsolo/Cell Ranger, Velocyto, and CellBender outputs into one h5ad.

This is a refactored alternative to gather_matrices.py intended for side-by-side
comparison. The original script is left unchanged.
"""

import argparse
import logging
import os
import signal
import sys
import warnings

import anndata
import scanpy as sc
import sctk as sk

signal.signal(signal.SIGPIPE, signal.SIG_DFL)
warnings.simplefilter(action="ignore", category=FutureWarning)


def _has_any(path, names):
    return any(os.path.exists(os.path.join(path, name)) for name in names)


def is_10x_mtx_dir(path):
    """
    Accept both gzipped and plain 10x/STARsolo matrix triplets.
    Also accept legacy 'genes.tsv(.gz)' in place of 'features.tsv(.gz)'.
    """
    return (
        _has_any(path, ("barcodes.tsv.gz", "barcodes.tsv"))
        and _has_any(path, ("matrix.mtx.gz", "matrix.mtx"))
        and _has_any(path, ("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"))
    )


def has_input(path):
    if path is None:
        return False
    if isinstance(path, str) and path.strip() in ("", "[]", "None"):
        return False
    return os.path.exists(path)


def is_h5(path):
    return os.path.isfile(path) and os.path.realpath(path).endswith(".h5")


def resolve_gene_input(cr_gene_input, gather_mode, gene_matrix_level):
    """
    Resolve which matrix path to read based on input type and desired level.

    gene_matrix_level:
      - filtered: prefer filtered matrix
      - raw: prefer raw matrix
      - auto: filtered for starsolo mode, raw for cellbender mode
    """
    if gene_matrix_level not in ("filtered", "raw", "auto"):
        raise ValueError("gene_matrix_level must be one of: filtered, raw, auto")

    if not os.path.exists(cr_gene_input):
        raise FileNotFoundError(cr_gene_input)

    if gene_matrix_level == "auto":
        prefer_raw = gather_mode == "cellbender"
    else:
        prefer_raw = gene_matrix_level == "raw"

    if os.path.isdir(cr_gene_input) and is_10x_mtx_dir(cr_gene_input):
        input_mode = "starsolo"
        if prefer_raw:
            raw_candidate = os.path.join(os.path.dirname(os.path.realpath(cr_gene_input)), "raw")
            if is_10x_mtx_dir(raw_candidate):
                return input_mode, raw_candidate
        return input_mode, cr_gene_input

    input_mode = "cellranger"
    if is_h5(cr_gene_input):
        if prefer_raw:
            raw_h5 = os.path.join(os.path.dirname(os.path.realpath(cr_gene_input)), "raw_feature_bc_matrix.h5")
            if os.path.exists(raw_h5):
                return input_mode, raw_h5
            raw_mtx = os.path.join(os.path.dirname(os.path.realpath(cr_gene_input)), "raw_feature_bc_matrix")
            if is_10x_mtx_dir(raw_mtx):
                return input_mode, raw_mtx
        return input_mode, os.path.realpath(cr_gene_input)

    if os.path.isdir(cr_gene_input):
        if prefer_raw:
            raw_h5 = os.path.join(os.path.dirname(os.path.realpath(cr_gene_input)), "raw_feature_bc_matrix.h5")
            if os.path.exists(raw_h5):
                return input_mode, raw_h5
            raw_mtx = os.path.join(os.path.dirname(os.path.realpath(cr_gene_input)), "raw_feature_bc_matrix")
            if is_10x_mtx_dir(raw_mtx):
                return input_mode, raw_mtx
        return input_mode, cr_gene_input

    raise ValueError(f"Unrecognized gene matrix input: {cr_gene_input}")


def resolve_gather_mode(gather_mode, cell_or_nuclei):
    if gather_mode in ("starsolo", "cellbender"):
        return gather_mode
    if gather_mode != "auto":
        raise ValueError("gather_mode must be one of: auto, starsolo, cellbender")

    label = (cell_or_nuclei or "cell").strip().lower()
    if label == "cell":
        return "cellbender"
    if label == "nuclei":
        return "starsolo"
    raise ValueError("gather_mode=auto requires cell_or_nuclei to be either cell or nuclei")


def read_gene_matrix(resolved_path, input_mode):
    if is_h5(resolved_path):
        ad = sc.read_10x_h5(resolved_path)
    else:
        ad = sc.read_10x_mtx(resolved_path)

    if input_mode == "cellranger":
        ad.obs_names = ad.obs_names.str.replace("-1$", "", regex=True)
    standardize_var_names(ad)
    return ad


def standardize_var_names(ad):
    if "gene_name" not in ad.var:
        ad.var["gene_name"] = ad.var_names
    ad.var_names_make_unique()
    return ad


def read_cellbender_matrix(path):
    ad = sk.read_cellbender(path)
    return standardize_var_names(ad)


def read_velocyto_matrix(path):
    ad = sk.read_velocyto(os.path.realpath(path))
    return standardize_var_names(ad)


def align_vars_to_reference(source_ad, reference_var_names, matrix_label):
    if not source_ad.var_names.is_unique:
        raise ValueError(f"{matrix_label} matrix has duplicated feature names after standardization.")

    missing = reference_var_names.difference(source_ad.var_names)
    extra = source_ad.var_names.difference(reference_var_names)
    if len(missing) or len(extra):
        raise ValueError(
            f"{matrix_label} features do not match main matrix "
            f"({len(missing)} missing, {len(extra)} extra)."
        )

    return source_ad[:, reference_var_names]


def align_obs(main_ad, raw_ad, velo_ad=None):
    selected = main_ad.obs_names.intersection(raw_ad.obs_names)
    if len(selected) == 0:
        raise ValueError("No overlapping cell barcodes between main matrix and raw matrix.")

    main_idx = main_ad.obs_names.get_indexer(selected)
    raw_idx = raw_ad.obs_names.get_indexer(selected)

    if (main_idx < 0).any() or (raw_idx < 0).any():
        raise ValueError("Barcode alignment failed for main/raw matrices.")

    if velo_ad is not None:
        selected = selected.intersection(velo_ad.obs_names)
        if len(selected) == 0:
            raise ValueError("No overlapping cell barcodes after including Velocyto matrix.")
        main_idx = main_ad.obs_names.get_indexer(selected)
        raw_idx = raw_ad.obs_names.get_indexer(selected)
        velo_idx = velo_ad.obs_names.get_indexer(selected)
        if (velo_idx < 0).any():
            raise ValueError("Barcode alignment failed for Velocyto matrix.")
        return selected, main_idx, raw_idx, velo_idx

    return selected, main_idx, raw_idx, None


def gather_matrices(
    cr_gene_filtered_mtx,
    cr_velo_filtered_mtx,
    cb_filtered_h5,
    gather_mode,
    cell_or_nuclei="cell",
    gene_matrix_level="auto",
):
    requested_gather_mode = gather_mode
    gather_mode = resolve_gather_mode(gather_mode, cell_or_nuclei)
    input_mode, resolved_gene_path = resolve_gene_input(
        cr_gene_filtered_mtx, gather_mode=gather_mode, gene_matrix_level=gene_matrix_level
    )
    cr_gene_ad = read_gene_matrix(resolved_gene_path, input_mode=input_mode)
    logging.info(
        "cr_gene matrix loaded (input=%s, gather_mode=%s, requested_gather_mode=%s, cell_or_nuclei=%s, level=%s, path=%s)",
        input_mode,
        gather_mode,
        requested_gather_mode,
        cell_or_nuclei,
        gene_matrix_level,
        resolved_gene_path,
    )

    cr_velo_ad = None
    if has_input(cr_velo_filtered_mtx):
        cr_velo_ad = read_velocyto_matrix(cr_velo_filtered_mtx)
        logging.info("cr_velo matrix loaded")

    cb_gene_ad = None
    if has_input(cb_filtered_h5):
        cb_gene_ad = read_cellbender_matrix(cb_filtered_h5)
        logging.info("cellbender matrix loaded")

    if cb_gene_ad is not None:
        if gather_mode == "starsolo":
            main_ad = cb_gene_ad[cb_gene_ad.obs_names.intersection(cr_gene_ad.obs_names)].copy()
        elif gather_mode == "cellbender":
            main_ad = cb_gene_ad
        else:
            raise ValueError("gather_mode must be 'starsolo' or 'cellbender'")
    else:
        main_ad = cr_gene_ad

    selected_obs, main_idx, raw_idx, velo_idx = align_obs(main_ad, cr_gene_ad, cr_velo_ad)
    cr_gene_ad = align_vars_to_reference(cr_gene_ad, main_ad.var_names, "raw")
    if cr_velo_ad is not None:
        cr_velo_ad = align_vars_to_reference(cr_velo_ad, main_ad.var_names, "Velocyto")

    layers = {
        "raw": cr_gene_ad.X[raw_idx, :],
    }

    if cr_velo_ad is not None:
        layers["spliced"] = cr_velo_ad.X[velo_idx, :]
        layers["unspliced"] = cr_velo_ad.layers["unspliced"][velo_idx, :]
        layers["ambiguous"] = cr_velo_ad.layers["ambiguous"][velo_idx, :]

    ad = anndata.AnnData(
        X=main_ad.X[main_idx, :],
        obs=main_ad.obs.iloc[main_idx].copy(),
        var=main_ad.var.copy(),
        layers=layers,
    )
    ad.obs_names = selected_obs

    for layer in ("spliced", "unspliced", "ambiguous"):
        if layer in ad.layers and hasattr(ad.layers[layer], "eliminate_zeros"):
            ad.layers[layer].eliminate_zeros()

    return ad


def main(args):
    logging.info(args)

    if not os.path.exists(args.cr_gene):
        raise FileNotFoundError(args.cr_gene)

    ad = gather_matrices(
        cr_gene_filtered_mtx=args.cr_gene,
        cr_velo_filtered_mtx=args.cr_velo,
        cb_filtered_h5=args.cb_h5,
        gather_mode=args.gather_mode,
        cell_or_nuclei=args.cell_or_nuclei,
        gene_matrix_level=args.gene_matrix_level,
    )

    if args.cell_or_nuclei is not None:
        ad.uns["cell_or_nuclei"] = args.cell_or_nuclei
    else:
        logging.info("No cell_or_nuclei provided; defaulting to 'cell'")
        ad.uns["cell_or_nuclei"] = "cell"

    output_h5ad = "gene_velo_cellbender.filtered.h5ad"
    ad.write(output_h5ad, compression="gzip")
    logging.info("done")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cr_gene", required=True, help="Path to gene matrix input")
    parser.add_argument("--cr_velo", default=None, const=None, nargs="?", help="Path to velocyto filtered folder")
    parser.add_argument("--cb_h5", default=None, const=None, nargs="?", help="Path to cellbender h5")
    parser.add_argument("--cell_or_nuclei", default=None, const=None, nargs="?", help="cell or nuclei label")
    parser.add_argument(
        "--gather_mode",
        default="cellbender",
        choices=["auto", "starsolo", "cellbender"],
        help="Cell set policy when CellBender is present. auto uses cellbender for cell and starsolo for nuclei.",
    )
    parser.add_argument(
        "--gene_matrix_level",
        default="auto",
        choices=["auto", "filtered", "raw"],
        help="Which expression matrix level to use for cr_gene input",
    )
    args = parser.parse_args()

    try:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s; %(levelname)s; %(funcName)s; %(message)s",
            datefmt="%y-%m-%d %H:%M:%S",
        )
        main(args)
    except KeyboardInterrupt:
        logging.warning("Interrupted")
        sys.exit(1)
