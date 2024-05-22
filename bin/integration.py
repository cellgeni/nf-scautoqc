#!/usr/bin/env python
"""Add metadata to pooled object

Usage: sample_qc.py --obj pooled_obj.h5ad 

Options:
  --obj        path of the pooled h5ad
  --batch      batch key
"""
def run_mde(
    data,
    n_neighbors=15,
    random_state=0,
    **kwargs,
):
    """
    Util to run :func:`pymde.preserve_neighbors` for visualization of scvi-tools embeddings.

    Parameters
    ----------
    data
        The data of shape (n_obs, k), where k is typically defined by one of the models
        in scvi-tools that produces an embedding (e.g., :class:`~scvi.model.SCVI`.)
    n_neighbors
        Number of nearest neighbors, if set to None, a sensible number will be chosen according
        to the number of observations in `data`.
    random_state
        Random seed passed to mde for reproducibility
    kwargs
        Keyword args to :func:`pymde.preserve_neighbors`
    Returns
    -------
    The pymde embedding, defaults to two dimensions.

    Notes
    -----
    This function is a modification of scvi.model.utils.mde().

    If you use this function in your research please cite:

    Agrawal, Akshay, Alnur Ali, and Stephen Boyd. "Minimum-distortion embedding." arXiv preprint arXiv:2103.02559 (2021).
    """
    try:
        import pymde
        import torch
    except ImportError:
        raise ImportError("Please install pymde package via `pip install pymde`")

    if isinstance(data, pd.DataFrame):
        data = data.values

    device = "cuda"

    _kwargs = dict(
        embedding_dim=2,
        constraint=pymde.Standardized(),
        repulsive_fraction=0.5,
        verbose=False,
        device=device,
        n_neighbors=n_neighbors,
    )
    _kwargs.update(kwargs)

    pymde.seed(random_state)
    mde = pymde.preserve_neighbors(data, **_kwargs)

    return mde

def mde_embed(mde):
    emb = mde.embed(verbose=False)

    emb = emb.cpu().numpy()

    return emb
    
import scanpy as sc
import pandas as pd
import scvi
import sctk as sk
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
import anndata
import gc
import argparse
import joblib

my_parser = argparse.ArgumentParser()
my_parser.add_argument("--obj", default=None, help="path of the pooled h5ad")
my_parser.add_argument("--batch", default=None, help="batch key")
args = my_parser.parse_args()

arches_params = dict(
    n_layers=2,
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
)

pooled_ad1 = sc.read(args.obj)

pooled_ad1 = pooled_ad1[~pooled_ad1.obs.stringent_doublet].copy()

gc.collect()

pooled_ad1.raw = pooled_ad1
pooled_ad1 = pooled_ad1[:, ~pooled_ad1.var["cc"]].copy()
gc.collect()

sc.pp.highly_variable_genes(
    pooled_ad1,
    flavor="seurat_v3",
    n_top_genes=7500,
    subset=True,
)

scvi.model.SCVI.setup_anndata(
    pooled_ad1,
    batch_key=args.batch,
    continuous_covariate_keys=["log1p_n_counts", "percent_mito"],
)

vae = scvi.model.SCVI(pooled_ad1, n_latent=20, **arches_params)

vae.train(
    train_size=0.9,
    early_stopping_patience=45,
    max_epochs=400,
    batch_size=256,
    limit_train_batches=100,
    use_gpu=True,
)

joblib.dump(
    vae,
    f"scvi_model.pkl",
)

plt.plot(vae.history["elbo_train"])[0].figure.savefig("elbo_training.png")

pooled_ad1.obsm["X_scvi"] = vae.get_latent_representation()

joblib.dump(
    pooled_ad1.obsm["X_scvi"],
    "Xscvi_embed.pkl"
)

mde_obj = run_mde(pooled_ad1.obsm["X_scvi"])

# conn, dist = mde_neighbors(mde_obj, pooled_ad1.n_obs)

# ad1.obsp["mde_connectivities"] = conn
# ad1.obsp["mde_distances"] = dist

# ad1.uns["neighbors_mde"] = {
#     "method": "mde",
#     "params": {"n_neighbors": 15},
#     "random_state": 0,
#     "use_rep": "X_scvi",
#     "connectivities_key": "mde_connectivities",
#     "distances_key": "mde_distances",
# }

pooled_ad1.obsm["X_mde"] = mde_embed(mde_obj)

sc.pp.neighbors(pooled_ad1, use_rep="X_scvi")

sc.tl.umap(pooled_ad1, min_dist=0.3)

aux_ad = anndata.AnnData(
    X=pooled_ad1.raw.X,
    obs=pooled_ad1.obs,
    var=pooled_ad1.raw.var,
    obsm=pooled_ad1.obsm,
    obsp=pooled_ad1.obsp,
    uns=pooled_ad1.uns,
)
del pooled_ad1

aux_ad.write(
    "scautoqc_integrated.h5ad",
    compression="gzip",
)
