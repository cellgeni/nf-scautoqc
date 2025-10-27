#!/usr/bin/env python
"""Add metadata to pooled object

Usage: sample_qc.py --obj pooled_obj.h5ad 

Options:
  --obj        path of the pooled h5ad
  --batch      batch key
  --n_top_genes number of top genes to use
"""
    
import scanpy as sc
import pandas as pd
import scvi
import matplotlib.pyplot as plt
import anndata
import gc
import argparse
import joblib
# import rapids_singlecell as rsc

my_parser = argparse.ArgumentParser()
my_parser.add_argument("--obj", default=None, help="path of the pooled h5ad")
my_parser.add_argument("--batch", default=None, help="batch key")
my_parser.add_argument("--n_top_genes", default=7500, help="number of top genes to use")
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
    n_top_genes=int(args.n_top_genes),
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
)

joblib.dump(
    vae,
    "scvi_model.pkl",
)

plt.plot(vae.history["elbo_train"])[0].figure.savefig("elbo_training.png")

pooled_ad1.obsm["X_scvi"] = vae.get_latent_representation()

joblib.dump(
    pooled_ad1.obsm["X_scvi"],
    "Xscvi_embed.pkl"
)

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
