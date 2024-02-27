#!/usr/bin/env python
#"""Gather starsolo and cellbender outputs into a single h5ad
#
#Usage: gather_matrices.py [options] 
#
#Options:
#  --debug            print debug information
#  --profile          print profile information
#  --dry              print starsolo and cellbender outputs to gather
#  --force            force overwrite if output exists
#  --cr_gene   path to 10x matrix filtered folder
#  --cr_velo    path to velocyto filtered folder
#  --cb_h5    path to cellbender h5
#"""


import logging
import signal
import sys
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os.path
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import sctk as sk

import h5py

def read_cellbender(
    input_h5,
    remove_zero=True,
    remove_nan=True,
    train_history=False,
    latent_gene_encoding=False,
    add_suffix=None,
):
    """
    Read cellbender output h5 generated from mtx input
    """
    import scipy.sparse as sp

    f = h5py.File(input_h5, "r")
    if "matrix" in f:
        mat = f["matrix"]
        feat = mat["features"]
        feat_name = feat["name"][()]
        vardict = {
            "gene_ids": feat["id"][()].astype(str),
            "feature_type": feat["feature_type"][:].astype(str),
        }
    elif "background_removed" in f:
        mat = f["background_removed"]
        vardict = {
            "gene_ids": mat["genes"][()].astype(str),
        }
        feat_name = mat["gene_names"][()]
    else:
        raise ValueError("The data doesn't look like cellbender output")
    n_var, n_obs = tuple(mat["shape"][()])
    if "metadata" in f:
        ad = read_cellbender_v3(f,
            feat=feat,
            feat_name=feat_name,
            vardict=vardict,
            n_var=n_var,
            n_obs=n_obs,
            remove_zero=remove_zero,
            remove_nan=remove_nan,
            train_history=train_history,
            latent_gene_encoding=latent_gene_encoding,
            add_suffix=add_suffix,)
    else:
        ad = sk.read_cellbender(input_h5)
    return ad

def read_cellbender_v3(f,
    feat,
    feat_name,
    vardict,
    n_var,
    n_obs,
    remove_zero=True,
    remove_nan=True,
    train_history=False,
    latent_gene_encoding=False,
    add_suffix=None,
):
    import numpy as np
    import anndata
    import scipy.sparse as sp
    import pandas as pd
    cols = ["cell_probability", "droplet_efficiency"]
    if f['droplet_latents']['barcode_indices_for_latents'].shape[0] < n_obs:
        bidx = f["droplet_latents"]["barcode_indices_for_latents"][()]
        obsdict = {}
        for x in cols:
            val = np.empty(n_obs)
            val.fill(np.nan)
            val[bidx] = f["droplet_latents"][x][()]
            obsdict[x] = val
        if latent_gene_encoding:
            lge = f["droplet_latents"]["gene_expression_encoding"][()]
            obsm = np.empty((n_obs, lge.shape[1]))
            obsm.fill(np.nan)
            obsm[bidx, :] = lge
    else:
        obsdict = {x: f["droplet_latents"][x] for x in cols}
        if latent_gene_encoding:
            obsm = f["droplet_latents"]["gene_expression_encoding"][()]
    barcodes = np.array(
        [b[:-2] if b.endswith("-1") else b for b in f["matrix"]["barcodes"][()].astype(str)]
    )
    ad = anndata.AnnData(
        X=sp.csr_matrix(
            (f["matrix"]["data"][()], f["matrix"]["indices"][()], f["matrix"]["indptr"][()]),
            shape=(n_obs, n_var),
        ),
        var=pd.DataFrame(vardict, index=feat_name.astype(str)),
        obs=pd.DataFrame(obsdict, index=barcodes),
        uns={
            "target_false_positive_rate": f['metadata']["target_false_positive_rate"][()],
            "test_elbo": list(f['metadata']['learning_curve_test_elbo']),
            "test_epoch": list(f['metadata']['learning_curve_test_epoch']),
            # "overall_change_in_train_elbo": list(f['metadata']["overall_change_in_train_elbo"]), this doesn't exist in default h5 v3
        }
        if train_history
        else {},
    )
    ad.var_names_make_unique()
    if latent_gene_encoding:
        ad.obsm["X_latent_gene_encoding"] = obsm

    mask_nan = np.isnan(ad.obs.cell_probability)
    mask_0 = ad.X.sum(axis=1).A1 <= 0

    mask_remove = np.zeros(n_obs).astype(bool)
    if remove_nan:
        mask_remove = mask_remove | mask_nan
    if remove_zero:
        mask_remove = mask_remove | mask_0
    idx_remove = np.where(mask_remove)[0]

    idx_sort = pd.Series(np.argsort(ad.obs_names))
    idx_sort = idx_sort[~idx_sort.isin(idx_remove)]

    ad1 = ad[idx_sort.values].copy()
    del ad

    if add_suffix:
        ad1.obs_names = ad1.obs_names.astype(str) + add_suffix

    return ad1

def gather_matrices(cr_gene_filtered_mtx, cr_velo_filtered_mtx, cb_filtered_h5):
    cr_gene_filtered_ad = sc.read_10x_mtx(cr_gene_filtered_mtx)
    logging.info("cr_gene_filtered_mtx done")
    cr_velo_filtered_ad = sk.read_velocyto(os.path.realpath(cr_velo_filtered_mtx))
    logging.info("cr_velo_filtered_mtx done")
    cb_gene_filtered_ad = read_cellbender(cb_filtered_h5)
    logging.info("cb_filtered_h5 done")

    common_cells = list(
        set(cr_gene_filtered_ad.obs_names.tolist())
        & set(cb_gene_filtered_ad.obs_names.tolist())
    )
    k_cr = cr_gene_filtered_ad.obs_names.isin(common_cells)
    k_cb = cb_gene_filtered_ad.obs_names.isin(common_cells)
    ad = anndata.AnnData(
        X=cb_gene_filtered_ad.X[np.where(k_cb)[0], :],
        obs=cb_gene_filtered_ad.obs[k_cb].copy(),
        var=cb_gene_filtered_ad.var.copy(),
        layers={
            "raw": cr_gene_filtered_ad.X[np.where(k_cr)[0], :],
            "spliced": cr_velo_filtered_ad.X[np.where(k_cr)[0]],
            "unspliced": cr_velo_filtered_ad.layers["unspliced"][np.where(k_cr)[0]],
            "ambiguous": cr_velo_filtered_ad.layers["ambiguous"][np.where(k_cr)[0]],
        }
    )
    for layer in ("spliced", "unspliced", "ambiguous"):
        ad.layers[layer].eliminate_zeros()
    return ad


def main(args):
    logging.debug(args)

    cr_gene_filtered_mtx = args.cr_gene
    cr_velo_filtered_mtx = args.cr_velo
    cb_filtered_h5 = args.cb_h5
    output_h5ad = "gene_velo_cellbender.filtered.h5ad"

    if not os.path.exists(cr_gene_filtered_mtx):
        raise FileNotFoundError(cr_gene_filtered_mtx)
    if not os.path.exists(cr_velo_filtered_mtx):
        raise FileNotFoundError(cr_velo_filtered_mtx)
    if not os.path.exists(cb_filtered_h5):
        raise FileNotFoundError(cb_filtered_h5)

    if args.dry:
        print(cr_gene_filtered_mtx, cr_velo_filtered_mtx, cb_filtered_h5)
    else:
        if os.path.exists(output_h5ad) and not args.force:
            logging.info(output_h5ad + " exists, skip gathering without --force")
        else:
            ad = gather_matrices(cr_gene_filtered_mtx, cr_velo_filtered_mtx, cb_filtered_h5)
            ad.write(output_h5ad, compression="gzip")
            
            logging.info("done")

    return 0


if __name__ == '__main__':
    import argparse
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument("--debug", default=None, help="print debug information")
    my_parser.add_argument("--profile", default=None, help="print profile information")
    my_parser.add_argument("--dry", default=None, help="print starsolo and cellbender outputs to gather")
    my_parser.add_argument("--force", default=None, help="force overwrite if output exists")
    my_parser.add_argument("--cr_gene", default=None, help="path to 10x matrix filtered folder")
    my_parser.add_argument("--cr_velo", default=None, help="path to velocyto filtered folder")
    my_parser.add_argument("--cb_h5", default=None, help="path to cellbender h5")
    args = my_parser.parse_args()
    try:
        if args.debug:
            logLevel = logging.DEBUG
        else:
            logLevel = logging.INFO
        logging.basicConfig(
            level=logLevel,
            format='%(asctime)s; %(levelname)s; %(funcName)s; %(message)s',
            datefmt='%y-%m-%d %H:%M:%S')
        if args.profile:
            import cProfile
            cProfile.run('main(args)')
        else:
            main(args)
    except KeyboardInterrupt:
        logging.warning('Interrupted')
        sys.exit(1)
