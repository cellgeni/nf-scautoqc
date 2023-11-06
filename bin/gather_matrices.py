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


def gather_matrices(cr_gene_filtered_mtx, cr_velo_filtered_mtx, cb_filtered_h5):
    cr_gene_filtered_ad = sc.read_10x_mtx(cr_gene_filtered_mtx)
    logging.info("cr_gene_filtered_mtx done")
    cr_velo_filtered_ad = sk.read_velocyto(os.path.realpath(cr_velo_filtered_mtx))
    logging.info("cr_velo_filtered_mtx done")
    cb_gene_filtered_ad = sk.read_cellbender(cb_filtered_h5)
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
