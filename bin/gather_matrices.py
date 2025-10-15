#!/usr/bin/env python
#"""Gather STARsolo (or Cell Ranger), Velocyto, CellBender  outputs into a single h5ad
#
#Usage: gather_matrices.py [options] 
#
#Options:
#  --cr_gene   path to 10x matrix filtered folder
#  --cr_velo    path to velocyto filtered folder
#  --cb_h5    path to cellbender h5#"""


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


def gather_matrices(cr_gene_filtered_mtx, cr_velo_filtered_mtx, cb_filtered_h5, gather_mode):
    barcodes_file = os.path.join(cr_gene_filtered_mtx, "barcodes.tsv.gz")
    features_file = os.path.join(cr_gene_filtered_mtx, "features.tsv.gz")
    matrix_file = os.path.join(cr_gene_filtered_mtx, "matrix.mtx.gz")

    if os.path.exists(barcodes_file) and os.path.exists(features_file) and os.path.exists(matrix_file):
        input_mode = 'starsolo'
        if gather_mode == 'cellbender':
            cr_gene_filtered_mtx = os.path.realpath(cr_gene_filtered_mtx)
            cr_gene_filtered_mtx = os.path.dirname(cr_gene_filtered_mtx) + "/raw"
        cr_gene_filtered_ad = sc.read_10x_mtx(cr_gene_filtered_mtx)
    else: # if cellranger
        input_mode = 'cellranger'
        if gather_mode == 'cellbender':
            cr_gene_filtered_mtx = os.path.realpath(cr_gene_filtered_mtx)
            cr_gene_filtered_mtx = os.path.dirname(cr_gene_filtered_mtx) + "/raw_feature_bc_matrix.h5"
        cr_gene_filtered_ad = sc.read_10x_h5(cr_gene_filtered_mtx)
        cr_gene_filtered_ad.obs_names = cr_gene_filtered_ad.obs_names.str.replace("-1$", "", regex=True)
    logging.info(f"cr_gene_filtered_mtx done (input: {input_mode}, gather_mode: {gather_mode})")

    if cr_velo_filtered_mtx is not None:
        # If Velocyto output exists, process it and merge into h5ad
        cr_velo_filtered_ad = sk.read_velocyto(os.path.realpath(cr_velo_filtered_mtx))[cr_gene_filtered_ad.obs_names].copy()
        logging.info("cr_velo_filtered_mtx done")

    if cb_filtered_h5 is not None:
        # If CellBender output is provided, process it and merge into h5ad
        cb_gene_filtered_ad = sk.read_cellbender(cb_filtered_h5)
        logging.info("cb_filtered_h5 done")

        if args.gather_mode == 'starsolo':
            common_cells = list(
                set(cr_gene_filtered_ad.obs_names.tolist())
                & set(cb_gene_filtered_ad.obs_names.tolist())
            )
        elif args.gather_mode == 'cellbender':
            common_cells = list(
                set(cb_gene_filtered_ad.obs_names.tolist())
            )

        k_cr = cr_gene_filtered_ad.obs_names.isin(common_cells)
        k_cb = cb_gene_filtered_ad.obs_names.isin(common_cells)

        if cr_velo_filtered_mtx is not None:
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
        else:
            ad = anndata.AnnData(
                X=cb_gene_filtered_ad.X[np.where(k_cb)[0], :],
                obs=cb_gene_filtered_ad.obs[k_cb].copy(),
                var=cb_gene_filtered_ad.var.copy(),
                layers={
                    "raw": cr_gene_filtered_ad.X[np.where(k_cr)[0], :],
                }
            )
            return ad
        
    else:  # If CellBender output is not provided
        if cr_velo_filtered_mtx is not None:
            ad = anndata.AnnData(
                X=cr_gene_filtered_ad.X,
                obs=cr_gene_filtered_ad.obs.copy(),
                var=cr_gene_filtered_ad.var.copy(),
                layers={
                    "raw": cr_gene_filtered_ad.X,
                    "spliced": cr_velo_filtered_ad.X,
                    "unspliced": cr_velo_filtered_ad.layers["unspliced"],
                    "ambiguous": cr_velo_filtered_ad.layers["ambiguous"],
                }
            )
            for layer in ("spliced", "unspliced", "ambiguous"):
                ad.layers[layer].eliminate_zeros()
            return ad
        else:
            # If neither CellBender nor velocyto data is provided
            ad = anndata.AnnData(
                X=cr_gene_filtered_ad.X,
                obs=cr_gene_filtered_ad.obs.copy(),
                var=cr_gene_filtered_ad.var.copy(),
                layers={
                    "raw": cr_gene_filtered_ad.X,
                }
            )
            return ad

def main(args):
    logging.info(args)

    cr_gene_filtered_mtx = args.cr_gene
    cr_velo_filtered_mtx = args.cr_velo
    cb_filtered_h5 = args.cb_h5
    gather_mode = args.gather_mode
    output_h5ad = "gene_velo_cellbender.filtered.h5ad"

    if not os.path.exists(cr_gene_filtered_mtx):
        raise FileNotFoundError(cr_gene_filtered_mtx)

    ad = gather_matrices(cr_gene_filtered_mtx, cr_velo_filtered_mtx, cb_filtered_h5, gather_mode)
    if args.cell_or_nuclei is not None:
        ad.uns['cell_or_nuclei'] = args.cell_or_nuclei
    else:
        print("No cell or nuclei information provided, assuming 'cell'")
        ad.uns['cell_or_nuclei'] = 'cell'

    ad.write(output_h5ad, compression="gzip")

    logging.info("done")

    return 0

if __name__ == '__main__':
    import argparse
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument("--cr_gene", default=None, help="path to 10x matrix filtered folder")
    my_parser.add_argument("--cr_velo", default=None, const=None, nargs='?', help="path to velocyto filtered folder")
    my_parser.add_argument("--cb_h5", default=None, const=None, nargs='?', help="path to cellbender h5")
    my_parser.add_argument("--cell_or_nuclei", default=None, const=None, nargs='?', help="path to cellbender h5")
    my_parser.add_argument("--gather_mode", default=None, const=None, nargs='?', help="gathering mode: original (starsolo-focused), cellbender (cellbender-focused)")
    args = my_parser.parse_args()
    try:
        logLevel = logging.INFO
        logging.basicConfig(
            level=logLevel,
            format='%(asctime)s; %(levelname)s; %(funcName)s; %(message)s',
            datefmt='%y-%m-%d %H:%M:%S')
        main(args)
    except KeyboardInterrupt:
        logging.warning('Interrupted')
        sys.exit(1)
