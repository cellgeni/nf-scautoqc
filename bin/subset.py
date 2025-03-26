#!/usr/bin/env python
"""Subset the object based on the cutoffs provided
Usage: subset.py --sample_id sample_id --cr_prefix /path/to/cellranger_outputs --limits_csv /path/to/cutoffs.csv

Options:
  --sample_id         sample id
  --cr_prefix         prefix for cell ranger outputs
  --limits_csv        csv file with subset cutoffs
"""


import logging
import signal
import sys
import os
from pathlib import Path

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import pandas as pd
import scanpy as sc

import sctk as sk


def calculate_qc(ad, run_scrublet=False):
    """Calculate QC metrics assuming:
    1) .X contains post-soup-removal counts
    2) .raw contains pre-soup-removal counts
    3) .spliced contains spliced counts
    4) .unspliced contains unspliced counts
    """
    sk.calculate_qc(ad, log1p=False)
    if 'raw' in ad.layers:
        sk.calculate_qc(ad, suffix="_raw", layer="raw", log1p=False)
        ad.obs["percent_soup"] = (1 - ad.obs["n_counts"] / ad.obs["n_counts_raw"]) * 100
    if 'spliced' in ad.layers and 'unspliced' in ad.layers:
        sk.calculate_qc(ad, suffix="_spliced", layer="spliced", log1p=False)
        sk.calculate_qc(ad, suffix="_unspliced", layer="unspliced", log1p=False)
        ad.obs["percent_spliced"] = (
            ad.obs["n_counts_spliced"]
            / (ad.obs["n_counts_spliced"] + ad.obs["n_counts_unspliced"])
            * 100
        )
    if run_scrublet:
        sk.run_scrublet(ad)


def subset_anndata(adata, cutoffs):
    for _, row in cutoffs.iterrows():
        metric, low, high = row['metric'], row['low'], row['high']
        if pd.notna(low):
            adata = adata[adata.obs[metric] >= low]
        if pd.notna(high):
            adata = adata[adata.obs[metric] <= high]
    return adata.copy()


def main(args):

    sid = args.sample_id
    root_dir = Path(args.cr_prefix) / sid

    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.startswith("filtered_feature_bc_matrix"):
                full_path = os.path.join(dirpath, filename)
                continue
    
    if full_path.endswith(".h5"):
        ad = sc.read_10x_h5(full_path)
    else:
        ad = sc.read(full_path)

    calculate_qc(ad)

    # filter out mouse cells and genes
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.startswith("gem_classification"):
                gem_class = pd.read_csv(os.path.join(dirpath, filename), index_col=0)
                ad.obs = pd.concat([ad.obs, gem_class], axis=1)
                ad = ad[ad.obs['call'] == 'GRCh38', ad.var['genome'] == 'GRCh38'].copy()
              
    # subset according to the cutoffs
    cutoffs = pd.read_csv(args.limits_csv)
    ad = subset_anndata(ad, cutoffs)

    ad.obs['sampleID'] = sid

    ad.var['gene_name'] = ad.var_names
    ad.var_names = ad.var['gene_ids']
    
    ad.write("subsetted.h5ad", compression="gzip")

    if (ad.X.sum(0) > 0).sum() < ad.shape[1] * 0.2:
        open(f'{sid}_no-scr', 'a').close()
    else:
        open(f'{sid}_yes-scr', 'a').close()

    return 0


if __name__ == "__main__":
    import argparse
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument("--sample_id", default=None, help="sample id")
    my_parser.add_argument("--cr_prefix", default=None, help="prefix for cell ranger outputs")
    my_parser.add_argument("--limits_csv", default=None, help="csv file with subset cutoffs")

    args = my_parser.parse_args()

    try:
        logLevel = logging.WARN
        logging.basicConfig(
            level=logLevel,
            format="%(asctime)s; %(levelname)s; %(funcName)s; %(message)s",
            datefmt="%y-%m-%d %H:%M:%S",
        )
        main(args)
    except KeyboardInterrupt:
        logging.warning("Interrupted")
        sys.exit(1)
