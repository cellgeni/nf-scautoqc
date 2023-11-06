#!/usr/bin/env python
"""Flag doublet by scrublet

Usage: program [options] --samp <sampleid> --input <input>

Options:
  --debug           print debug information
  --profile         print profile information
  --filter <str>    only consider cells that have "True" values in this variable of obs
  --groupby <str>   run scrublet for each group specified by this variable of obs
  --samp            sample id
  --input <input>   input h5ad
"""


import logging
import signal
import sys
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import gc
import pandas as pd
import scanpy as sc
import sctk as sk

import argparse

def _find_rep(adata, use_rep, return_var_names=True):
    if isinstance(use_rep, str):
        if use_rep == "raw" and adata.raw:
            X = adata.raw.X
            var_names = adata.raw.var_names.values
        elif use_rep == "X":
            X = adata.X
            var_names = adata.var_names.values
        elif use_rep in adata.layers.keys():
            X = adata.layers[use_rep]
            var_names = adata.var_names.values
        elif use_rep in adata.obsm.keys():
            X = adata.obsm[use_rep]
            var_names = np.array([f"{use_rep}{i+1}".replace("X_", "") for i in range(X.shape[1])])
        elif use_rep in adata.obs.keys():
            X = adata.obs[use_rep].values.reshape((adata.n_obs, 1))
            var_names = np.array([use_rep])
        else:
            raise ValueError("Invalid `use_rep` provided.")
    elif isinstance(use_rep, np.ndarray) and use_rep.shape[0] == adata.n_obs:
        x = use_rep
        var_names = np.arange(x.shape[1])
    elif isinstance(use_rep, pd.DataFrame) and use_rep.shape[0] == adata.n_obs:
        x = use_rep.values
        var_names = use_rep.columns.values
    else:
        raise ValueError("Invalid `use_rep` provided.")
    if return_var_names:
        return (X, var_names)
    return X

def main(args):
    logging.debug(args)

    input_h5ad = args.input
    filter_var = args.filter
    groupby = args.groupby

    ad = sc.read(input_h5ad)
    if (ad.X.sum(0) > 0).sum() < ad.shape[1] * 0.2:
        exit()
    if groupby and groupby in ad.obs.columns:
        groups = ad.obs[groupby].cat.categories
        dfs = []
        for i, grp in enumerate(groups):
            print(grp)
            k_grp = ad.obs[groupby] == grp
            if filter_var and filter_var in ad.obs.columns:
                filter_value = ad.obs[filter_var]
                if filter_value.dtype.kind == "b":
                    k_filter = filter_value.values
                elif filter_value.dtype.kind == "O":
                    k_filter = filter_value.values == "True"
                else:
                    raise TypeError(f"{filter_var} is not a boolean or categorical")
                k = k_grp & k_filter
            else:
                k = k_grp
            if k.sum() > 10:
                ad1 = ad[k]
                try:
                    grp_df = sk.run_scrublet(ad1, inplace=False)
                    grp_df.to_csv(output_csv.replace(".csv", f".{grp}.csv"))
                    dfs.append(grp_df)
                except:
                    logging.warn(f"{grp}: scrublet failed")
                del ad1
                gc.collect()
        df = pd.concat(dfs)
    else:
        if filter_var and filter_var in ad.obs.columns:
            filter_value = ad.obs[filter_var]
            if filter_value.dtype.kind == "b":
                k = filter_value.values
            elif filter_value.dtype.kind == "O":
                k = filter_value.values == "True"
            else:
                raise TypeError(f"{filter_var} is not a boolean or categorical")
            ad1 = ad[k]
        else:
            ad1 = ad

        df = sk.run_scrublet(ad1, inplace=False)

    df.index.name = args.samp
    df.to_csv("gene_velo_cellbender.good_qc_cluster_mito80.scrublet.csv")

    return 0


if __name__ == '__main__':
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument("--input", default=None, help="input h5ad")
    my_parser.add_argument("--samp", default=None, help="input h5ad")
    my_parser.add_argument("--filter", default=None, help="only consider cells that have True values in this variable of obs")
    my_parser.add_argument("--groupby", default=None, help="run scrublet for each group specified by this variable of obs")
    my_parser.add_argument("--debug", default=None, help="print debug information")
    my_parser.add_argument("--profile", default=None, help="print profile information")

    args = my_parser.parse_args()

    try:
        if args.debug:
            logLevel = logging.DEBUG
        else:
            logLevel = logging.WARN
        logging.basicConfig(
            level=logLevel,
            format="%(asctime)s; %(levelname)s; %(funcName)s; %(message)s",
            datefmt="%y-%m-%d %H:%M:%S"        )
        if args.profile:
            import cProfile
            cProfile.run("main(args)")
        else:
            main(args)
    except KeyboardInterrupt:
        logging.warning("Interrupted")
        sys.exit(1)
