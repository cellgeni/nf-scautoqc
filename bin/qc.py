#!/usr/bin/env python
"""Perform automatic QC

Usage: qc.py [options] <sample_id>

Options:
  --debug             print debug information
  --profile           print profile information
  --qc_metrics <str>  comma-separated list of QC metrics [default: log1p_n_counts,log1p_n_genes,percent_mito,percent_ribo,percent_hb,percent_top50,percent_soup,percent_spliced]
  --models <str>      comma-separated <name>:<model.pkl> pairs giving celltypist models to use [default: ctp_pred:Immune_All_High.pkl]
  --clst_res <float>  resolution for QC clustering [default: 0.2]
  --min_frac <float>  min frac of pass_auto_filter for a cluster to be called good [default: 0.5]
  --out_path <path>   path of the output files
  --sample_id         sample id
  --mode <str>        mode of qc operation, automatic-qc (default) or filtering
"""


import logging
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42

import celltypist
import scanpy as sc
import sctk as sk


celltypist_models = {
    "cecilia22_predH": "Immune_All_High.pkl",
    "cecilia22_predL": "Immune_All_Low.pkl",
    "elmentaite21_pred": "Cells_Intestinal_Tract.pkl",
    "suo22_pred": "Pan_Fetal_Human.pkl",
}


def calculate_qc(ad, run_scrublet=True):
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


def prepare_metric_pairs(qc_metrics):
    metric_pairs = [("log1p_n_counts", qm) for qm in qc_metrics[1:]]
    if "percent_mito" in qc_metrics and "percent_ribo" in qc_metrics:
        metric_pairs.append(("log(percent_mito)", "percent_ribo"))
    if 'spliced' in ad.layers and 'unspliced' in ad.layers:
        if "percent_soup" in qc_metrics and "percent_spliced" in qc_metrics:
            metric_pairs.append(("log(percent_soup)", "percent_spliced"))
        if "percent_spliced" in qc_metrics and "scrublet_score" in qc_metrics:
            metric_pairs.append(("percent_spliced", "scrublet_score"))
    return metric_pairs


def run_mito_qc_loop(ad, qc_metrics, metrics_csv, res, threshold):
    mito_thresholds = [20, 50, 80]

    for max_mito in mito_thresholds:
        if metrics_csv is not None:
            metrics = metrics_csv.copy()
            metrics.at['percent_mito', 'max'] = max_mito
        else:
            if 'spliced' in ad.layers and 'unspliced' in ad.layers:
                metrics = {
                    "n_counts": (1000, None, "log", "min_only", 0.1),
                    "n_genes": (100, None, "log", "min_only", 0.1),
                    "percent_mito": (0.1, max_mito, "log", "max_only", 0.1),
                    "percent_spliced": (50, 97.5, "log", "both", 0.1),
                }
            else:
                metrics = {
                    "n_counts": (1000, None, "log", "min_only", 0.1),
                    "n_genes": (100, None, "log", "min_only", 0.1),
                    "percent_mito": (0.1, max_mito, "log", "max_only", 0.1),
                }
        sk._pipeline.cellwise_qc(
            ad,
            metrics,
            cell_qc_key=f"good_qc_cell_mito{max_mito}"
        )
        sk._pipeline.multi_resolution_cluster_qc(
            ad,
            metrics=metrics,
            cell_qc_key=f"good_qc_cell_mito{max_mito}",
            key_added=f"good_qc_cluster_mito{max_mito}",
            consensus_call_key=f"consensus_passed_qc_mito{max_mito}",
            consensus_frac_key=f"consensus_fraction_mito{max_mito}",
        )
        ad.obs[f"pass_auto_filter_mito{max_mito}"] = ad.obs[f"good_qc_cluster_mito{max_mito}"]

    ad.obs["pass_auto_filter"] = 100
    for max_mito in mito_thresholds[::-1]:
        ad.obs.loc[
            ad.obs[f"pass_auto_filter_mito{max_mito}"], "pass_auto_filter"
        ] = max_mito

    ad.obs["good_qc_cluster"] = 100
    for max_mito in mito_thresholds[::-1]:
        ad.obs.loc[
            ad.obs[f"good_qc_cluster_mito{max_mito}"], "good_qc_cluster"
        ] = max_mito


def generate_qc_plots(ad, qc_metrics, metric_pairs, models, ctp_name):
    ad.uns["qc_cluster_colors"] = sk._plot.make_palette(
        ad.obs["qc_cluster"].cat.categories.size
    )
    sk.set_figsize((3, 3))
    metric_ufig = sc.pl.embedding(
        ad,
        basis="umap_qc",
        color=qc_metrics + ["good_qc_cluster", "qc_cluster"],
        ncols=int(np.ceil((len(qc_metrics) + 2) / 2)),
        show=False,
        vmax=[
            80
            if qm in ("percent_mito", "percent_hb", "percent_soup")
            else 0.5
            if qm == "scrublet_score"
            else None
            for qm in qc_metrics
        ]
        + [None, None],
        size=10,
        return_fig=True,
    )

    metric_vfig = sk.plot_qc_violin(
        ad, groupby="qc_cluster", metrics=qc_metrics, figsize=(3, 3), return_fig=True
    )

    good_cluster_sfig = sk.plot_qc_scatter(
        ad,
        color_by="good_qc_cluster",
        metric_pairs=metric_pairs,
        use_hexbin=False,
        figsize=(3.5, 3),
        wspace=0.3,
        vmin=0,
        vmax=100,
        s=10 ** (2 - np.round(np.log10(ad.n_obs))),
        return_fig=True,
    )
    good_cluster_sfig.suptitle("good qc cluster")

    ad.obs["pass_default"] = (
        (ad.obs.n_counts >= 1e3) & (ad.obs.n_genes >= 200) & (ad.obs.percent_mito < 20)
    )
    ad.obs["pass_default"] = ad.obs["pass_default"].astype("category")

    cell_size = 10 ** (2 - np.round(np.log10(ad.n_obs)))

    pass_default_sfig = sk.plot_qc_scatter(
        ad,
        metric_pairs=metric_pairs,
        color_by="pass_default",
        s=cell_size,
        figsize=(3.5, 3),
        return_fig=True,
    )
    pass_default_sfig.suptitle("pass default filter")

    pass_auto_sfig = sk.plot_qc_scatter(
        ad,
        metric_pairs=metric_pairs,
        color_by="pass_auto_filter",
        s=cell_size,
        figsize=(3.5, 3),
        vmin=0,
        vmax=100,
        return_fig=True,
    )
    pass_auto_sfig.suptitle("pass auto filter")

    ctp_prob_sfig = None
    ctp_pred_ufig = None
    if models is not None:
        ctp_prob_sfig = sk.plot_qc_scatter(
            ad,
            metric_pairs=metric_pairs,
            color_by=f"{ctp_name}_prob",
            s=cell_size,
            figsize=(3.5, 3),
            return_fig=True,
        )
        ctp_prob_sfig.suptitle(f"{ctp_name} prob")

        sk.set_figsize((3, 3))
        sc.pl.embedding(
            ad,
            basis="umap_qc",
            color=[f"{ctp_name}_prob", f"{ctp_name}_uncertain", ctp_name],
            size=30,
            legend_fontsize=12,
            wspace=0.4,
            show=False,
        )
        ctp_pred_ufig = plt.gcf()

    sc.pl.embedding(
        ad,
        basis="umap_qc",
        color=["good_qc_cluster", "pass_auto_filter", "pass_default", "qc_cluster"],
        size=30,
        legend_fontsize=12,
        show=False,
    )
    qc_cluster_ufig = plt.gcf()

    return (
        metric_vfig,
        metric_ufig,
        good_cluster_sfig,
        pass_default_sfig,
        pass_auto_sfig,
        ctp_prob_sfig,
        ctp_pred_ufig,
        qc_cluster_ufig,
    )


def run_QC(
    ad,
    qc_metrics=None,
    models=None,
    metrics_csv=None,
    res=0.2,
    threshold=0.5,
):
    if qc_metrics is None:
        qc_metrics = [
            "log1p_n_counts",
            "log1p_n_genes",
            "percent_mito",
            "percent_ribo",
            "percent_hb",
            "percent_top50",
            "percent_soup",
            "percent_spliced",
        ]

    metric_pairs = prepare_metric_pairs(qc_metrics)

    if models is not None:
        ctp_name = list(models.keys())[0]
        for name, mod in models.items():
            run_celltypist(ad, mod, min_prob=0.3, key_added=name)
    else:
        ctp_name = None

    calculate_qc(ad, run_scrublet=("scrublet_score" in qc_metrics))

    run_mito_qc_loop(ad, qc_metrics, metrics_csv, res, threshold)

    return generate_qc_plots(ad, qc_metrics, metric_pairs, models, ctp_name)


def run_celltypist(ad, model, min_prob=0.3, key_added="ctp_pred"):
    aux_ad = ad.copy()
    sc.pp.normalize_total(aux_ad, target_sum=1e4)
    sc.pp.log1p(aux_ad)
    try:
        pred = celltypist.annotate(aux_ad, model=model, majority_voting=True)
        ad.obs[key_added] = pred.predicted_labels.majority_voting
    except ValueError:
        pred = celltypist.annotate(aux_ad, model=model, majority_voting=False)
        ad.obs[key_added] = pred.predicted_labels.predicted_labels
    ad.obs[f"{key_added}_prob"] = pred.probability_matrix.max(axis=1)
    ad.obs[f"{key_added}_uncertain"] = ad.obs[key_added].astype(str)
    ad.obs.loc[
        ad.obs[f"{key_added}_prob"] < min_prob, f"{key_added}_uncertain"
    ] = "Uncertain"
    ad.obs[f"{key_added}_uncertain"] = ad.obs[f"{key_added}_uncertain"].astype(
        "category"
    )
    del aux_ad


def parse_model_option(model_str):
    models = {}
    for model in model_str.split(","):
        if ":" in model:
            name, mod = model.split(":")
        else:
            name = "ctp_pred"
            mod = model
        models[name] = mod

    return models


def process_sample(ad, qc_metrics, models, metrics_csv, clst_res, min_frac):

    qc_figs = run_QC(
        ad,
        qc_metrics=qc_metrics,
        models=models,
        metrics_csv=metrics_csv,
        res=clst_res,
        threshold=min_frac,
    )
    return qc_figs


def main(args):
    logging.debug(args)

    sid = args.sample_id
    if args.models is not None:
        models = parse_model_option(args.models)
    else:
        models = args.models
    clst_res = float(args.clst_res)
    min_frac = float(args.min_frac)
    metrics_csv = pd.read_csv(args.metrics_csv, index_col=0)

    input_h5ad = args.out_path 
    
    ad = sc.read(input_h5ad)

    if args.qc_metrics is not None:
        qc_metrics = args.qc_metrics.split(",")
    else:
        if 'spliced' in ad.layers and 'unspliced' in ad.layers:
            qc_metrics = None
        else:
            qc_metrics = [
                    "log1p_n_counts",
                    "log1p_n_genes",
                    "percent_mito",
                    "percent_ribo",
                    "percent_hb",
                    "percent_top50",
                    "percent_soup",
                ]

    (
        metric_vfig,
        metric_ufig,
        good_cluster_sfig,
        pass_default_sfig,
        pass_auto_sfig,
        ctp_prob_sfig,
        ctp_pred_ufig,
        qc_cluster_ufig,
    ) = process_sample(ad, qc_metrics, models, metrics_csv, clst_res, min_frac)

    metric_vfig.savefig(
        f"{sid}.qc_plot.metric_vfig.png", bbox_inches="tight"
    )
    metric_ufig.savefig(
        f"{sid}.qc_plot.metric_ufig.png", bbox_inches="tight"
    )
    good_cluster_sfig.savefig(
        f"{sid}.qc_plot.good_sfig.png", bbox_inches="tight"
    )
    pass_default_sfig.savefig(
        f"{sid}.qc_plot.default_sfig.png", bbox_inches="tight"
    )
    pass_auto_sfig.savefig(
        f"{sid}.qc_plot.auto_sfig.png", bbox_inches="tight"
    )
    if models is not None:
        ctp_prob_sfig.savefig(
            f"{sid}.qc_plot.ctp_sfig.png", bbox_inches="tight"
        )
        ctp_pred_ufig.savefig(
            f"{sid}.qc_plot.ctp_ufig.png", bbox_inches="tight"
        )
    qc_cluster_ufig.savefig(
        f"{sid}.qc_plot.cluster_ufig.png", bbox_inches="tight"
    )

    ad.obs['sampleID'] = sid
    
    ad.uns['scautoqc_ranges'] = ad.uns['scautoqc_ranges'].applymap(lambda x: x.item() if hasattr(x, "item") else x)

    ad.write(f"{sid}_postqc.h5ad", compression="gzip")

    if ad.obs['good_qc_cluster_mito80'].mean() < 0.25 or (ad.X.sum(0) > 0).sum() < ad.shape[1] * 0.2:
        open(f'{sid}_no-scr', 'a').close()
    else:
        open(f'{sid}_yes-scr', 'a').close()

    return 0


if __name__ == "__main__":
    import argparse

    def nullable_string(val):
        if not val:
            return None
        return val

    my_parser = argparse.ArgumentParser()
    my_parser.add_argument("--debug", default=None, help="print debug information")
    my_parser.add_argument("--profile", default=None, help="print profile information")
    my_parser.add_argument("--qc_metrics", default=None, help="comma-separated list of QC metrics [default: log1p_n_counts,log1p_n_genes,percent_mito,percent_ribo,percent_hb,percent_top50,percent_soup,percent_spliced]")
    my_parser.add_argument("--metrics_csv", type=nullable_string, nargs='?', help="csv file of metric cutoffs")
    my_parser.add_argument("--models", default=None, help="comma-separated <name>:<model.pkl> pairs giving celltypist models to use [default: ctp_pred:Immune_All_High.pkl]")
    my_parser.add_argument("--clst_res", default=None, help="resolution for QC clustering [default: 0.2]")
    my_parser.add_argument("--min_frac", default=None, help="min frac of pass_auto_filter for a cluster to be called good [default: 0.5]")
    my_parser.add_argument("--out_path", default=None, help="path of the output files")
    my_parser.add_argument("--sample_id", default=None, help="sample id")

    args = my_parser.parse_args()

    # args = docopt(__doc__)
    # args = {k.lstrip("-<").rstrip(">"): args[k] for k in args}
    try:
        if args.debug:
            logLevel = logging.DEBUG
        else:
            logLevel = logging.WARN
        logging.basicConfig(
            level=logLevel,
            format="%(asctime)s; %(levelname)s; %(funcName)s; %(message)s",
            datefmt="%y-%m-%d %H:%M:%S",
        )
        if args.profile:
            import cProfile

            cProfile.run("main(args)")
        else:
            main(args)
    except KeyboardInterrupt:
        logging.warning("Interrupted")
        sys.exit(1)
