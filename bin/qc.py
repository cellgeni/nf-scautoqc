#!/usr/bin/env python
"""Perform automatic QC.

This script runs single-cell QC and produces figures and an updated AnnData file.

Notes
-----
- The CLI is preserved for pipeline compatibility; new optional flags are additive.
- Figures and outputs follow the existing naming scheme based on sample_id.
"""
import logging
import signal
import sys
from pathlib import Path

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib_venn import venn3

rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42

import celltypist
import scanpy as sc
import sctk as sk

# -------------------------
# Constants and utilities
# -------------------------

MITO_THRESHOLDS = (20, 50, 80)
DEFAULT_QC_METRICS = (
    "log1p_n_counts",
    "log1p_n_genes",
    "percent_mito",
    "percent_ribo",
    "percent_hb",
    "percent_top50",
    "percent_soup",
    "percent_spliced",
)

def nullable_string(val):
    """Return None for empty/falsey CLI values, otherwise passthrough."""
    if not val:
        return None
    return val

def parse_model_option(model_str):
    """Parse comma-separated name:model pairs into a mapping.

    Examples: "ctp_pred:Immune_All_High.pkl,other:Foo.pkl".
    If a value lacks a name prefix, use 'ctp_pred'.
    """
    models = {}
    for model in filter(None, (m.strip() for m in model_str.split(","))):
        if ":" in model:
            name, mod = model.split(":", 1)
        else:
            name = "ctp_pred"
            mod = model
        models[name.strip()] = mod.strip()
    return models


def calculate_qc(ad, run_scrublet=True):
    """Calculate QC metrics and annotate the AnnData object in-place.

    Assumes:
    - .X contains post-soup-removal counts
    - ad.layers may include: 'raw' (pre-soup), 'spliced', 'unspliced'
    """
    sk.calculate_qc(ad, log1p=False)
    if "raw" in ad.layers:
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


def run_celltypist(ad, model, min_prob=0.3, key_added="ctp_pred"):
    """Run Celltypist on a normalized/log1p copy and add predictions to `ad.obs`."""
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
    ad.obs.loc[ad.obs[f"{key_added}_prob"] < min_prob, f"{key_added}_uncertain"] = "Uncertain"
    ad.obs[f"{key_added}_uncertain"] = ad.obs[f"{key_added}_uncertain"].astype("category")
    del aux_ad


def prepare_metric_pairs(ad, qc_metrics):
    """Return common metric pairs for scatter/hexbin QC plots."""
    metric_pairs = [("log1p_n_counts", qm) for qm in qc_metrics[1:]]
    if "percent_mito" in qc_metrics and "percent_ribo" in qc_metrics:
        metric_pairs.append(("log(percent_mito)", "percent_ribo"))
    if ('spliced' in ad.layers) and ('unspliced' in ad.layers):
        if "percent_soup" in qc_metrics and "percent_spliced" in qc_metrics:
            metric_pairs.append(("log(percent_soup)", "percent_spliced"))
        if "percent_spliced" in qc_metrics and "scrublet_score" in qc_metrics:
            metric_pairs.append(("percent_spliced", "scrublet_score"))
    return metric_pairs

def run_qc_mito_loop_original(ad, qc_metrics, metrics_custom, threshold):
    """Original QC approach looping over multiple mitochondrial thresholds."""

    sk._pipeline.generate_qc_clusters(ad, metrics=qc_metrics)

    for max_mito in MITO_THRESHOLDS:
        if metrics_custom is None:
            if ad.uns['cell_or_nuclei'] == 'cell':
                metrics = {
                    "n_counts": (1000, None, "log", "min_only", 0.1),
                    "n_genes": (100, None, "log", "min_only", 0.1),
                    "percent_mito": (0.1, max_mito, "log", "max_only", 0.1),
                    "percent_spliced": (50, 97.5, "log", "both", 0.1),
                }
            elif ad.uns['cell_or_nuclei'] == 'nuclei':
                metrics = {
                    "n_counts": (300, None, "log", "min_only", 0.1),
                    "n_genes": (100, None, "log", "min_only", 0.1),
                    "percent_mito": (0.1, max_mito, "log", "max_only", 0.1),
                    "percent_spliced": (0.01, 100, "log", "max_only", 0.1),
                }
        else:
            metrics = {k: tuple(v.values()) for k, v in metrics_custom.to_dict(orient="index").items()}
            metrics["percent_mito"] = (metrics["percent_mito"][0], max_mito, *metrics["percent_mito"][2:])
        if not (('spliced' in ad.layers) and ('unspliced' in ad.layers)):
            del metrics['percent_spliced']

        sk._pipeline.cellwise_qc(
            ad,
            metrics,
            cell_qc_key=f"good_qc_cell_mito{max_mito}"
        )
        sk._pipeline.clusterwise_qc(
            ad,
            cell_qc_key=f"good_qc_cell_mito{max_mito}",
            key_added=f"good_qc_cluster_mito{max_mito}",
            threshold=threshold
        )
        ad.obs[f"pass_auto_filter_mito{max_mito}"] = ad.obs[f"good_qc_cluster_mito{max_mito}"]

    ad.obs["pass_auto_filter"] = 100
    for max_mito in reversed(MITO_THRESHOLDS):
        ad.obs.loc[ad.obs[f"pass_auto_filter_mito{max_mito}"], "pass_auto_filter"] = max_mito

    ad.obs["good_qc_cluster"] = 100
    for max_mito in reversed(MITO_THRESHOLDS):
        ad.obs.loc[ad.obs[f"good_qc_cluster_mito{max_mito}"], "good_qc_cluster"] = max_mito

def run_qc_multi_res(ad, qc_metrics, metrics_custom, threshold):
    """Multi-resolution QC approach with consensus across clusterings."""

    sk._pipeline.generate_qc_clusters(ad, metrics=qc_metrics)

    if metrics_custom is None:
        if ad.uns['cell_or_nuclei'] == 'cell':
            metrics = {
                "n_counts": (1000, None, "log", "min_only", 0.1),
                "n_genes": (100, None, "log", "min_only", 0.1),
                "percent_mito": (0.1, 20, "log", "max_only", 0.1),
                "percent_spliced": (50, 97.5, "log", "both", 0.1),
            }
        elif ad.uns['cell_or_nuclei'] == 'nuclei':
            metrics = {
                "n_counts": (300, None, "log", "min_only", 0.1),
                "n_genes": (100, None, "log", "min_only", 0.1),
                "percent_mito": (0.1, 20, "log", "max_only", 0.1),
                "percent_spliced": (0.01, 100, "log", "max_only", 0.1),
            }
    else:
        metrics = {k: tuple(v.values()) for k, v in metrics_custom.to_dict(orient="index").items()}
    if not (('spliced' in ad.layers) and ('unspliced' in ad.layers)):
        del metrics['percent_spliced']

    sk._pipeline.cellwise_qc(ad, metrics)
    sk._pipeline.multi_resolution_cluster_qc(ad, metrics=qc_metrics)

    # Figure out consensus threshold that is the closest match to the single cell level calls
    best_jaccard = 0.0
    consensus_threshold = 0.0
    for this_threshold in np.unique(ad.obs["consensus_fraction"]):
        #threshold the calls
        ad.obs["consensus_passed_qc"] = (ad.obs["consensus_fraction"] >= this_threshold)
        #compute a jaccard of the current thresholding versus the cell level calls
        #negate to compute the jaccard of the cells failing qc
        this_jaccard = np.sum(~ad.obs["cell_passed_qc"] & ~ad.obs["consensus_passed_qc"])/np.sum(~ad.obs["cell_passed_qc"] | ~ad.obs["consensus_passed_qc"])
        if this_jaccard > best_jaccard:
            best_jaccard = this_jaccard
            consensus_threshold = this_threshold
    print("Best consensus overlap found for threshold "+str(consensus_threshold))
    ad.obs["consensus_passed_qc"] = ad.obs["consensus_fraction"] >= consensus_threshold

def run_qc_combined(ad, qc_metrics, metrics_custom, threshold):
    """Combined original + multi-resolution QC across mitochondrial thresholds."""

    for max_mito in MITO_THRESHOLDS:
        if metrics_custom is None:
            if ad.uns['cell_or_nuclei'] == 'cell':
                metrics = {
                    "n_counts": (1000, None, "log", "min_only", 0.1),
                    "n_genes": (100, None, "log", "min_only", 0.1),
                    "percent_mito": (0.1, max_mito, "log", "max_only", 0.1),
                    "percent_spliced": (50, 97.5, "log", "both", 0.1),
                }
            elif ad.uns['cell_or_nuclei'] == 'nuclei':
                metrics = {
                    "n_counts": (300, None, "log", "min_only", 0.1),
                    "n_genes": (100, None, "log", "min_only", 0.1),
                    "percent_mito": (0.1, max_mito, "log", "max_only", 0.1),
                    "percent_spliced": (0.01, 100, "log", "max_only", 0.1),
                }
        else:
            metrics = {k: tuple(v.values()) for k, v in metrics_custom.to_dict(orient="index").items()}
            metrics["percent_mito"] = (metrics["percent_mito"][0], max_mito, *metrics["percent_mito"][2:])
        if not (('spliced' in ad.layers) and ('unspliced' in ad.layers)):
            del metrics['percent_spliced']

        sk._pipeline.cellwise_qc(
            ad,
            metrics,
            cell_qc_key=f"good_qc_cell_mito{max_mito}"
        )
        sk._pipeline.multi_resolution_cluster_qc(
            ad,
            metrics=qc_metrics,
            cell_qc_key=f"good_qc_cell_mito{max_mito}",
            key_added=f"good_qc_cluster_mito{max_mito}",
            consensus_call_key=f"consensus_passed_qc_mito{max_mito}",
            consensus_frac_key=f"consensus_fraction_mito{max_mito}",
            threshold=threshold
        )
        ad.obs[f"pass_auto_filter_mito{max_mito}"] = ad.obs[f"good_qc_cluster_mito{max_mito}"]
    ad.obs["pass_auto_filter"] = 100
    for max_mito in reversed(MITO_THRESHOLDS):
        ad.obs.loc[ad.obs[f"pass_auto_filter_mito{max_mito}"], "pass_auto_filter"] = max_mito

    ad.obs["good_qc_cluster"] = 100
    for max_mito in reversed(MITO_THRESHOLDS):
        ad.obs.loc[ad.obs[f"good_qc_cluster_mito{max_mito}"], "good_qc_cluster"] = max_mito

    for max_mito in MITO_THRESHOLDS:
        # Figure out consensus threshold that is the closest match to the single cell level calls
        best_jaccard = 0.0
        consensus_threshold = 0.0
        for this_threshold in np.unique(ad.obs[f"consensus_fraction_mito{max_mito}"]):
            #threshold the calls
            ad.obs["consensus_passed_qc"] = (ad.obs[f"consensus_fraction_mito{max_mito}"] >= this_threshold)
            #compute a jaccard of the current thresholding versus the cell level calls
            #negate to compute the jaccard of the cells failing qc
            this_jaccard = np.sum(~ad.obs[f"good_qc_cell_mito{max_mito}"] & ~ad.obs["consensus_passed_qc"])/np.sum(~ad.obs[f"good_qc_cell_mito{max_mito}"] | ~ad.obs["consensus_passed_qc"])
            if this_jaccard > best_jaccard:
                best_jaccard = this_jaccard
                consensus_threshold = this_threshold
        print("Best consensus overlap found for threshold "+str(consensus_threshold))
        ad.obs[f"consensus_passed_qc_mito{max_mito}"] = (ad.obs[f"consensus_fraction_mito{max_mito}"] >= consensus_threshold)

    ad.obs["consensus_pass_auto_filter"] = 100  # sentinel
    for max_mito in reversed(MITO_THRESHOLDS):  # start with easy (80) → hard (20)
        ad.obs.loc[ad.obs[f"consensus_passed_qc_mito{max_mito}"], "consensus_pass_auto_filter"] = max_mito


def generate_qc_plots(ad, qc_metrics, qc_mode, metric_pairs, ctp_models, ctp_name):
    """Generate and return QC figures for later saving."""
    ad.uns["qc_cluster_colors"] = sk._plot.make_palette(
        ad.obs["qc_cluster"].cat.categories.size
    )
    sk.set_figsize((3, 3))
    if qc_mode == 'multires':
        good_qc_cluster = "consensus_passed_qc"
        additional_metrics = ['cluster_passed_qc', 'consensus_fraction', good_qc_cluster, 'qc_cluster']
    elif qc_mode == 'combined':
        good_qc_cluster = "consensus_pass_auto_filter"
        additional_metrics = [good_qc_cluster, 'qc_cluster']
    else:
        good_qc_cluster = "good_qc_cluster"
        additional_metrics = [good_qc_cluster, "qc_cluster"]
    metric_ufig = sc.pl.embedding(
        ad,
        basis="umap_qc",
        color=qc_metrics + additional_metrics,
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
        color_by=good_qc_cluster,
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
    if ctp_models is not None:
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

    if qc_mode == 'multires':  
        sc.pl.embedding(
            ad,
            basis="umap_qc",
            color=["consensus_passed_qc"],
            size=30,
            legend_fontsize=12,
            show=False,
        )
    else:
        sc.pl.embedding(
            ad,
            basis="umap_qc",
            color=["good_qc_cluster", "pass_auto_filter", "pass_default", "qc_cluster"],
            size=30,
            legend_fontsize=12,
            show=False,
        )
    qc_cluster_ufig = plt.gcf()

    qc_consensus_ufig = None
    qc_consensus_venn = None
    if qc_mode == 'multires':
        sc.pl.embedding(
            ad, 
            basis="umap_qc", 
            color=["cell_passed_qc", "cluster_passed_qc", "consensus_fraction", "consensus_passed_qc"],
            size=30,
            legend_fontsize=12,
            show=False,
        )
        qc_consensus_ufig = plt.gcf()

        cell_failed_qc = set(ad.obs_names[~ad.obs["cell_passed_qc"]])
        cluster_failed_qc = set(ad.obs_names[~ad.obs["cluster_passed_qc"]])
        consensus_failed_qc = set(ad.obs_names[~ad.obs["consensus_passed_qc"]])
        venn3([cell_failed_qc, cluster_failed_qc, consensus_failed_qc], set_labels=["cell", "cluster", "consensus"])
        qc_consensus_venn = plt.gcf()

    return (
        metric_vfig,
        metric_ufig,
        good_cluster_sfig,
        pass_default_sfig,
        pass_auto_sfig,
        ctp_prob_sfig,
        ctp_pred_ufig,
        qc_cluster_ufig,
        qc_consensus_ufig,
        qc_consensus_venn
    )
def run_qc(ad, ctp_models=None, qc_mode="original", metrics_custom=None, threshold=0.5):
    """End-to-end QC run: compute metrics, run QC, and build figures."""

    qc_metrics = list(DEFAULT_QC_METRICS)
    if not (('spliced' in ad.layers) and ('unspliced' in ad.layers)):
        # Drop percent_spliced if velocity layers are absent
        qc_metrics.pop()

    metric_pairs = prepare_metric_pairs(ad, qc_metrics)

    ctp_name = None
    if ctp_models is not None:
        logging.info(f"- Running Celltypist for {len(ctp_models.items())} models")
        ctp_name = list(ctp_models.keys())[-1]
        for name, mod in ctp_models.items():
            logging.info(f"- Running [{name}]")
            run_celltypist(ad, mod, min_prob=0.3, key_added=name)
    else:
        ctp_name = None

    logging.info("Calculating QC metrics")
    calculate_qc(ad, run_scrublet=("scrublet_score" in qc_metrics))

    logging.info(f"Starting QC: mode={qc_mode}")
    if qc_mode == "original":
        run_qc_mito_loop_original(ad, qc_metrics, metrics_custom, threshold)
    elif qc_mode == "multires":
        run_qc_multi_res(ad, qc_metrics, metrics_custom, threshold)
    elif qc_mode == "combined":
        run_qc_combined(ad, qc_metrics, metrics_custom, threshold)
    else:
        raise ValueError(f"Unknown QC mode: {qc_mode}")

    logging.info("Generating QC plots")
    return generate_qc_plots(ad, qc_metrics, qc_mode, metric_pairs, ctp_models, ctp_name)


def process_sample(ad, ctp_models, qc_mode, metrics_custom, min_frac):
    """Wrapper around run_qc to keep the original structure."""
    return run_qc(
        ad=ad,
        ctp_models=ctp_models,
        qc_mode=qc_mode,
        metrics_custom=metrics_custom,
        threshold=min_frac,
    )


def main(args):
    """Entry point for CLI execution."""
    logging.info(args)

    celltypist_models = {
        "gut": {
            "cecilia22_predH": "Immune_All_High.pkl",
            "cecilia22_predL": "Immune_All_Low.pkl",
            "elmentaite21_pred": "Cells_Intestinal_Tract.pkl",
            "suo22_pred": "Pan_Fetal_Human.pkl",
            "megagut_pred": "/nfs/cellgeni/tickets/tic-2456/actions/MegaGut_Human.pkl"
        }
    }

    sid = args.sample_id
    if args.celltypist_model and ":" in args.celltypist_model:
        logging.info(f"Using custom celltypist models: {args.celltypist_model}")
        ctp_models = parse_model_option(args.celltypist_model)
    elif args.celltypist_model:
        # Treat as predefined set key if provided
        ctp_models = celltypist_models.get(args.celltypist_model)
        if ctp_models is None:
            logging.info(f"Unknown predefined celltypist set '{args.celltypist_model}'; skipping celltypist")
    else:
        ctp_models = None
    
    if args.min_frac is not None:
        min_frac = float(args.min_frac)
    else:
        min_frac = float(0.5)

    if args.metrics_csv is not None:
        logging.info(f"Using custom metrics from {args.metrics_csv}")
        metrics_custom = pd.read_csv(Path(args.metrics_csv), index_col=0)
    else:
        logging.info("Using default metrics")
        metrics_custom = None
    
    logging.info(f"Loading AnnData object from {args.gath_obj}")

    ad = sc.read(Path(args.gath_obj))

    qc_mode = args.qc_mode

    logging.info("Running QC")
    (
        metric_vfig,
        metric_ufig,
        good_cluster_sfig,
        pass_default_sfig,
        pass_auto_sfig,
        ctp_prob_sfig,
        ctp_pred_ufig,
        qc_cluster_ufig,
        qc_consensus_ufig,
        qc_consensus_venn
    ) = process_sample(ad, ctp_models, qc_mode, metrics_custom, min_frac)

    logging.info("Saving QC plots")
    def _save(fig, filename):
        fig.savefig(filename, bbox_inches="tight")

    _save(metric_vfig, f"{sid}.qc_plot.metric_vfig.png")
    _save(metric_ufig, f"{sid}.qc_plot.metric_ufig.png")
    _save(good_cluster_sfig, f"{sid}.qc_plot.good_sfig.png")
    _save(pass_default_sfig, f"{sid}.qc_plot.default_sfig.png")
    _save(pass_auto_sfig, f"{sid}.qc_plot.auto_sfig.png")
    if ctp_models is not None and ctp_prob_sfig is not None and ctp_pred_ufig is not None:
        _save(ctp_prob_sfig, f"{sid}.qc_plot.ctp_sfig.png")
        _save(ctp_pred_ufig, f"{sid}.qc_plot.ctp_ufig.png")
    _save(qc_cluster_ufig, f"{sid}.qc_plot.cluster_ufig.png")

    if qc_mode == 'multires' and qc_consensus_ufig is not None and qc_consensus_venn is not None:
        _save(qc_consensus_ufig, f"{sid}.qc_plot.consensus_ufig.png")
        _save(qc_consensus_venn, f"{sid}.qc_plot.consensus_venn.png")

    ad.obs['sampleID'] = sid
    
    logging.info(f"Saving the final ranges to {sid}_metrics.csv")
    ad.uns['scautoqc_ranges'] = ad.uns['scautoqc_ranges'].applymap(lambda x: x.item() if hasattr(x, "item") else x)
    ad.uns['scautoqc_ranges'].index.name = f'metrics|{sid}'
    ad.uns['scautoqc_ranges'].to_csv(f"{sid}_metrics.csv")

    ad.uns['qc_mode'] = qc_mode

    logging.info(f"Saving AnnData object to {sid}_postqc.h5ad")
    ad.write(f"{sid}_postqc.h5ad", compression="gzip")

    logging.info("Checking whether this sample has good QC clusters in general")
    condition_qc = 'consensus_passed_qc' if qc_mode == 'multires' else 'good_qc_cluster_mito80'
    if (ad.obs[condition_qc].mean() < 0.25) or (ad.X.sum(0) > 0).sum() < ad.shape[1] * 0.2:
        Path(f'{sid}_no-scr').touch()
    else:
        Path(f'{sid}_yes-scr').touch()

    return 0


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run sc-auto QC on a gathered AnnData object.")
    parser.add_argument("--qc_mode", default="original", help="qc mode: original (megagut), multires, combined [default: original]")
    parser.add_argument("--metrics_csv", type=nullable_string, nargs='?', help="CSV file of metric cutoffs")
    parser.add_argument("--celltypist_model", default=None, help="comma-separated <name>:<model.pkl> pairs or a predefined set key (e.g., 'gut')")
    parser.add_argument("--min_frac", default=None, help="min frac of pass_auto_filter for a cluster to be called good [default: 0.5]")
    parser.add_argument("--gath_obj", default=None, help="path to the AnnData object from gather_matrices step")
    parser.add_argument("--sample_id", default=None, help="sample id")

    args = parser.parse_args()

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
