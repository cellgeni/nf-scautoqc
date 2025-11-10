# nf-scautoqc

nf-scautoqc is the Nextflow implementation of [scAutoQC pipeline](https://teichlab.github.io/sctk/notebooks/automatic_qc.html) used in [Oliver et al, 2024](https://doi.org/10.1038/s41586-024-07571-1). 

## How to run:

The recommended way to use nextflow is to run it in a screen session. These steps can be directly used in Sanger's FARM, but you can modify each step according to the environment you're working on or the job scheduler your HPC uses:

1. Start a screen session: `screen -S nf_run1`
2. Start a small interactive job for nextflow: `bsub -G cellgeni -n1 -R"span[hosts=1]" -Is -q long -R"select[mem>2000] rusage[mem=2000]" -M2000 bash`
3. Modify one of RESUME scripts in examples folder (pre-made Nextflow run scripts)
4. Run the RESUME scripts you modified: `./RESUME-scautoqc-all`
5. You can leave your screen session and let it run in the background: `Ctrl+A, D`

## Files:

* `main.nf` - the Nextflow pipeline that executes scAutoQC pipeline.
* `nextflow.config` - the configuration script that allows the processes to be submitted to IBM LSF on Sanger's HPC and ensures correct environment is set via singularity container (this is an absolute path). Global default parameters are also set in this file and some contain absolute paths.
* `examples/` - a folder that includes pre-made Nextflow run scripts for each workflow:
  * `RESUME-scautoqc-all` 
  * `RESUME-scautoqc-onlyqc`
  * `RESUME-scautoqc-afterqc`
  * `RESUME-scautoqc-untilintegrate`
  * `RESUME-scautoqc-onlyintegrate`
  * `RESUME-scautoqc-subset` 
* `bin/` - a folder that includes Python scripts used in the pipeline:
  * `gather_matrices.py` - gathers matrices from STARsolo, Velocyto and CellBender outputs (used in step 1).
  * `qc.py` - runs automatic QC workflow (used in step 2).
  * `subset.py` - subsets the input object (used in step 2a).
  * `flag_doublet.py` - runs scrublet to find doublets (used in step 3).
  * `pool_all.py` - combines all of the output objects after QC step (used in step 4).
  * `finalize_qc.py` - calculates overall QC scores per sample, adds scrublet scores (and metadata if available), and filters out bad-QC cells (used in step 5).
  * `finalize_qc_basic.py` - adds scrublet scores but doesn't remove any cells or samples (used in 5a).
  * `integration.py` - runs scVI integration (used in step 6).
* `genes_list/` - a folder that includes cell cycle, immunoglobulin and T cell receptor genes.
* `Dockerfile` - a dockerfile to reproduce the environment used to run the pipeline.

## Workflow

![](images/scautoqc-diagram-light.png#gh-light-mode-only)
![](images/scautoqc-diagram-dark.png#gh-dark-mode-only)  

This pipeline implements the scAutoQC workflow, performing automated quality control, doublet detection, and optional integration for single-cell RNA-seq data. It processes samples individually through QC steps, pools them, adds metadata and doublet flags, and finally integrates the data using scVI. Different run modes allow flexibility in executing specific parts of the workflow.

### Inputs

The pipeline requires the following primary inputs, typically configured via command-line parameters or within a RESUME script:

* INPUT 1: Legacy format (sample IDs + prefixes)
  *   A file listing the samples to be processed (`--SAMPLEFILE`).
  *   Paths to primary gene expression data:
      *   STARsolo output directories (`--ss_prefix`), with output type selected by `--ss_out` (`Gene` or `GeneFull`; default: `GeneFull`).
      *   Or Cell Ranger output directories (`--cr_prefix`) — primarily for `subset` mode or when STARsolo data is unavailable.
  *   Path to CellBender HDF5 outputs (`--cb_prefix`). If provided, used for the main expression matrix and cell/gene metadata; otherwise STARsolo/Cell Ranger outputs are used.
* INPUT 2: New CSV format (sample IDs + absolute paths)
  * Provide a CSV via `--SAMPLEFILE` with these columns:
    - required: `sampleID`
    - one of: `path_to_starsolo` or `path_to_cellranger`
    - optional: `path_to_cellbender`, `cell_or_nuclei`
  * `cell_or_nuclei` is a user‑provided annotation (values: `cell` or `nuclei`) used to select default QC thresholds; it is not inferred. If omitted, all samples are treated as `cell`.
  * STARsolo output type to be used: `Gene` or `GeneFull` (`--ss_out`, default: `GeneFull`).

*   Gather mode (`--gather_mode`): `starsolo` or `cellbender` (default: `cellbender`).
*   QC mode (`--qc_mode`): `original`, `multires`, or `combined` (default: `original`).
*   CellTypist model (`--celltypist_model`): `gut` preset or a path to a custom model.
*   (optional) Cell‑level metadata CSV to add after pooling (`--metadata`).
*   (optional) CSV with custom QC thresholds for the GMM (`--metrics_csv`).
*   CSV for manual QC cutoffs when running in `subset` mode (`--limits_csv`) (see `example_cutoffs.csv`).

### Outputs

The pipeline generates several outputs, organized into a results directory (e.g., `scautoqc-results-<project_tag>`):

*   ***[output 1]:*** Individual H5AD objects per sample with raw, spliced, unspliced, and ambiguous layers (`1_gathered_objects/`).
*   ***[output 2]:*** Individual H5AD objects per sample after QC, containing QC metrics and filtering flags (`2_qc_objects/`). In `subset` mode, this contains QC metrics but no automatic filtering flags.
*   ***[output 3]:*** QC plots generated for each individual sample during the `run_qc` step (`2_qc_plots_individual/`). Not generated in `subset` mode.
*   ***[output 4]:*** CSV files containing predicted doublet scores from Scrublet for each sample (`3_doublet_scores/`).
*   ***[output 5]:*** QC thresholds calculated in QC step for all samples (`qc_thresholds.csv`). 
*   ***[output 6]:*** The pooled H5AD object with the QC scores, the metadata (if provided) and doublet scores/flags (`scautoqc_pooled.h5ad`).
*   ***[output 7]:*** Filtered version of the previous object (`scautoqc_pooled_filtered.h5ad`).
*   ***[output 8]:*** Overall QC plots generated after pooling and metadata addition (`5_qc_plots_overall/`). Structure differs slightly in `subset` mode.
 *   ***[output 9]:*** A CSV file summarizing sample-level QC pass rates based on different mitochondrial thresholds (`sample_passqc_df.csv`). Not generated in `subset` mode.
*   ***[output 10]:*** The final, integrated H5AD object after running scVI (`scautoqc_integrated.h5ad`). Generated only if integration steps are run.
*   ***[output 11]:*** ELBO plot from the scVI model training (`5_qc_plots_overall/`). Generated only if integration steps are run.

### Run Modes

![](images/workflow_modes.png)  
The default version of the pipeline runs all the steps shown the diagram above. This pipeline has six run modes as shown above:
* `all`: runs all steps (1-2-3-4-5-6)
* `only_qc`: runs the steps until pooling including doublet finding (1-2-3)
* `after_qc`: runs the steps starting from pooling (4-5-6)
* `until_integrate`: runs the steps until integration (1-2-3-4-5)
* `only_integrate`: runs the integration step only (6)
* `subset`: runs the pipeline differently than the other modes (2a-3-4-5a)  
  * `subset_object` step replaces the `run_qc` step. This mode doesn't run automatic QC; instead, it calculates QC metrics and subsets according to the cutoffs you provide in the `--limits_csv` parameter in the 2nd step (2a).
  * `finalize_qc_basic` step replaces the `finalize_qc` step. The only difference in this step is that it doesn't remove any cells or samples (5a).

The parameters needed for all run modes are already specified in different RESUME scripts in `examples` folder, and also can be found below:

<details>
<summary><strong>INPUT MODE: Legacy format</strong></summary>

<details>
<summary>Workflow: all</summary>

```
# to run all the steps
nextflow run cellgeni/nf-scautoqc -r main \
  -entry all \            # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --metadata /path/to/metadata/file \
  --ss_prefix /path/to/starsolo-results \
  --cb_prefix /path/to/cellbender-results \
  --ss_out GeneFull \         # to specify which STARsolo output folder to use (Gene or GeneFull)
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --batch_key sampleID \  # batch key to use in scVI integration
  --ansi-log false \
  -resume
```
</details>


<details>
<summary>Workflow: only_qc</summary>

```
# to run all the steps before pooling
nextflow run cellgeni/nf-scautoqc -r main \
  -entry only_qc \            # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --metadata /path/to/metadata/file \
  --ss_prefix /path/to/starsolo-results \
  --cb_prefix /path/to/cellbender-results \
  --ss_out Gene \         # to specify which STARsolo output folder to use (Gene or GeneFull)
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --batch_key sampleID \  # batch key to use in scVI integration
  --ansi-log false \
  -resume
```

</details>


<details>
<summary>Workflow: after_qc</summary>

```
# to run after qc steps 
nextflow run cellgeni/nf-scautoqc -r main \
  -entry after_qc \            # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --postqc_path /path/to/postqc/objects \
  --scrublet_path /path/to/scrublet/csvs \
  --metadata /path/to/metadata/file \
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --batch_key sampleID \  # batch key to use in scVI integration
  --ansi-log false \
  -resume
```

</details>


<details>
<summary>Workflow: until_integrate</summary>

```
# to run steps until integration
nextflow run cellgeni/nf-scautoqc -r main \
  -entry until_integrate \            # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --metadata /path/to/metadata/file \
  --ss_prefix /path/to/starsolo-results \
  --cb_prefix /path/to/cellbender-results \
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --ansi-log false \
  -resume
```

</details>


<details>
<summary>Workflow: only_integrate</summary>

```
# to run the integration step only
nextflow run cellgeni/nf-scautoqc -r main \
  -entry only_integrate \            # to choose run mode
  --path_for_scvi /path/to/object/to/integrate \
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --batch_key sampleID \  # batch key to use in scVI integration
  --ansi-log false \
  -resume
```

</details>


<details>
<summary>Workflow: subset</summary>

```
# to run all the steps without automatic qc
nextflow run cellgeni/nf-scautoqc -r main \
  -entry subset \                             # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --metadata /path/to/metadata/file \
  --cr_prefix /path/to/cellranger/folder \    # to specify the cellranger output path
  --limits_csv /path/to/limits/file \         # to specify the cutoffs used for subsetting
  --project_tag project1 \         # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --ansi-log false \
  -resume
```

</details>

</details>


<details>
<summary><strong>INPUT MODE: New format</strong></summary>

<details>
<summary>Workflow: all</summary>

```
# to run all the steps
nextflow run cellgeni/nf-scautoqc -r main \
  -entry all \            # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --metadata /path/to/metadata/file \
  --celltypist_model gut \
  --gather_mode cellbender \
  --qc_mode original \
  --ss_out GeneFull \     # to specify which STARsolo output folder to use (Gene or GeneFull)
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --batch_key sampleID \  # batch key to use in scVI integration
  --ansi-log false \
  -resume
```
</details>


<details>
<summary>Workflow: only_qc</summary>

```
# to run all the steps before pooling
nextflow run cellgeni/nf-scautoqc -r main \
  -entry only_qc \            # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --metadata /path/to/metadata/file \
  --ss_prefix /path/to/starsolo-results \
  --cb_prefix /path/to/cellbender-results \
  --ss_out GeneFull \         # to specify which STARsolo output folder to use (Gene or GeneFull)
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --batch_key sampleID \  # batch key to use in scVI integration
  --ansi-log false \
  -resume
```

</details>


<details>
<summary>Workflow: after_qc</summary>

```
# to run after qc steps 
nextflow run cellgeni/nf-scautoqc -r main \
  -entry after_qc \            # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --postqc_path /path/to/postqc/objects \
  --scrublet_path /path/to/scrublet/csvs \
  --metadata /path/to/metadata/file \
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --batch_key sampleID \  # batch key to use in scVI integration
  --ansi-log false \
  -resume
```

</details>


<details>
<summary>Workflow: until_integrate</summary>

```
# to run steps until integration
nextflow run cellgeni/nf-scautoqc -r main \
  -entry until_integrate \            # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --metadata /path/to/metadata/file \
  --ss_prefix /path/to/starsolo-results \
  --cb_prefix /path/to/cellbender-results \
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --ansi-log false \
  -resume
```

</details>


<details>
<summary>Workflow: only_integrate</summary>

```
# to run the integration step only
nextflow run cellgeni/nf-scautoqc -r main \
  -entry only_integrate \            # to choose run mode
  --path_for_scvi /path/to/object/to/integrate \
  --project_tag test1 \   # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --batch_key sampleID \  # batch key to use in scVI integration
  --ansi-log false \
  -resume
```

</details>


<details>
<summary>Workflow: subset</summary>

```
# to run all the steps without automatic qc
nextflow run cellgeni/nf-scautoqc -r main \
  -entry subset \                             # to choose run mode
  --SAMPLEFILE /path/to/sample/file \
  --metadata /path/to/metadata/file \
  --cr_prefix /path/to/cellranger/folder \    # to specify the cellranger output path
  --limits_csv /path/to/limits/file \         # to specify the cutoffs used for subsetting
  --project_tag project1 \         # to specify the run to add to the end of output folder (e.g. scautoqc-results-test1)
  --ansi-log false \
  -resume
```

</details>

</details>

### 1. `gather_matrices`  

The inputs for the first step are determined according to STARsolo output which is specified via `--ss_out`. By default, the "GeneFull" folder is used; "Gene" can be selected to focus on reads not mapped to introns.

This step can use up to three inputs:
  * STARsolo output folder selected by `--ss_out` (`Gene` or `GeneFull`)
  * STARsolo “Velocyto” folder (to add spliced/unspliced/ambiguous layers when available)
  * CellBender output in H5 format (if provided, used for main expression matrix and metadata)

`gather_matrices` combines these into one h5ad per sample with multiple layers: raw, spliced, unspliced, ambiguous. Main expression matrix and cell/gene metadata are taken from CellBender when present; otherwise STARsolo/Cell Ranger are used. The raw layer is taken from the STARsolo folder selected by `--ss_out`. Spliced, unspliced and ambiguous layers are added from the “Velocyto” folder when available, irrespective of using Gene or GeneFull. The cell barcode set can come from STARsolo+CellBender (`gather_mode=starsolo`) or only from CellBender (`gather_mode=cellbender`, default).

This step can also use Cell Ranger inputs if `--cr_prefix` is provided or the CSV header includes `path_to_cellranger` instead of `path_to_starsolo`; Cell Ranger runs will not include Velocyto layers.

This step produces:  
* ***[output 1]:*** h5ad object with different layers

### 2. `run_qc`

This step requires the output of `gather_matrices` step which is the h5ad object with one or more layers.  

`run_qc` step optionally runs CellTypist when the `celltypist_model` parameter is provided. If `celltypist_model` is set to "gut", the pipeline loads a curated set of gut-focused CellTypist models. Alternatively, `celltypist_model` can be a path to a custom CellTypist model. The "gut" preset includes the following models:
*  **cecilia22_predH:** Immune populations across 20 tissues (32 cell types) (ref: [Domínguez-Conde et al, 2022](https://doi.org/10.1126/science.abl5197))
*  **cecilia22_predL:** Immune subpopulations across 20 tissues (98 cell types) (ref: [Domínguez-Conde et al, 2022](https://doi.org/10.1126/science.abl5197))
*  **elmentaite21_pred:** Intestinal cells from fetal, pediatric and adult gut studies (134 cell types) (ref: [Elmentaite et al, 2021](https://doi.org/10.1038/s41586-021-03852-1))
*  **suo22_pred:** Fetal stromal and immune populations (138 cell types) (ref: [Suo et al, 2022](https://doi.org/10.1126/science.abo0510))
*  **megagut_pred:** Pan‑GI study covering all cell types (89 cell types) (ref: [Oliver et al, 2024](https://doi.org/10.1038/s41586-024-07571-1)).

After optional CellTypist annotation, `run_qc` applies the automatic QC to each sample. Three QC modes are supported and selectable via `qc_mode`:
* `original` mode uses the main automatic QC logic (summarised [here](https://teichlab.github.io/sctk/notebooks/automatic_qc.html)). It (i) computes the eight QC metrics, (ii) builds a QC embedding (PCA → neighbors → UMAP → clustering), then (iii) loops over mitochondrial upper bounds (20, 50, 80). For each bound it refits Gaussian Mixture Models (GMMs) for four core metrics (n_counts, n_genes, percent_mito, percent_spliced) (default thresholds for these four core metrics differ if the sample is single-nuc `nuclei` or single-cell `cell`) to derive sample‑specific pass ranges, assigns per‑cell pass/fail across all required metrics, and calls clusters good if ≥50% of their cells pass (via `min_frac`). The loop produces cellwise and clusterwise pass flags per mito threshold; final summary columns (`pass_auto_filter`, `good_qc_cluster`) encode the most stringent (lowest mito) threshold each cell/cluster still satisfies (smaller value = stricter pass).
* `multires` mode keeps the same metric computation and initial QC embedding but fixes the mitochondrial bound at 20 (no mito loop). It runs the GMM cellwise QC once to obtain `cell_passed_qc`, then performs multi‑resolution clustering (grid 0.1–1.0) applying clusterwise QC at each resolution. From these runs it records, for every cell, a `consensus_fraction` (fraction of resolutions where its cluster passes) and selects a consensus cut‑off that maximizes the Jaccard overlap of failing cells versus the cellwise GMM calls (explicit search over unique consensus_fraction values). Cells with consensus_fraction ≥ chosen threshold become `consensus_passed_qc`; a degeneracy guard forces all False if the fractions carry no signal, and an additional `keep_multires` flag rescues high‑agreement cells (≥0.90) or cluster passes. This yields per‑cell QC flags less sensitive to one arbitrary clustering resolution.
* `combined` mode nests the multi‑resolution procedure inside the mitochondrial loop (20, 50, 80). For each mito bound it refits the GMMs (adjusting percent_mito upper limit), performs cellwise QC (`good_qc_cell_mitoX`), runs multi‑resolution cluster QC (producing `good_qc_cluster_mitoX`, `consensus_fraction_mitoX`), searches a consensus threshold (same Jaccard strategy) to define `consensus_passed_qc_mitoX`, and applies a consensus fraction floor (0.90) plus degeneracy checks. Pass propagation merges stricter passes upward (20 → 50 → 80) and collapses to numeric summaries (`pass_auto_filter`, `consensus_pass_auto_filter`) representing the most stringent mitochondrial window retained. This delivers robust QC calls across both clustering resolution and mitochondrial burden.

A few more details on `run_qc` step:
- Thresholds for the four core metrics (n_counts, n_genes, percent_mito, percent_spliced) to use in GMM are predefined in this step. Different thresholds are used if a sample is single-nuc `nuclei `or single-cell `cell`.
- Behaviour of GMM can be modified using `gmm_cutoff` option: `inner` is the default one and uses closest points, and `outer` uses furthest points.


This step produces [output 2] and [output 3]: 
* ***[output 2]:*** an h5ad object with QC-related columns,  
* ***[output 3]:*** QC plots for each sample:  
  * UMAP plots:
    * coloured by each metric and good_qc_cluster and qc_cluster
    * coloured by probability scores, uncertain and predicted labels using the first model (cecilia22_predH)
    * coloured by good_qc_cluster, pass_auto_filter, pass_default_filter and qc_cluster
  * violin plot for each metric grouped by QC clusters  
  * scatter plots for n_counts vs each metric:
    * coloured by good_qc_cluster
    * coloured by pass_default_filter
    * coloured by pass_auto_filter
    * coloured by the probability of a CellTypist prediction using the first model (cecilia22_predH)

### 3. `find_doublets`  

This step requires the output of `run_qc` step which is the h5ad object with postqc columns.  

`find_doublets` step runs [scrublet](https://github.com/swolock/scrublet) on the h5ad object and annotates the doublet scores to cells. This step runs in parallel with step 4 for efficiency.

This step produces:  
* ***[output 4]:*** CSV file with scrublet scores for each cell barcode for each sample.

### 4. `pool_all`  

This step requires the outputs of `run_qc` step from all the samples.  

`pool_all` step combines all of the objects produced in `run_qc` step in a single h5ad object.

This step produces:
* ***[output 5]:*** CSV file with QC cutoffs from all samples as CSV

### 5. `finalize_qc`  

This step requires the h5ad output from `pool_all` and the scrublet CSV outputs from `find_doublets`.

`finalize_qc` computes sample‑level QC summaries, adds scrublet scores (and optional cell‑level metadata), and removes bad‑QC cells/samples according to the selected QC mode.

This step produces [output 6], [output 7], [output 8]:
* ***[output 6]:*** pooled h5ad object before filtering
* ***[output 7]:*** pooled h5ad object after filtering
* ***[output 8]:*** QC plots: 
  * scatter plot of total_cell_count vs passQC_count80 coloured by each sample
  * histogram of log of samples that have passQC_count20, passQC_count50 and passQC_count80 columns True
* ***[output 9]:*** CSV with percentages of the cells passed QC for each test

### 6. `integrate`  

This step requires the h5ad object from `finalize_qc` step.
`integrate` removes stringent doublets and applies scVI integration to all samples using "sampleID" as the batch key (by default), "log1p_n_counts" and "percent_mito" columns as continuous covariates. The following preprocessing is applied:  
* Stringent doublets (cells with doublet score > 0.3 and BH‑adjusted p < 0.05) are removed.  
* 7500 highly variable genes are selected (configurable via `--n_top_genes`).  
* Cell‑cycle genes are removed.  
* Latent space dimensionality is 20; batch size is 256.
Covariates used by scVI include log1p_n_counts and percent_mito.

This step produces:
* ***[output 10]:*** final integrated h5ad object
* ***[output 11]:*** ELBO plot from scVI training

## Different steps for `subset` mode

### 2a. `subset_object` (`subset` mode only)

This step requires the Cell Ranger outputs specified by the `--cr_prefix` parameter and a CSV file containing the subsetting cutoffs provided via the `--limits_csv` parameter.

The `subset_object` step is exclusive to the `subset` mode and replaces the `run_qc` step from the main pipeline. Unlike the main pipeline, it does not execute the automatic QC algorithm, which is a core component of the scAutoQC pipeline. Instead, it calculates QC metrics and subsets the data based on the cutoffs defined in the `--limits_csv` parameter.

This step produces:
* ***[output 2]:*** an h5ad object with QC metrics.

### 5a. `finalize_qc_basic` (`subset` mode only)

This step requires the h5ad output from the `pool_all` step and the scrublet CSV outputs from the `find_doublets` step.

The `add_metadata_basic` step is also exclusive to the `subset` mode and replaces the `add_metadata` step from the main pipeline. The key differences are that it does not perform QC scoring per sample and does not remove any cells or samples.

## Original workflow scheme

![](images/scautoqc-original-diagram.png)

# Changelog

For a version history/changelog, please see the [CHANGELOG file](CHANGELOG.md).