# Changelog

### 25-314
This release introduces significant changes to the pipeline.  
⚠️ Note: h5ad outputs from this version might **not be backward‑compatible** with previous scAutoQC releases.

#### ➤ New Input Format
Two input methods are now supported:
1. Legacy method  
   - `SAMPLEFILE`: list of sample IDs  
   - `ss_prefix`: STARsolo output directory  
   - `cb_prefix`: Cellbender output directory  
   - `cr_prefix`: Cell Ranger output directory (alternative to STARsolo)
2. New method  
   - CSV file with **sample IDs and absolute paths** (see `example_input.csv`)  
   - The CSV may include a `cell_or_nuclei` column; its value is propagated to the AnnData (`ad.uns['cell_or_nuclei']`) and is used to select default QC thresholds. This is an explicit annotation from the metadata, not automatic detection.

This new format increases flexibility and simplifies setup.

#### ➤ Major Workflow & Script Changes
**gather_matrices**
- Introduced gather modes:  
  - `starsolo`: uses barcodes from STARsolo + Cellbender (original behavior).  
  - `cellbender`: uses Cellbender barcodes only (new default).
- STARsolo default folder is now `"GeneFull"` (was `"Gene"`). 
- Velocyto layers are now included whenever available, regardless of using Gene vs GeneFull (previously GeneFull runs could miss these layers).
- Reads the `cell_or_nuclei` annotation from the new CSV input and stores it in the output object for downstream QC thresholding, otherwise it assumes all samples as `cell`.
- Removed the local `read_cellbender` shim; the pipeline now relies on the updated `sctk` implementation which supports current CellBender outputs.
- Increased default memory allocation (8 GB → 16 GB) and removed unused flags.

**run_qc**
- Added QC modes (select via `--qc_mode`):  
  - `original`: auto‑QC with mitochondrial thresholds (used in Pan‑GI Atlas).  
  - `multires`: auto‑QC with multi‑resolution clustering and consensus selection.  
  - `combined`: runs the multi‑resolution consensus inside the mitochondrial‑threshold loop.
- Separate default thresholds for single‑cell vs single‑nuc based on the `cell_or_nuclei` annotation.
- New option to choose GMM cutoff strategy: `gmm_cutoff = inner|outer`.
- Support for custom thresholds via `--metrics_csv`.
- Integrated CellTypist models (use `--celltypist_model <tissue>` or `--celltypist_model /path/to/model.pkl`).
- Refactoring with clearer logging and explicit keys stored in `AnnData.obs/uns`.

**subset_object**
- More robust input discovery
- The output h5ad name now includes the sample ID.

**pool_all**
- Now outputs a consolidated `qc_thresholds.csv`.  
- The pooled h5ad object is kept internal and passed directly to the next step.

**add_metadata → finalize_qc**
- Step renamed for clarity.  
- QC filtering applied according to the selected QC mode.  
- Output changes:  
  - `scautoqc_pooled.h5ad`: unfiltered object with all metadata (moved from `pool_all`).  
  - `scautoqc_pooled_filtered.h5ad`: renamed from `scautoqc_pooled_doubletflagged_metaadded.h5ad`; contains the filtered dataset (same criteria as before).

**add_metadata_basic → finalize_qc_basic**
- Output changes:  
  - `scautoqc_pooled_basic.h5ad`: unfiltered object with all metadata and doublet scores.

**integrate**
- New parameter for the number of top genes used in integration (`n_top_genes`).  
- New `from_scautoqc` parameter to allow integrating arbitrary input objects (when false, it does not remove cell‑cycle genes or stringent doublets).  

#### ➤ System Updates
- Updated Dockerfile and Singularity image for improved stability and compatibility.
- Smarter memory requests.
- Added tag support so Nextflow now shows which sample each process runs.
- Improved output file/folder names.
- Added timestamps to timeline, report and trace files.

---

### 25-101

#### ➤ Bug Fixes
- `qc.py`: correctly applies QC metrics based on the presence of spliced/unspliced layers.

#### ➤ Stability & Resilience
- `add_metadata`: no longer fails when metadata input is absent.
- `add_metadata_basic`: fully supports RESUME for interrupted runs.

---

### 25-091
This release introduces a new workflow and multiple pipeline enhancements.  

#### ➤ New Workflow
- **`subset` workflow**: Allows running the pipeline by subsetting objects using predefined cutoffs instead of the automatic QC workflow (steps 2a-3-4-5a).  
  - More details in [README](README.md#different-steps-for-subset-mode).

#### ➤ Pipeline Improvements
- Added support for **Cell Ranger outputs** (in addition to STARsolo).  
- Updated **Singularity image** for better compatibility.  
- Renamed certain **output files** for clarity.  
- Improved **RESUME functionality** for reliability.  
- Introduced **smart memory allocation** for `pool_all` and `add_metadata` based on input size.  
- Optimised **resource allocation** across processes.  
- Full support for **symbolic links** as input.

#### ➤ Script Optimizations
- Removed unused code, characters, and packages.  
- Fixed hardcoded paths for better flexibility.  
- Reduced memory usage in `pool_all`.

---

### 25-064
#### ➤ New Workflows
- **`until_integrate`**: Run steps 1–5 (until integration).  
- **`only_integrate`**: Run only step 6 (integration).  

#### ➤ Script Updates
- Output folder names renamed for consistency.

---

### 24-143
#### ➤ New Workflow
- **`only_qc`**: Run the pipeline until pooling.  
  - Useful for processing sample sets separately, then combining them later with `after_qc`.

#### ➤ Improvements
- Pipeline can run without **Cellbender outputs** (uses STARsolo only if available).  
- Reverted removal of metadata columns in `integration.py`.  
- Added new **RESUME scripts**.  
- `after_qc` workflow can now use only the **samples listed** (instead of all in the folder).  
- **Reports folder** now named consistently with results folder.  
- Fixed **plot outputs** in `run_qc` step.

---

### 24-088
#### ➤ Features
- Added support for **single-nuclei samples**.  

#### ➤ Script Updates
- `integration.py`: now removes columns created in earlier steps.  
- RESUME scripts reorganized.

---

### 24-064
#### ➤ New Workflow
- **`after_qc`**: Enables re-using previously processed samples with scAutoQC.  

#### ➤ Improvements
- Added new **RESUME script** for `after_qc`.  
- `integration.py` optimized for efficiency.  
- Removed parameters: `ss_matrix`, `covar_keys` (may be reintroduced later).

---

### 24-060
- **Bug fixes.**

---

### 24-059
- **Typo fixes.**

---

### 24-058
#### ➤ Script Improvements
- Unified output filenames across scripts.  

**gather_matrices.py**
- Added support for **Cellbender v3 outputs**.  

**qc.py**
- Moved `sampleID` annotation step from `pool_all.py`.  

**pool_all.py**
- Added **RESUME functionality** (Nextflow).  
- Removed `sampleID` annotation (moved to `qc.py`).  

**add_scrublet_meta.py**
- Added cleanup step to remove empty metadata columns.  

**integration.py**
- Added support for **batch key** and **categorical covariates**.  
  - Default batch key: `sampleID`.  
  - Covariates: none by default.  

**main.py**
- Fixed output generation issue.  
- Modified `find_doublets` to run only on good QC samples.  
- Added **`batch_key`** and **`covar_keys`** parameters for scVI integration.  

---

### 24-050
- **First release.**