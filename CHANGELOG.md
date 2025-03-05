## Changelog

### 25-064
* Added two new workflows:
  * `until_integrate` workflow makes it easier to run the steps until integration (1-2-3-4-5)
  * `only_integrate` workflow makes it easier to run the integration step only (6)
* Improvements and changes in scripts:
  * Folder names in the outputs were renamed.


### 24-143
* <ins>**New workflow:**</ins> `only_qc`
  * It is now easier to run the pipeline until the pooling step. 
  * This can be used to process different sets of samples in different times, then all the outputs from this workflow  can be used together with `after_qc` mode.
* Improvements and changes in scripts:
  * It is now possible to run the pipeline without Cellbender output, STARsolo output is used only in this case.
  * Removing metadata columns in integration.py was reverted back
  * New RESUME scripts have been added
  * `after_qc` workflow can now use the samples in the sample list only rather than all the objects in the folder.
  * Reports folder is now named similar to the results folder.
  * Outputting plots in `run_qc` step now works as expected.

### 24-088
* Added support for single-nuc samples
* Improvements in scripts
  * integration.py now removes the columns which were created in previous steps
  * RESUME scripts have been reorganised

### 24-064
* <ins>**New workflow:**</ins> `after_qc`
  * It is now easier to work with the samples which has been processed with scAutoQC pipeline before. 
* Improvements in nextflow pipeline and python scripts
  * Created new RESUME script for afterqc workflow.
  * integration.py should now work more efficiently.
  * ss_matrix and covar_keys parameters have been removed (they will be considered for the future releases). 

### 24-060
* Bug fixes

### 24-059
* Typo fix

### 24-058
* Improvements in python scripts
  * Changed output filenames in all scripts
  * gather_matrices.py
    - added support for Cellbender v3 outputs
  * qc.py
    - moved sampleID annotation step from pool_all.py
  * pool_all.py
    - added support for resume functionality of nextflow
    - removed sampleID annotation step (moved to qc.py)
  * add_scrublet_meta.py
    - added an extra step to remove empty columns in input metadata
  * integration.py
    - added support to specify batch key and categorical covariate keys (default: sampleID for batch key, empty for covariates)
* Improvements in main.py
  * fixed output generation problem
  * modified find_doublets process so it only runs for good QC samples
  * added batch_key parameter to specify the column from cell metadata for scVI integration (default: sampleID)
  * added covar_keys parameter to specify the columns from cell metadata for scVI integration (no default)
* Updated workflow figure in README.md

### 24-050
* First release