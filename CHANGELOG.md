## Changelog

### v0.5.0
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
* Updates in README
  * The workflow diagram has been recreated.
  * The workflow modes have been described with a new figure.
  * Outputs from each step have been described in detail.
  * Text in some steps were revised.

### v0.4.0
* Added support for single-nuc samples
* Improvements in scripts
  * integration.py now removes the columns which were created in previous steps
  * RESUME scripts have been reorganised

### v0.3.0
* <ins>**New workflow:**</ins> `after_qc`
  * It is now easier to work with the samples which has been processed with scAutoQC pipeline before. 
* Improvements in nextflow pipeline and python scripts
  * Created new RESUME script for afterqc workflow.
  * integration.py should now work more efficiently.
  * ss_matrix and covar_keys parameters have been removed (they will be considered for the future releases). 

### v0.2.2
* Bug fixes

### v0.2.1
* Typo fix

### v0.2.0
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

### v0.1.0
* First release