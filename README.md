# nf-scautoqc

nf-scautoqc is the Nextflow implementation of scAutoQC pipeline used in Oliver et al, 2023. 

## Contents of Repo:
* `main.nf` - the Nextflow pipeline that executes STARsolo.
* `nextflow.config` - the configuration script that allows the processes to be submitted to IBM LSF on Sanger's HPC and ensures correct environment is set via singularity container (this is an absolute path). Global default parameters are also set in this file and some contain absolute paths.
* `RESUME-starsolo` - an example run script that executes the pipeline it has four hardcoded argument: `/path/to/sample/file`, `/path/to/starsolo-results`, `/path/to/cellbender-results`, (optional: `/path/to/metadata/file`) that need to be changed based on your local set up.
* `bin/gather_matrices.py` - a Python script that gathers matrices from STARsolo, Velocyto and Cellbender outputs.
* `bin/qc.py` - a Python script that runs automatic QC workflow.
* `bin/flag_doublet.py` - a Python script that runs scrublet to find doublets.
* `bin/pool_all.py` - a Python script that combines all of the output objects after QC step.
* `bin/add_scrublet_meta.py` - a Python script that adds scrublet scores (and metadata if available).
* `bin/integration.py` - a Python script that runs scVI integration.
* `genes_list/` - a folder that includes cell cycle, immunoglobulin and T cell receptor genes.
* `Dockerfile` - a dockerfile to reproduce the environment used to run the pipeline.

## Workflow

![](scautoqc-diagram.png)

### 1. `gather_matrices`  

This step requires three inputs: 
* STARsolo output folder named "Gene"
* STARsolo output folder named "Velocyto"
* Cellbender output in h5 format (if Cellbender output doesn't exist, change mode to cb+normal) (option to change mode will be added in the future)  

`gather_matrices` step combines the matrices from three inputs into one h5ad object with four layers (raw, spliced, unspliced, ambiguous). Main expression matrix, cell and gene metadata are retrieved from Cellbender output. Raw matrix is retrieved from the expression matrix of STARsolo output folder named Gene. Spliced, unspliced and ambiguous matrices are all retrieved from the expression matrices of STARsolo output folder named Velocyto.

### 2. `run_qc`

This step requires the output of `gather_matrices` step which is the h5ad object with four layers.  

`run_qc` step uses main automatic QC workflow which is summarised [here](https://teichlab.github.io/sctk/notebooks/automatic_qc.html).

### 3. `find_doublets`  

This step requires the output of `run_qc` step which is the h5ad object with postqc columns.  

`find_doublets` step runs [scrublet](https://github.com/swolock/scrublet) on the h5ad object and annotates the doublet scores to cells. This step runs simultaneously with step 4 for efficiency.

### 4. `pool_all`  

This step requires the outputs of `run_qc` step from all the samples.  

`pool_all` step combines all of the objects produced in `run_qc` step in a single h5ad object.

### 5. `add_metadata`  

This step requires the h5ad output from `pool_all` and the scrublet csv outputs from `find_doublets` steps.  

`add_metadata` step gives scores of different QC metrics of each sample, and adds scrublet scores (and the cell metadata according to samples if provided) to the h5ad object. (QC-related plots are generated but are not copied to output directory, this will be fixed in the future.)  

### 6. `integrate`  

This step requires the h5ad object from `add_metadata` step.
`integrate` step removes stringent doublets (doublet score higher than 0.3, and bh score lower than 0.05) applies scVI integration to all samples by using "sampleID" as a batch key, and "log1p_n_counts" and "percent_mito" columns as categorical covariates. The final integrated object is given as the output of all of this pipeline. The steps below are applied before running integration:  
* Stringent doublets are removed.  
* 7500 higly variable genes are chosen.  
* All cell cycle genes are removed.  
* The dimensionality of the latent space is chosen as 20.
* The batch size is chosen as 256.

## Original workflow scheme

![](scautoqc-original-diagram.png)
