#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process gather_matrices {

  tag "${samp}"

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/1_gathered_objects", mode: 'copy', saveAs: {filename -> "${samp}.${filename}"}

  input:
  val(samp)
  path(cr_gene), stageAs: "filteredSS"
  path(cr_velo), stageAs: "filteredVelo"
  path(cb_h5)
  val(cell_or_nuclei)

  output:
  tuple val(samp), path("*gene_velo_cellbender.filtered.h5ad"), emit: obj

  script:
  """
  python ${projectDir}/bin/gather_matrices.py --cr_gene ${cr_gene} --cr_velo ${cr_velo} --cb_h5 ${cb_h5} --cell_or_nuclei ${cell_or_nuclei} --gather_mode ${params.gather_mode}
  """
}

process run_qc {

  tag "${samp}"

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/2_qc_objects", pattern: '*.h5ad', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/2_qc_plots_individual", pattern: '*.png', mode: 'copy'

  input:
  tuple val(samp), path(gath_out)

  output:
  tuple val(samp), path("*_postqc.h5ad"), path("*.csv"), emit: samp_obj
  path("*.png")

  script:
  """
  python ${projectDir}/bin/qc.py --sample_id ${samp} --metrics_csv ${params.metrics_csv} --celltypist ${params.celltypist_model} --qc_mode ${params.qc_mode} --gath_obj ${gath_out} --gmm_cutoff ${params.gmm_cutoff}
  """
}

process subset_object {

  tag "${samp}"

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/2_qc_objects", pattern: '*.h5ad', mode: 'copy', saveAs: {filename -> "${samp}_${filename}"}

  input:
  val(samp)

  output:
  tuple val(samp), path("*_subsetted.h5ad"), emit: samp_obj

  script:
  """
  python ${projectDir}/bin/subset.py --sample_id ${samp} --cr_prefix ${params.cr_prefix} --limits_csv ${params.limits_csv}
  """
}

process find_doublets {

  tag "${samp}"

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/3_doublet_scores", pattern: '*.csv', mode: 'copy'

  input:
  tuple val(samp), path(qc_out), path(qc_thres)

  output:
  tuple val(samp), path("*_scrublet.csv")

  script:
  def filter_column = ""
  if (params.qc_mode == "multires") {
    filter_column = "cluster_passed_qc"
  } else {
    filter_column = "good_qc_cluster_mito80"  
  }
  
  """
  python ${projectDir}/bin/flag_doublet.py --filter ${filter_column} --samp ${samp} --input ${qc_out}
  """

}


process pool_all {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.csv', mode: 'copy'

  memory {
      def numSamples = samp.size()
      def baseMemory = 10.GB
      def perSampleMemory = 1.GB * numSamples
      def requestedMemory = [baseMemory, perSampleMemory].max()
      return (requestedMemory * 1.4 * task.attempt) as nextflow.util.MemoryUnit
  }

  input:
  val(samp)
  path(qc_out)
  path(ranges_out)
  val(numSamples)

  output:
  path("scautoqc_pooled0.h5ad"), emit: obj
  path("qc_thresholds.csv"), optional: true
  val(numSamples), emit: numSamples

  script:
  """
  python ${projectDir}/bin/pool_all.py --samples ${samp.join(",")} --objects ${qc_out.join(",")} --ranges ${ranges_out.join(",")}
  """
}

process finalize_qc {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.h5ad', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/5_6_qc_plots_overall", pattern: '*.png', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.csv', mode: 'copy'

  memory {
      def inputSizeGB = pool_out.size() / (1024 * 1024 * 1024)
      def baseMemory = 8.GB
      def memoryMultiplier = 3.5  // Conservative multiplier between observed values
      def requestedMemory = baseMemory + (inputSizeGB * memoryMultiplier).GB
      return (requestedMemory * task.attempt) as nextflow.util.MemoryUnit
  }

  input:
  path(pool_out)
  path(scr_out)
  val(numSamples)

  output:
  path("scautoqc_pooled_filtered.h5ad"), emit: obj
  path("scautoqc_pooled.h5ad")
  path("*.png")
  path("sample_passqc_df.csv")

  script:
  """
  export BASE_DIR=${projectDir}
  python ${projectDir}/bin/finalize_qc.py --obj ${pool_out} --scr ${scr_out.join(",")} --meta ${params.metadata} --qc_mode ${params.qc_mode}
  """
}

process finalize_qc_basic {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.h5ad', mode: 'copy'

  memory {
      def inputSizeGB = pool_out.size() / (1024 * 1024 * 1024)
      def baseMemory = 8.GB
      def memoryMultiplier = 3.5  // Conservative multiplier between observed values
      def requestedMemory = baseMemory + (inputSizeGB * memoryMultiplier).GB
      return (requestedMemory * task.attempt) as nextflow.util.MemoryUnit
  }
    
  input:
  path(pool_out)
  path(scr_out)
  val(numSamples)

  output:
  path("scautoqc_pooled_basic.h5ad"), emit: obj

  script:
  """
  export BASE_DIR=${projectDir}
  python ${projectDir}/bin/finalize_qc_basic.py --obj ${pool_out} --scr ${scr_out.join(",")} --meta ${params.metadata}
  """
}

process integrate {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}", pattern: '*.h5ad', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/5_6_qc_plots_overall", pattern: '*.png', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/6_models", pattern: '*.pkl', mode: 'copy'

  memory {
      def inputSizeGB = qc2_out.size() / (1024 * 1024 * 1024)
      def requestedMemory = (inputSizeGB * 15).GB
      return (requestedMemory * task.attempt) as nextflow.util.MemoryUnit
  }

  input:
  path(qc2_out)

  output:
  path("scautoqc_integrated.h5ad")
  path("*.png")
  path("*.pkl")

  script:
  """
  python ${projectDir}/bin/integration.py --obj ${qc2_out} --batch ${params.batch_key} --n_top_genes ${params.n_top_genes} --from_scautoqc ${params.from_scautoqc}
  """
}

// Function to detect CSV format and create standardized channels
def createInputChannels(samplefile) {
    println "Using ${params.ss_out} output"
    // Read first line to detect format
    def firstLine = file(samplefile).readLines()[0]
    def hasHeaders = firstLine.contains('sampleID')     
    if (hasHeaders) {
        // New format: CSV with full column headers
        println "Detected the CSV format with headers"
        return Channel.fromPath(samplefile)
            .splitCsv(header: true)
            .multiMap { row ->
                samp: row.sampleID
                cr_gene: firstLine.contains('path_to_cellranger') ? 
                    "${row.path_to_cellranger}/filtered_feature_bc_matrix.h5" : 
                    "${row.path_to_starsolo}/output/${params.ss_out}/filtered/"
                cr_velo: row.path_to_starsolo ? 
                    "${row.path_to_starsolo}/output/Velocyto/raw/" : []
                cb_h5: row.path_to_cellbender ? 
                    "${row.path_to_cellbender}/${params.cellbender_input}" : []
                ss_gene: row.path_to_starsolo ? 
                    "${row.path_to_starsolo}/output/${params.ss_out}/filtered/" : []
                cell_or_nuclei: row.cell_or_nuclei ?: 'cell'
            }
    } else {
        // Legacy format: CSV with just sampleIDs + prefixes
        println "Detected the legacy format (CSV with sample IDs + prefixes)"
        return Channel.fromPath(samplefile)
            .splitCsv(header: false)
            .flatten()
            .map { it ->
                def prefix = params.ss_prefix ? params.ss_prefix : params.cr_prefix
                def resolvedPath = "readlink -f ${prefix}/${it}".execute().text.trim()
                [it, resolvedPath]
            }
            .multiMap { it, resolvedPath ->
                samp: it
                cr_gene: "${params.ss_prefix}" == "" ? 
                    "${resolvedPath}/filtered_feature_bc_matrix.h5" : 
                    "${resolvedPath}/output/${params.ss_out}/filtered/"
                cr_velo: "${params.ss_prefix}" == "" ? [] : 
                    "${resolvedPath}/output/Velocyto/raw/"
                cb_h5: "${params.cb_prefix}" == "" ? [] : 
                    "${params.cb_prefix}/${it}/${params.cellbender_input}"
                ss_gene: "${params.ss_prefix}" == "" ? [] : 
                    "${resolvedPath}/output/${params.ss_out}/filtered/"
                cell_or_nuclei: "${params.cell_or_nuclei}" // default value for legacy format
            }
    }
}


workflow all {
  def samples = createInputChannels("${params.SAMPLEFILE}")

  gather_matrices(samples.samp, samples.cr_gene, samples.cr_velo, samples.cb_h5, samples.cell_or_nuclei)
  run_qc(gather_matrices.out.obj)
  find_doublets(run_qc.out.samp_obj)

  def samp_collected = run_qc.out.samp_obj.collect(){ it[0] }
  def numSamples = samp_collected.map { it.size() }

  pool_all(samp_collected, run_qc.out.samp_obj.collect() { it[1] }, run_qc.out.samp_obj.collect() { it[2] }, numSamples)
  finalize_qc(pool_all.out.obj, find_doublets.out.collect(){ it[1] }, pool_all.out.numSamples)
  integrate(finalize_qc.out.obj)
}

workflow only_qc {
  def samples = createInputChannels("${params.SAMPLEFILE}")

  gather_matrices(samples.samp, samples.cr_gene, samples.cr_velo, samples.cb_h5, samples.cell_or_nuclei)
  run_qc(gather_matrices.out.obj)
  find_doublets(run_qc.out.samp_obj)
}

workflow after_qc {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv (header: false) 
       .flatten()
       .set {samples}
  Channel.fromPath("${params.postqc_path}/*.h5ad")
       .flatten()
       .set {objects}
  Channel.fromPath("${params.scrublet_path}/*.csv")
       .flatten()
       .set {scrublets}
  pool_all(run_qc.out.samp_obj.collect(){ it[0] }, run_qc.out.samp_obj.collect() { it[1] }, run_qc.out.samp_obj.collect() { it[2] })
  finalize_qc(pool_all.out, find_doublets.out.collect(){ it[1] }, pool_all.out.numSamples)
  integrate(finalize_qc.out.obj)
}

workflow until_integrate {
  def samples = createInputChannels("${params.SAMPLEFILE}")

  gather_matrices(samples.samp, samples.cr_gene, samples.cr_velo, samples.cb_h5, samples.cell_or_nuclei)
  run_qc(gather_matrices.out.obj)
  find_doublets(run_qc.out.samp_obj)

  def samp_collected = run_qc.out.samp_obj.collect(){ it[0] }
  def numSamples = samp_collected.map { it.size() }

  pool_all(samp_collected, run_qc.out.samp_obj.collect() { it[1] }, run_qc.out.samp_obj.collect() { it[2] }, numSamples)
  finalize_qc(pool_all.out.obj, find_doublets.out.collect(){ it[1] }, pool_all.out.numSamples)
}

workflow only_integrate {
  if (params.path_for_scvi.endsWith('.csv')) {
    println "Detected the CSV file with multiple objects to integrate"
    obj_scvi = Channel.fromPath("${params.path_for_scvi}")
      .splitCsv(header: false)
      .flatten()
      .map { row -> file(row) }
  } else {
    obj_scvi = Channel.fromPath("${params.path_for_scvi}")
  }
  integrate(obj_scvi)
}

workflow subset {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv (header: false) 
       .flatten()
       .set {samples}
  subset_object(samples)

  find_doublets(subset_object.out.samp_obj)

  def samp_collected = subset_object.out.samp_obj.collect(){ it[0] }
  def numSamples = samp_collected.map { it.size() }

  pool_all(samp_collected, subset_object.out.samp_obj.collect() { it[1] }, subset_object.out.samp_obj.collect() { it[2] }, numSamples)
  finalize_qc_basic(pool_all.out, find_doublets.out.collect(){ it[1] }, pool_all.out.numSamples)
}

workflow subset_new {
  params.rng = "$projectDir/bin/subset.py"
  opt_file = file(params.rng)

  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv (header: false) 
       .flatten()
       .set {samples}
  subset_object(samples)

  for_doubs = subset_object.out.samp_obj.map { it + [ opt_file ] }
  find_doublets(for_doubs)
  
  def samp_collected = subset_object.out.samp_obj.collect(){ it[0] }
  def numSamples = samp_collected.map { it.size() }
  
  pool_all(samp_collected, subset_object.out.samp_obj.collect() { it[1] }, opt_file, numSamples)
  finalize_qc_basic(pool_all.out.obj, find_doublets.out.collect(){ it[1] }, pool_all.out.numSamples)
}