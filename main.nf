#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process gather_matrices {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/1_gathered_objects", mode: 'copy', saveAs: {filename -> "${samp}.${filename}"}

  input:
  val(samp)
  path(cr_gene), stageAs: "filteredSS"
  path(cr_velo), stageAs: "filteredVelo"
  path(cb_h5)

  output:
  tuple val(samp), path("*gene_velo_cellbender.filtered.h5ad"), emit: obj

  script:
  """
  python ${baseDir}/bin/gather_matrices.py --cr_gene ${cr_gene} --cr_velo ${cr_velo} --cb_h5 ${cb_h5} --gather_mode ${params.gather_mode}
  """
}

process run_qc {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/2_qc_objects", pattern: '*.h5ad', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/2_qc_plots_individual", pattern: '*.png', mode: 'copy'

  input:
  tuple val(samp), path(gath_out)

  output:
  tuple val(samp), path("*_postqc.h5ad"), path("*-scr"), path("*.csv"), emit: samp_obj
  path("*.png")

  script:
  """
  python ${baseDir}/bin/qc.py --sample_id ${samp} --metrics_csv ${params.metrics_csv} --celltypist ${params.celltypist_model} --qc_mode ${params.qc_mode} --gath_obj ${gath_out}
  """
}

process subset_object {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/2_qc_objects", pattern: '*.h5ad', mode: 'copy', saveAs: {filename -> "${samp}_${filename}"}

  input:
  val(samp)

  output:
  tuple val(samp), path("subsetted.h5ad"), path("*-scr"), emit: samp_obj

  script:
  """
  python ${baseDir}/bin/subset.py --sample_id ${samp} --cr_prefix ${params.cr_prefix} --limits_csv ${params.limits_csv}
  """
}

process find_doublets {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/3_doublet_scores", pattern: '*.csv', mode: 'copy'

  input:
  tuple val(samp), path(qc_out), path(scr_bool), path(qc_thres)

  output:
  tuple val(samp), path("*_scrublet.csv")

  when:
  scr_bool.name.endsWith('yes-scr')

  script:
  def filter_column = ""
  if (params.qc_mode == "multires") {
    filter_column = "cluster_passed_qc"
  } else {
    filter_column = "good_qc_cluster_mito80"  
  }
  
  """
  python ${baseDir}/bin/flag_doublet.py --filter ${filter_column} --samp ${samp} --input ${qc_out}
  """

}


process pool_all {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.csv', mode: 'copy'

  memory = { 
        def numSamples = samp.size()  // Use the size of the input list
        return 2.GB * numSamples * task.attempt
    }

  input:
  val(samp)
  path(qc_out)
  path(ranges_out)

  output:
  path("scautoqc_pooled0.h5ad"), emit: obj
  path("qc_thresholds.csv")

  script:
  """
  python ${baseDir}/bin/pool_all.py --samples ${samp.join(",")} --objects ${qc_out.join(",")} --ranges ${ranges_out.join(",")}
  """
}

process finalize_qc {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.h5ad', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/5_qc_plots_overall", pattern: '*.png', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.csv', mode: 'copy'

  memory = { 5.GB * (pool_out.size() / (1024 * 1024 * 1024)) * task.attempt }

  input:
  path(pool_out)
  path(scr_out)

  output:
  path("scautoqc_pooled_filtered.h5ad"), emit: obj
  path("scautoqc_pooled.h5ad")
  path("*.png")
  path("sample_passqc_df.csv")

  script:
  """
  export BASE_DIR=${baseDir}
  python ${baseDir}/bin/finalize_qc.py --obj ${pool_out} --scr ${scr_out.join(",")} --meta ${params.metadata} --qc_mode ${params.qc_mode}
  """
}

process finalize_qc_basic {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.h5ad', mode: 'copy'

  memory = { 5.GB * (pool_out.size() / (1024 * 1024 * 1024)) * task.attempt }
  
  input:
  path(pool_out)
  path(scr_out)

  output:
  path("scautoqc_pooled_basic.h5ad"), emit: obj

  script:
  """
  export BASE_DIR=${baseDir}
  python ${baseDir}/bin/finalize_qc_basic.py --obj ${pool_out} --scr ${scr_out.join(",")} --meta ${params.metadata}
  """
}

process integrate {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}", pattern: '*.h5ad', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/5_qc_plots_overall", pattern: '*.png', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/models", pattern: '*.pkl', mode: 'copy'

  input:
  path(qc2_out)

  output:
  path("scautoqc_integrated.h5ad")
  path("*.png")
  path("*.pkl")

  script:
  """
  python ${baseDir}/bin/integration.py --obj ${qc2_out} --batch ${params.batch_key} --n_top_genes ${params.n_top_genes} --from_scautoqc ${params.from_scautoqc}
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
                cell_or_nuclei: 'cell' // default value for legacy format
            }
    }
}


workflow all {
  def samples = createInputChannels("${params.SAMPLEFILE}")

  gather_matrices(samples.samp, samples.cr_gene, samples.cr_velo, samples.cb_h5)
  run_qc(gather_matrices.out.obj)
  find_doublets(run_qc.out.samp_obj)
  pool_all(run_qc.out.samp_obj.collect(){ it[0] }, run_qc.out.samp_obj.collect() { it[1] }, run_qc.out.samp_obj.collect() { it[3] })
  finalize_qc(pool_all.out.obj, find_doublets.out.collect(){ it[1] })
  integrate(finalize_qc.out.obj, params.batch_key)
}

workflow only_qc {
  def samples = createInputChannels("${params.SAMPLEFILE}")

  gather_matrices(samples.samp, samples.cr_gene, samples.cr_velo, samples.cb_h5)
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
  pool_all(run_qc.out.samp_obj.collect(){ it[0] }, run_qc.out.samp_obj.collect() { it[1] }, run_qc.out.samp_obj.collect() { it[3] })
  finalize_qc(pool_all.out, find_doublets.out.collect(){ it[1] })
  integrate(finalize_qc.out.obj, params.batch_key)
}

workflow until_integrate {
  def samples = createInputChannels("${params.SAMPLEFILE}")

  gather_matrices(samples.samp, samples.cr_gene, samples.cr_velo, samples.cb_h5)
  run_qc(gather_matrices.out.obj)
  find_doublets(run_qc.out.samp_obj)
  pool_all(run_qc.out.samp_obj.collect(){ it[0] }, run_qc.out.samp_obj.collect() { it[1] }, run_qc.out.samp_obj.collect() { it[3] })
  finalize_qc(pool_all.out.obj, find_doublets.out.collect(){ it[1] })
}

workflow only_integrate {
  obj_scvi = Channel.fromPath("${params.path_for_scvi}")
  integrate(obj_scvi, params.batch_key)
}

workflow subset {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv (header: false) 
       .flatten()
       .set {samples}
  subset_object(samples)
  find_doublets(subset_object.out.samp_obj)
  pool_all(run_qc.out.samp_obj.collect(){ it[0] }, run_qc.out.samp_obj.collect() { it[1] }, run_qc.out.samp_obj.collect() { it[3] })
  finalize_qc_basic(pool_all.out, find_doublets.out.collect(){ it[1] })
}