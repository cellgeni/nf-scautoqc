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
  python ${baseDir}/bin/gather_matrices.py --cr_gene ${cr_gene} --cr_velo ${cr_velo} --cb_h5 ${cb_h5} --ss_out ${params.ss_out}
  """
}

process run_qc {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/2_qc_objects", pattern: '*.h5ad', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/2_qc_plots_individual", pattern: '*.png', mode: 'copy'

  input:
  tuple val(samp), path(gath_out)

  output:
  tuple val(samp), path("*_postqc.h5ad"), path("*-scr"), emit: samp_obj
  path("*.png")

  script:
  """
  python ${baseDir}/bin/qc.py --clst_res ${params.cluster_res} --min_frac ${params.min_frac} --models ${params.models} --sample_id ${samp} --out_path ${gath_out} --ss_out ${params.ss_out}
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
  tuple val(samp), path(qc_out), path(scr_bool)

  output:
  tuple val(samp), path("*_scrublet.csv")

  when:
  scr_bool.name.endsWith('yes-scr')

  script:
  """
  python ${baseDir}/bin/flag_doublet.py --filter good_qc_cluster_mito80 --samp ${samp} --input ${qc_out}
  """
}


process pool_all {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.h5ad', mode: 'copy'

  memory = { 
        def numSamples = samp.size()  // Use the size of the input list
        return 2.GB * numSamples * task.attempt
    }

  input:
  val(samp)
  path(qc_out)

  output:
  path("scautoqc_pooled.h5ad")

  script:
  """
  python ${baseDir}/bin/pool_all.py --samples ${samp.join(",")} --objects ${qc_out.join(",")} --ss_out ${params.ss_out}
  """
}

process add_metadata {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.h5ad', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/5_qc_plots_overall", pattern: '*.png', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.csv', mode: 'copy'

  memory = { 3.6.GB * (pool_out.size() / (1024 * 1024 * 1024)) * task.attempt }

  input:
  path(pool_out)
  path(scr_out)
  path(meta_path)

  output:
  path("scautoqc_pooled_doubletflagged_metaadded.h5ad"), emit: obj
  path("*.png")
  path("sample_passqc_df.csv")

  script:
  """
  export BASE_DIR=${baseDir}
  python ${baseDir}/bin/add_scrublet_meta.py --obj ${pool_out} --scr ${scr_out.join(",")} --meta ${meta_path}
  """
}

process add_metadata_basic {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.h5ad', mode: 'copy'

  memory = { 3.6.GB * (pool_out.size() / (1024 * 1024 * 1024)) * task.attempt }
  
  input:
  path(pool_out)
  val(scr_out)
  val(meta_path)

  output:
  path("scautoqc_pooled_doubletflagged_metaadded_basic.h5ad"), emit: obj

  script:
  """
  export BASE_DIR=${baseDir}
  python ${baseDir}/bin/add_scrublet_meta_basic.py --obj ${pool_out} --scr ${scr_out} --meta ${meta_path}
  """
}

process integrate {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}", pattern: '*.h5ad', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/5_qc_plots_overall", pattern: '*.png', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/models", pattern: '*.pkl', mode: 'copy'

  input:
  path(qc2_out)
  val(batch_key)

  output:
  path("scautoqc_integrated.h5ad")
  path("*.png")
  path("*.pkl")

  script:
  """
  python ${baseDir}/bin/integration.py --obj ${qc2_out} --batch ${batch_key}
  """
}

workflow all {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv(header: false) 
       .flatten()
       .map { it ->
          def resolvedPath = "readlink -f ${params.ss_prefix}/${it}".execute().text.trim()
          [it, resolvedPath]
      }
       .multiMap { it, resolvedPath -> 
           samp: it
           cr_gene: "${params.ss_prefix}" == "" ? "${resolvedPath}/filtered_feature_bc_matrix.h5" : "${resolvedPath}/output/${params.ss_out}/filtered/"
           cr_velo: "${params.ss_prefix}" == "" ? [] : "${resolvedPath}/output/Velocyto/filtered/"
           cb_h5: "${params.cb_prefix}" == "" ? [] : "${params.cb_prefix}/${it}/cellbender_out_filtered.h5"
       }
       .set {samples}
  gather_matrices(samples.samp, samples.cr_gene, samples.cr_velo, samples.cb_h5)
  run_qc(gather_matrices.out.obj)
  find_doublets(run_qc.out.samp_obj)
  pool_all(run_qc.out.samp_obj.collect(){ it[0] },run_qc.out.samp_obj.collect() { it[1] })
  add_metadata(pool_all.out, find_doublets.out.collect(){ it[1] }, params.metadata)
  integrate(add_metadata.out.obj, params.batch_key)
}

workflow only_qc {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv(header: false) 
       .flatten()
       .map { it ->
          def resolvedPath = "readlink -f ${params.ss_prefix}/${it}".execute().text.trim()
          [it, resolvedPath]
      }
       .multiMap { it, resolvedPath -> 
           samp: it
           cr_gene: "${params.ss_prefix}" == "" ? "${resolvedPath}/filtered_feature_bc_matrix.h5" : "${resolvedPath}/output/${params.ss_out}/filtered/"
           cr_velo: "${params.ss_prefix}" == "" ? [] : "${resolvedPath}/output/Velocyto/filtered/"
           cb_h5: "${params.cb_prefix}" == "" ? [] : "${params.cb_prefix}/${it}/cellbender_out_filtered.h5"
       }
       .set {samples}
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
  pool_all(run_qc.out.samp_obj.collect(){ it[0] },run_qc.out.samp_obj.collect() { it[1] })
  add_metadata(pool_all.out, find_doublets.out.collect(){ it[1] }, params.metadata)
  integrate(add_metadata.out.obj, params.batch_key)
}

workflow until_integrate {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv(header: false) 
       .flatten()
       .map { it ->
          def resolvedPath = "readlink -f ${params.ss_prefix}/${it}".execute().text.trim()
          [it, resolvedPath]
      }
       .multiMap { it, resolvedPath -> 
           samp: it
           cr_gene: "${params.ss_prefix}" == "" ? "${resolvedPath}/filtered_feature_bc_matrix.h5" : "${resolvedPath}/output/${params.ss_out}/filtered/"
           cr_velo: "${params.ss_prefix}" == "" ? [] : "${resolvedPath}/output/Velocyto/filtered/"
           cb_h5: "${params.cb_prefix}" == "" ? [] : "${params.cb_prefix}/${it}/cellbender_out_filtered.h5"
       }
       .set {samples}
  gather_matrices(samples.samp, samples.cr_gene, samples.cr_velo, samples.cb_h5)
  run_qc(gather_matrices.out.obj)
  find_doublets(run_qc.out.samp_obj)
  pool_all(run_qc.out.samp_obj.collect(){ it[0] },run_qc.out.samp_obj.collect() { it[1] })
  add_metadata(pool_all.out, find_doublets.out.collect(){ it[1] }, params.metadata)
}

workflow only_integrate {
  integrate(params.path_for_scvi, params.batch_key)
}

workflow subset {
  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv (header: false) 
       .flatten()
       .set {samples}
  subset_object(samples)
  find_doublets(subset_object.out.samp_obj)
  pool_all(subset_object.out.samp_obj.collect(){ it[0] },subset_object.out.samp_obj.collect() { it[1] })
  add_metadata_basic(pool_all.out, find_doublets.out.collect(){ it[1] }, params.metadata)
}