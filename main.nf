#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process gather_matrices {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/gathered_h5ad", mode: 'copy', saveAs: {filename -> "${samp}.${filename}"}

  input:
  val(samp)
  path(cr_gene), stageAs: "filteredSS"
  path(cr_velo), stageAs: "filteredVelo"
  path(cb_h5)

  output:
  tuple val(samp), path("*gene_velo_cellbender.filtered.h5ad"), emit: obj

  script:
  """
  python ${baseDir}/bin/gather_matrices.py --cr_gene ${cr_gene} --cr_velo ${cr_velo} --cb_h5 ${cb_h5}
  """
}

process run_qc {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/postqc_h5ad", pattern: '*.h5ad', mode: 'copy', saveAs: {filename -> "${samp}_${filename}"}
  // publishDir "${launchDir}/scautoqc-results-${params.project_tag}/qc_plots/${samp}", pattern: '*.png', mode: 'copy'

  input:
  tuple val(samp), path(gath_out)

  output:
  tuple val(samp), path("gene_velo_cellbender.post_qc.h5ad"), emit: samp_obj

  script:
  """
  python ${baseDir}/bin/qc.py --clst_res ${params.cluster_res} --min_frac ${params.min_frac} --models ${params.models} --sample_id ${samp} --out_path ${gath_out}
  """
}

process find_doublets {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/scrublet_out", pattern: '*.csv', mode: 'copy', saveAs: {filename -> "${samp}_${filename}"}

  input:
  tuple val(samp), path(qc_out)

  output:
  tuple val(samp), path("gene_velo_cellbender.good_qc_cluster_mito80.scrublet.csv")

  script:
  """
  python ${baseDir}/bin/flag_doublet.py --filter good_qc_cluster_mito80 --samp ${samp} --input ${qc_out}
  """
}


process pool_all {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/pooled_h5ad", pattern: '*.h5ad', mode: 'copy'

  input:
  val(samp)
  val(qc_out)

  output:
  path("pooled.gene_velo_cellbender.post_qc.h5ad")

  script:
  """
  python ${baseDir}/bin/pool_all.py --samples ${samp} --objects ${qc_out}
  """
}

process add_metadata {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/pooled_h5ad", pattern: '*.h5ad', mode: 'copy'
  // publishDir "${launchDir}/scautoqc-results-${params.project_tag}/qc_plots", pattern: '*.png', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/", pattern: '*.csv', mode: 'copy'

  input:
  path(pool_out)
  val(scr_out)
  val(meta_path)

  output:
  path("pooled.gene_cellbender.good_qc_cluster_mito80.doublet_flagged.h5ad")

  script:
  """
  python ${baseDir}/bin/add_scrublet_meta.py --obj ${pool_out} --scr ${scr_out} --meta ${meta_path}
  """
}

process integrate {

  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/pooled_h5ad", pattern: '*.h5ad', mode: 'copy'
  publishDir "${launchDir}/scautoqc-results-${params.project_tag}/qc_plots", pattern: '*.png', mode: 'copy'

  input:
  path(qc2_out)

  output:
  path("pooled_healthy.gene_cellbender.good_qc_cluster_mito80.stringent_doublet_removed.hvg7500_noCC.scvi_output.lv20_batch256.h5ad")

  script:
  """
  python ${baseDir}/bin/integration.py --obj ${qc2_out}
  """
}

workflow {

  Channel.fromPath("${params.SAMPLEFILE}")
       .splitCsv (header: false) 
       .flatten()
   .multiMap { it ->
           samp: it
           cr_gene: "${params.ss_prefix}/${it}/output/Gene/filtered/"
           cr_velo: "${params.ss_prefix}/${it}/output/Velocyto/filtered/"
           cb_h5:   "${params.cb_prefix}/${it}/cellbender_out_filtered.h5"           }
       .set {samples}
  gather_matrices(samples.samp, samples.cr_gene, samples.cr_velo, samples.cb_h5)
  run_qc(gather_matrices.out.obj)
  find_doublets(run_qc.out.samp_obj)
  pool_all(run_qc.out.samp_obj.collect{ it[0] }.map { it.join(',') },run_qc.out.samp_obj.collect{ it[1] }.map { it.join(',') })
  add_metadata(pool_all.out, find_doublets.out.collect{ it[1] }.map { it.join(',') }, params.metadata)
  integrate(add_metadata.out)
}

