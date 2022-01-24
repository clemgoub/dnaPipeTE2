version = "4.1.2"
container_url = "lbmc/r-base:${version}"

params.analysis_out = ""
process analysis {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(tsv_id), path(tsv)
    tuple val(annotation_id), path(annotation)
    tuple val(cluster_id), path(cluster), path(bak_clstr)
    tuple val(transcript_id), path(fasta), path(gtf)
    
  output:
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
"""
ls -l
which R
"""
}