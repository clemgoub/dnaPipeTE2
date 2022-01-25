version = "4.1.2"
container_url = "lbmc/r-base:${version}"

process analysis {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$tsv_id"
  publishDir "results/", mode: 'copy'

  input:
    tuple val(tsv_id), path(tsv)
    tuple val(annotation_id), path(annotation)
    tuple val(cluster_id), path(cluster), path(bak_clstr)
    tuple val(transcript_id), path(fasta), path(gtf)

  output:
    tuple val(tsv_id), path("files/"), emit: report

  script:
"""
mkdir -p files
mv ${tsv} ${annotation} ${cluster} ${bak_clstr} ${fasta} ${gtf} files/
which R
"""
}