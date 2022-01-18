version = "4.8.1--h2e03b76_5"
container_url = "quay.io/biocontainers/cd-hit:${version}"

params.clustering = "-d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -bak 1"
params.clustering_out = ""

process clustering {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.clustering_out != "") {
    publishDir "results/${params.clustering_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), path("Trin.Clustered.clstr"), path("Trin.Clustered.bak.clstr"), emit: cluster
    tuple val(file_id), path("Trin.Clustered"), emit: fasta 

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }

"""
cd-hit-est \
  -i ${fasta} \
  -o Trin.Clustered \
  ${params.clustering} \
  -T ${task.cpus} \
  -M 0 
"""
}