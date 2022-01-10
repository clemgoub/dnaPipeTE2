version = "0.6.0"
container_url = "lbmc/rasusa:${version}"


params.coverage = ""
params.genome_size = ""
params.seed = ""
params.sampling_out = ""

process sampling {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "$file_id"
  if (params.sampling_out != "") {
    publishDir "results/${params.sampling_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fastq)

  output:
    tuple val(file_id), path("sub_*.fastq.gz"), emit: fastq

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

  if (fastq.size() == 2)
"""
rasusa \
  --seed ${params.seed} \
  --coverage ${params.coverage} \
  --genome-size ${params.genome_size} \
  -i ${fastq[0]} ${fastq[1]} \
  -o sub_${fastq[0].simpleName}.fastq.gz sub_${fastq[1].simpleName}.fastq.gz
"""
  else
"""
rasusa \
  --seed ${params.seed} \
  --coverage ${params.coverage} \
  --genome-size ${params.genome_size} \
  -i ${fastq} \
  -o sub_${fastq.simpleName}.fastq.gz
"""
}