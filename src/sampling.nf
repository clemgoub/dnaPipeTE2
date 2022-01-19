version = "0.6.0"
container_url = "lbmc/rasusa:${version}"


params.coverage = ""
params.genome_size = ""
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

  def corverage = params.coverage / 2.0;

  if (fastq.size() == 2)
"""
rasusa \
  --coverage ${params.coverage} \
  --genome-size ${params.genome_size} \
  -i ${fastq[0]} ${fastq[1]} \
  -o sub_${fastq[0].simpleName}.fastq.gz sub_${fastq[1].simpleName}.fastq.gz
"""
  else
"""
rasusa \
  --coverage ${coverage} \
  --genome-size ${params.genome_size} \
  -i ${fastq} \
  -o sub_${fastq.simpleName}.fastq.gz
"""
}

process sampling_enriched {
  container = "${container_url}"
  label "small_mem_mono_cpus"
  tag "$file_id"
  if (params.sampling_out != "") {
    publishDir "results/${params.sampling_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(file_id_enriched), path(fastq_enriched)

  output:
    tuple val(file_id), path("enriched_sub_*.fastq.gz"), emit: fastq

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
mkfifo ${fastq[0].simpleName}.pipe ${fastq[1].simpleName}.pipe
rasusa \
  --coverage ${params.coverage} \
  --genome-size ${params.genome_size} \
  -i ${fastq[0]} ${fastq[1]} \
  -O g \
  -o ${fastq[0].simpleName}.pipe ${fastq[1].simpleName}.pipe &
cat ${fastq[0].simpleName}.pipe \
  ${fastq_enriched[0]} \
  > enriched_sub_${fastq[0].simpleName}.fastq.gz &
cat ${fastq[1].simpleName}.pipe \
  ${fastq_enriched[1]} \
  > enriched_sub_${fastq[1].simpleName}.fastq.gz &
wait
"""
  else
"""
mkfifo ${fastq.simpleName}.pipe
rasusa \
  --coverage ${params.coverage} \
  --genome-size ${params.genome_size} \
  -i ${fastq} \
  -O g \
  -o ${fastq.simpleName}.pipe &
cat ${fastq.simpleName}.pipe \
  ${fastq_enriched} \
  > enriched_sub_${fastq.simpleName}.fastq.gz
"""
}