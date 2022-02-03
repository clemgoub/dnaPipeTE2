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
  -b 100kb \
  -i ${fastq[0]} ${fastq[1]} \
  -o sub_${fastq[0].simpleName}_R1.fastq.gz sub_${fastq[1].simpleName}_R2.fastq.gz
"""
  else
"""
rasusa \
  -b 100kb \
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
rasusa \
  --coverage ${params.coverage} \
  --genome-size ${params.genome_size} \
  -i ${fastq[0]} ${fastq[1]} \
  -O g \
  -o sub2_${fastq[0].simpleName}_R1.fastq.gz sub2_${fastq[1].simpleName}_R2.fastq.gz
cat sub2_${fastq[0].simpleName}_R1.fastq.gz \
  ${fastq_enriched[0]} \
  > enriched_sub_${fastq[0].simpleName}_R1.fastq.gz
cat sub2_${fastq[1].simpleName}_R2.fastq.gz \
  ${fastq_enriched[1]} \
  > enriched_sub_${fastq[1].simpleName}_R2.fastq.gz
rm sub2_${fastq[0].simpleName}_R1.fastq.gz sub2_${fastq[1].simpleName}_R2.fastq.gz
"""
  else
"""
rasusa \
  --coverage ${params.coverage} \
  --genome-size ${params.genome_size} \
  -i ${fastq} \
  -O g \
  -o sub2_${fastq.simpleName}.fastq.gz
cat sub2_${fastq.simpleName}.fastq.gz \
  ${fastq_enriched} \
  > enriched_sub_${fastq.simpleName}.fastq.gz
rm sub2_${fastq.simpleName}.fastq.gz
"""
}