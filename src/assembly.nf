version = "2.13.2--ha140323_0"
container_url = "quay.io/biocontainers/trinity:${version}"

params.min_glue = 1
params.min_contig_length = 200
params.assembly_out = ""

process assembly {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.assembly_out != "") {
    publishDir "results/${params.assembly_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fastq)

  output:
    tuple val(file_id), path("*.fasta"), emit: fasta

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
  mkdir output
  Trinity \
    --seqType fq \
    --max_memory ${task.memory} \
    --left ${fastq[0]} \
    --right ${fastq[1]} \
    --CPU ${task.cpus} \
    --min_glue ${params.min_glue} \
    --min_contig_length ${params.min_contig_length} \
    --output output
mv output/Trinity.fasta Trinity.fasta
"""
  else
"""
  mkdir output
  Trinity \
    --seqType fq \
    --max_memory ${task.memory} \
    --single ${fastq} \
    --CPU ${task.cpus} \
    --min_glue ${params.min_glue} \
    --min_contig_length ${params.min_contig_length} \
    --output output
mv output/Trinity.fasta Trinity.fasta
"""
}