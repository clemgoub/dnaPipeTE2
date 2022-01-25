version = "0.44.0"
container_url = "lbmc/kallisto:${version}"

workflow quantification {
  take:
    fasta
    fastq
  main:
    index_fasta(fasta)
    mapping_fastq(
      index_fasta.out.index.collect(),
      fastq
    )
  emit:
    counts = mapping_fastq.out.counts
    tsv = mapping_fastq.out.tsv
}

params.index_fasta = "-k 31 --make-unique"
params.index_fasta_out = ""
process index_fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta), path(gtf)

  output:
    tuple val(file_id), path("*.index*"), emit: index
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
"""
kallisto index ${params.index_fasta} -i ${fasta.baseName}.index ${fasta} \
2> ${fasta.baseName}_kallisto_index_report.txt
"""
}

params.mapping_fastq = "--bias --bootstrap-samples 100"
params.mapping_fastq_out = ""
process mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.mapping_fastq_out != "") {
    publishDir "results/${params.mapping_fastq_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("${file_prefix}"), emit: counts
  tuple val(file_id), path("${file_prefix}/*.tsv"), emit: tsv
  tuple val(file_id), path("*_report.txt"), emit: report

  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }

  if (reads.size() == 2)
  """
  mkdir ${file_prefix}
  kallisto quant -i ${index} -t ${task.cpus} \
  ${params.mapping_fastq} -o ${file_prefix} \
  ${reads[0]} ${reads[1]} &> ${file_prefix}_kallisto_mapping_report.txt
  grep "n_processed" ${file_prefix}/run_info.json \
    | sed 's/.*"n_processed":/read_processed\t0\t/' \
    | sed 's/,/\t0/' \
    >> ${file_prefix}/abundance.tsv
  """
  else
  """
  mkdir ${file_prefix}
  kallisto quant -i ${index} -t ${task.cpus} --single \
  ${params.mapping_fastq} -o ${file_prefix} \
  ${reads[0]} &> ${file_prefix}_kallisto_mapping_report.txt
  grep "n_processed" ${file_prefix}/run_info.json \
    | sed 's/.*"n_processed":/read_processed\t0\t/' \
    | sed 's/,/\t0/' \
    >> ${file_prefix}/abundance.tsv
  """
}
