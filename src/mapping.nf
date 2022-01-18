version = "2.2.1"
container_url = "lbmc/hisat2:${version}"

params.mapping_out = ""

workflow extract {
  take:
    fasta
    fastq
  main:
    index_fasta(fasta)
    extract_mapping_fastq(
      index_fasta.out.index,
      fastq
    )

  emit:
    fastq = extract_mapping_fastq.out.fastq
}

workflow mapping {
  take:
    fasta
    fastq
  main:
    index_fasta(fasta)
    mapping_fastq(
      index_fasta.out.index,
      fastq
    )

  emit:
    bam = mapping_fastq.out.bam
}

params.index_fasta_out = ""
process index_fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.index_fasta_out != "") {
    publishDir "results/${params.index_fasta_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), path("*.ht2*"), emit: index
    tuple val(file_id), path("*_report.txt"), emit: report

  script:
"""
hisat2-build -p ${task.cpus} \
  ${fasta} \
  ${fasta.simpleName} &> \
  ${fasta.simpleName}_hisat2_index_report.txt

if grep -q "Error" ${fasta.simpleName}_bowtie2_index_report.txt; then
  exit 1
fi
"""
}

params.mapping_fastq = ""
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
  tuple val(file_id), path("*.bam"), emit: bam
  path "*_report_tmp.txt", emit: report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.ht2.*/) {
        index_id = ( index_file =~ /(.*)\.1\.ht2.*/)[0][1]
    }
  }
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

  if (reads.size() == 2)
  """
  hisat2 ${params.mapping_fastq} \
    -p ${task.cpus} \
    -x ${index_id} \
    -1 ${reads[0]} \
    -2 ${reads[1]} 2> \
    ${file_prefix}_ht2_mapping_report_tmp.txt \
    | samtools view -@ ${task.cpus} -b - \
    | samtools sort -@ ${task.cpus} -o ${file_prefix}.bam -

  if grep -q "Error" ${file_prefix}_ht2_mapping_report_tmp.txt; then
    exit 1
  fi
  """
  else
  """
  hisat2 ${params.mapping_fastq} \
    -p ${task.cpus} \
    -x ${index_id} \
    -U ${reads} 2> \
    ${file_prefix}_ht2_mapping_report_tmp.txt \
    | samtools view -@ ${task.cpus} -b - \
    | samtools sort -@ ${task.cpus} -o ${file_prefix}.bam -

  if grep -q "Error" ${file_prefix}_ht2_mapping_report_tmp.txt; then
    exit 1
  fi
  """
}

params.extract_mapping_fastq_out = ""
process extract_mapping_fastq {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.mapping_fastq_out != "") {
    publishDir "results/${params.extract_mapping_fastq_out}", mode: 'copy'
  }

  input:
  tuple val(index_id), path(index)
  tuple val(file_id), path(reads)

  output:
  tuple val(file_id), path("selected_*.fastq.gz"), emit: fastq 
  path "*_report_tmp.txt", emit: report

  script:
  index_id = index[0]
  for (index_file in index) {
    if (index_file =~ /.*\.1\.ht2.*/) {
        index_id = ( index_file =~ /(.*)\.1\.ht2.*/)[0][1]
    }
  }
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

  if (reads.size() == 2)
  """
  mkfifo ${reads[0].simpleName}.pipe ${reads[1].simpleName}.pipe
  hisat2 ${params.mapping_fastq} \
    -p ${task.cpus} \
    -x ${index_id} \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    2> ${file_prefix}_ht2_mapping_report_tmp.txt \
    | samtools view -b -F 4 - \
    | samtools collate -u -O -@ ${task.cpus} - \
    | samtools fastq -1 ${reads[0].simpleName}.pipe -2 ${reads[1].simpleName}.pipe -0 /dev/null -s /dev/null -n \
    2> ${file_prefix}_samtool_report_tmp.txt &
  gzip -c ${reads[0].simpleName}.pipe > selected_${reads[0].simpleName}.fastq.gz &
  gzip -c ${reads[1].simpleName}.pipe > selected_${reads[1].simpleName}.fastq.gz &
  wait

  if grep -q "Error" ${file_prefix}_ht2_mapping_report_tmp.txt; then
    exit 1
  fi
  """
  else
  """
  mkfifo ${reads.simpleName}.pipe
  hisat2 ${params.mapping_fastq} \
    -p ${task.cpus} \
    -x ${index_id} \
    -U ${reads} \
    2> ${file_prefix}_ht2_mapping_report_tmp.txt \
    | samtools view -b -F 4 - \
    | samtools collate -u -O -@ ${task.cpus} - \
    | samtools fastq -0 ${reads.simpleName}.pipe -n 
    2> ${file_prefix}_samtool_report_tmp.txt &
  gzip -c ${reads[0].simpleName}.pipe > selected_${reads[0].simpleName}.fastq.gz &
  wait

  if grep -q "Error" ${file_prefix}_ht2_mapping_report_tmp.txt; then
    exit 1
  fi
  """
}