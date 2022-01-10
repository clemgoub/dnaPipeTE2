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
    tuple val(group_id), val(file_id), path(fastq)

  output:
    tuple val(group_id), path("trinity_output_${file_prefix}"), emit: folder

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break;
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break;
    default:
      file_prefix = file_id
    break;
  };
  def memory = "${task.memory}" - ~/\s*GB/

  if (fastq.size() == 2)
"""
  mkdir trinity_output_${file_prefix}
  Trinity \
    --no_run_chrysalis \
    --seqType fq \
    --max_memory ${memory}G \
    --left ${fastq[0]} \
    --right ${fastq[1]} \
    --CPU ${task.cpus} \
    --min_glue ${params.min_glue} \
    --min_contig_length ${params.min_contig_length} \
    --output trinity_output_${file_prefix}
"""
  else
"""
  mkdir trinity_output_${file_prefix}
  Trinity \
    --no_run_chrysalis \
    --seqType fq \
    --max_memory ${memory}G \
    --single ${fastq} \
    --CPU ${task.cpus} \
    --min_glue ${params.min_glue} \
    --min_contig_length ${params.min_contig_length} \
    --output trinity_output_${file_prefix}
"""
}

process merge {
  container = "lbmc/bioawk:1.0"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.assembly_out != "") {
    publishDir "results/${params.assembly_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(folders)

  output:
    tuple val(file_id), path("trinity_output_merge"), emit: folder

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break;
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break;
    default:
      file_prefix = file_id
    break;
  };
  def both_fas = [];
  def inchworm_ds_fa = [];
  def jellyfish_kmers_25_asm_fa = [];
  def jellyfish_kmers_25_asm_fa_histo = [];
  def r1_norm_fq = [];
  def r2_norm_fq = [];
  def left_norm_fq = [];
  def right_norm_fq = [];
  for (folder in folders) {
    both_fas.add(folder + "/both.fa");
    inchworm_ds_fa.add(folder + "/inchworm.DS.fa");
    jellyfish_kmers_25_asm_fa.add(folder + "/jellyfish.kmers.25.asm.fa");
    jellyfish_kmers_25_asm_fa_histo.add(folder + "/jellyfish.kmers.25.asm.fa.histo");
    r1_norm_fq.add(folder + "/insilico_read_normalization/sub_" + file_prefix + "_R1.fastq.gz.normalized_K25_maxC200_minC1_maxCV10000.fq");
    r2_norm_fq.add(folder + "/insilico_read_normalization/sub_" + file_prefix + "_R2.fastq.gz.normalized_K25_maxC200_minC1_maxCV10000.fq");
    left_norm_fq.add(folder + "/insilico_read_normalization/left.norm.fq");
    right_norm_fq.add(folder + "/insilico_read_normalization/right.norm.fq");
  }
"""
  cp -RL ${folders[0]} trinity_output_merge

  cat ${both_fas.join(" ")} \
    | awk '/^>/ { f = !(\$0 in a); a[\$0] } f==1 { print \$0 }' \
    > trinity_output_merge/both.fa
  
  grep ">" trinity_output_merge/both.fa \
    | wc -l \
    > trinity_output_merge/both.fa.read_count

  cat ${inchworm_ds_fa.join(" ")} \
    | awk 'BEGIN{seq_number=0};  !/^>/ {print \$0}; /^>/ { seq_number++; sub(">a[0-9]+;", ";", \$0); print ">a" seq_number \$0}' \
    > trinity_output_merge/inchworm.ds.fa
    
  # yolo kmer count with awk
  cat ${jellyfish_kmers_25_asm_fa.join(" ")} \
    | bioawk -c fastx '{kmer[\$seq]+=\$name} END{for (k in kmer) {print ">" k; print kmer[k]}}' \
    > trinity_output_merge/jellyfish.kmers.25.asm.fa

  cat ${jellyfish_kmers_25_asm_fa_histo.join(" ")} \
    | awk 'BEGIN{kmer_count = 0}; {kmer_count += \$2}; END{print kmer_count}' \
    > trinity_output_merge/inchworm.kmer_count

  cat ${r1_norm_fq.join(" ")} > trinity_output_merge/insilico_read_normalization/sub_${file_prefix}_R1.fastq.gz.normalized_K25_maxC200_minC1_maxCV10000.fq
  cat ${r2_norm_fq.join(" ")} > trinity_output_merge/insilico_read_normalization/sub_${file_prefix}_R2.fastq.gz.normalized_K25_maxC200_minC1_maxCV10000.fq
  cat ${left_norm_fq.join(" ")} > trinity_output_merge/insilico_read_normalization/left.norm.fq
  cat ${right_norm_fq.join(" ")} > trinity_output_merge/insilico_read_normalization/right.norm.fq

  cat ${jellyfish_kmers_25_asm_fa_histo.join(" ")} \
    | awk '{a[\$1] += \$2}; END{for (kmer in a){print kmer " " a[kmer]}}' \
    > trinity_output_merge/jellyfish.kmers.25.asm.fa.histo
"""
}

process final_assembly {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.assembly_out != "") {
    publishDir "results/${params.assembly_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(folder_id), path(folder)

  output:
    tuple val(file_id), path("trinity_output/"), emit: folder
    tuple val(file_id), path("trinity_output.Trinity.fasta"), emit: fasta
    tuple val(file_id), path("trinity_output/salmon_outdir/quant.sf"), emit: quant

  script:

  switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break;
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break;
    default:
      file_prefix = file_id
    break;
  };
  def memory = "${task.memory}" - ~/\s*GB/

  if (fastq.size() == 2)
"""
  cp -RL ${folder} trinity_output
  Trinity \
    --seqType fq \
    --max_memory ${memory}G \
    --left ${fastq[0]} \
    --right ${fastq[1]} \
    --CPU ${task.cpus} \
    --min_glue ${params.min_glue} \
    --min_contig_length ${params.min_contig_length} \
    --output trinity_output
"""
  else
"""
  cp -RL ${folder} trinity_output
  Trinity \
    --no_run_chrysalis \
    --seqType fq \
    --max_memory ${memory}G \
    --single ${fastq} \
    --CPU ${task.cpus} \
    --min_glue ${params.min_glue} \
    --min_contig_length ${params.min_contig_length} \
    --output trinity_output
"""
}