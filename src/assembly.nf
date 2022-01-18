version = "2.13.2--ha140323_0"
container_url = "quay.io/biocontainers/trinity:${version}"

params.sample = 3
params.min_glue = 1
params.min_contig_length = 200
params.assembly_out = ""

params.coverage = ""
params.genome_size = ""
include {
  sampling;
  sampling_enriched
} from './sampling.nf' addParams(
  coverage: params.coverage,
  genome_size: params.size
);

include {
  extract
} from './mapping.nf';

workflow assembly {
  take:
    fastq
  main:
    sampling(fastq.first())
    enrichment.recurse(fastq.first(), sampling.out.fastq).times(params.sample)
    complete_assembly(enrichment.out.fastq_enriched.last())

  emit:
    fasta = complete_assembly.out.fasta
}

workflow enrichment {
  take:
    fastq
    fastq_enriched
  main:
    sampling_enriched(fastq, fastq_enriched)
    pre_assembly(sampling_enriched.out.fastq)
    extract(
      pre_assembly.out.fasta,
      fastq
    )
  emit:
  fastq = fastq
  fastq_enriched = extract.out.fastq
}

process pre_assembly {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.assembly_out != "") {
    publishDir "results/${params.assembly_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fastq)

  output:
    tuple val(file_id), path("trinity_output_${file_prefix}/inchworm.DS.fa"), emit: fasta

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

process reads_enrichment {
  container = "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.assembly_out != "") {
    publishDir "results/${params.assembly_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fastq)
    tuple val(names_file_id), path(read_names)

  output:
    tuple val(file_id), path("enriched_*"), emit: fastq

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
  if (fastq.size() == 2)
"""
# generate a list of read names
awk '{print \$2}' ${read_names} \
  | grep "/1" \
  | sed 's|>||' \
  | sed 's|/1||' \
  > name.lst

# we remove every thing after the first space in the read names
zcat ${fastq[0]} \
  | sed -E 's|@([^ ]*).*|@\\1|' \
  | gzip -c \
  | seqtk subseq - name.lst \
  | gzip -c \
  > enriched_${fastq[0]}

zcat ${fastq[1]} \
  | sed -E 's|@([^ ]*).*|@\\1|' \
  | gzip -c \
  | seqtk subseq - name.lst \
  | gzip -c \
  > enriched_${fastq[1]}
"""
  else
"""
awk '{print \$2}' ${read_names} \
  | grep "/1" \
  | sed 's|>||' \
  | sed 's|/1||' \
  > name.lst

# we remove every thing after the first space in the read names
zcat ${fastq} \
  | sed -E 's|@([^ ]*).*|@\\1|' \
  | gzip -c \
  | seqtk subseq - name.lst \
  | gzip -c \
  > enriched_${fastq}
"""
}

process complete_assembly {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.assembly_out != "") {
    publishDir "results/${params.assembly_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fastq)

  output:
    tuple val(file_id), path("trinity_output_${file_prefix}/"), emit: folder
    tuple val(file_id), path("trinity_output_${file_prefix}.Trinity.fasta"), emit: fasta
    tuple val(file_id), path("trinity_output_${file_prefix}.Trinity.fasta.gene_trans_map"), emit: gene_map 
    tuple val(file_id), path("trinity_output_${file_prefix}/salmon_outdir/quant.sf"), emit: quant

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
    --seqType fq \
    --max_memory ${memory}G \
    --single ${fastq} \
    --CPU ${task.cpus} \
    --min_glue ${params.min_glue} \
    --min_contig_length ${params.min_contig_length} \
    --output trinity_output_${file_prefix}
"""
}