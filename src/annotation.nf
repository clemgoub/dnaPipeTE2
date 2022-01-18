version = "1.4"
container_url = "dfam/tetools:${version}"

params.annotation = ""
params.annotation_out = ""

workflow annotation {
  take:
    fasta
    dfam_db
    custom_db
  main:
    repeatmasker(fasta)
    repeatmasker_extented(
      fasta,
      dfam_db.collect()
    )
    repeatmasker_custom(
      fasta,
      custom_db.collect()
    )
    channel.empty()
      .mix(repeatmasker.out.annot)
      .mix(repeatmasker_extented.out.annot)
      .mix(repeatmasker_extented.out.annot)
      .last()
      .set { rm_annot }

  emit:
    rm_annot = rm_annot
}

process repeatmasker {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.annotation_out != "") {
    publishDir "results/${params.annotation_out}", mode: 'copy'
  }

  input:

    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), val("*"), emit: annot

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
RepeatMasker \
  -pa ${task.cpus} \
  -s \
  -no_is \
  -e hmmer \
  ${fasta}
"""
}

process repeatmasker_extented {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.annotation_out != "") {
    publishDir "results/${params.annotation_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)
    tuple val(db_id), path(db)

  output:
    tuple val(file_id), val("*"), emit: annot

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
mkdir libraries
mv ${db} libraries/
RepeatMasker \
  -pa ${task.cpus} \
   -s -no_is -e hmmer -libdir libraries/
  -s \
  -no_is \
  -e hmmer \
  ${fasta}
"""
}

process repeatmasker_custom {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.annotation_out != "") {
    publishDir "results/${params.annotation_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)
    tuple val(db_id), path(db)

  output:
    tuple val(file_id), val("*"), emit: annot

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
RepeatMasker \
  -pa ${task.cpus} \
  -s \
  -no_is \
  -lib ${db} \
  ${fasta}
"""
}