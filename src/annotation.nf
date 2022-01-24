version = "1.4"
container_url = "dfam/tetools:${version}"

params.annotation = ""
params.dfam_db = false
params.custom_db = ""
params.repeatmasker_threshold = 0.0
params.annotation_out = ""

workflow annotation {
  take:
    fasta
    dfam_db
    custom_db
  main:
    if (params.custom_db != "") {
      repeatmasker_custom(
        fasta,
        custom_db.collect()
      )
      repeatmasker_custom.out.annotation
        .set{ annotation }
    } else {
      if (params.dfam_db) {
        repeatmasker_extented(
          fasta,
          dfam_db.collect()
        )
        repeatmasker_extented.out.annotation
          .set{ annotation }
      } else {
        repeatmasker(fasta)
        repeatmasker.out.annotation
          .set{ annotation }
      }
    }
    parse_repeatmasker(
      annotation
    )

  emit:
    annotation = parse_repeatmasker.out.annotation
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
    tuple val(file_id), path("*"), emit: annotation

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
    tuple val(file_id), path("*"), emit: annotation

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
  def unzip_db = "mv ${db} libraries/RepeatMaskerLib.h5";
  if (db =~ /^.*\.gz$/) {
    unzip_db = "gunzip -c ${db} > libraries/RepeatMaskerLib.h5";
  }
"""
mkdir libraries
${unzip_db}
RepeatMasker \
  -pa ${task.cpus} \
   -s -no_is -e hmmer -libdir libraries/ \
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
    tuple val(file_id), path("*"), emit: annotation

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

process parse_repeatmasker {
  container = "python:3.9"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.annotation_out != "") {
    publishDir "results/${params.annotation_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(annotation)

  output:
    tuple val(file_id), path("*"), emit: annotation

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
#! /usr/local/bin/python
with open("Trin.Clustered.out", 'r') as trinity_handle:
    line_number = 0
    trinity_out = list()
    for line in trinity_handle:
        line_number += 1
        if line_number > 3:
            line = line.split()
            # we swap the start and left column if we are revese
            if line[8] == "C":
                tmp = line[13]
                line[13] = line[11][1:-1]
                line[11] = tmp
            else:
                line[13] = line[13][1:-1]
            trinity_out_line = list()
            trinity_out_line.append(line[4])
            # size of the dnaPipeTE contig
            trinity_out_line.append(int(line[6]) + int(line[7][1:-1]))
            # percent of hit on the query
            trinity_out_line.append(float(int(line[6]) - int(line[5])) / float(int(line[6]) + int(line[7][1:-1])))
            # ET name
            trinity_out_line.append(line[9])
            # class name
            trinity_out_line.append(line[10])
            # target size
            trinity_out_line.append(int(line[12]) + int(line[13]))
            # query position
            trinity_out_line.append("["+line[11]+"-"+line[12]+"]")
            # percent of hit on the target
            trinity_out_line.append(float(int(line[12]) - int(line[11])) / float(int(line[12]) + int(line[13])))
            trinity_out_line.append(int(line[0]))
            trinity_out.append(list(trinity_out_line))
    print(str(line_number)+" line read, sorting...")
    trinity_out = sorted(trinity_out, key=lambda x: (x[0], -x[8]))
    prev_contig = ""
    print("sort done, filtering...")
    with open("one_RM_hit_per_Trinity_contigs", 'w') as output, open("Best_RM_annot_80", 'w') as output_80_80, open("Best_RM_annot_partial", 'w') as output_partial:
        line_number = 0
        line_number_80 = 0
        line_number_partial = 0
        for trinity_out_line in trinity_out:
            if trinity_out_line[0] != prev_contig:
                prev_contig = trinity_out_line[0]
                if float(trinity_out_line[2]) >= float(${params.repeatmasker_threshold}) :
                    for i in trinity_out_line[:-1]:
                        output.write(str(i)+"\\t")
                    output.write("\\n")
                    line_number += 1
                    if float(trinity_out_line[2]) >= 0.80:
                        if float(trinity_out_line[7]) >= 0.80:
                            for i in trinity_out_line[:-1]:
                                output_80_80.write(str(i)+"\\t")
                            output_80_80.write("\\n")
                            line_number_80 += 1
                        if float(trinity_out_line[7]) < 0.80:
                            for i in trinity_out_line[:-1]:
                                output_partial.write(str(i)+"\\t")
                            output_partial.write("\\n")
                            line_number_partial += 1
    print(str(line_number)+" lines in one_RM_hit_per_Trinity_contigs")
    print(str(line_number_80)+" lines in Best_RM_annot_80")
    print(str(line_number_partial)+" lines in Best_RM_annot_partial")
    print("Done")
"""
}