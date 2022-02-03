version = "4.8.1--h2e03b76_5"
container_url = "quay.io/biocontainers/cd-hit:${version}"

workflow clustering {
  take:
    fasta

  main:
    cluster_fasta(fasta)
    cluster_member(cluster_fasta.out.cluster)
    extract_cluster_member(
      fasta.collect(),
      cluster_member.out.names
        .flatten()
        .map {it -> ["multiple", it]}
        .mix(
          cluster_member.out.mono_cluster
          .map {it -> ["mono", it]}
        )
    )
    align_cluster(
      extract_cluster_member.out.fasta
        .filter {it[0] == "multiple"}
    )
    consensus_cluster(
      align_cluster.out.align
    )
    merge_cluster_sequences(
      fasta.collect(),
      consensus_cluster.out.fasta
        .mix(
          extract_cluster_member.out.fasta
            .filter {it[0] == "mono"}
        )
        .map{ it -> it[1] }
        .collect()
    )

  emit:
    fasta = merge_cluster_sequences.out.fasta
    cluster = cluster_fasta.out.cluster
}

params.clustering = "-d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -bak 1"
params.clustering_out = ""

process cluster_fasta {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.clustering_out != "") {
    publishDir "results/${params.clustering_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), path("Trin.Clustered.clstr"), path("Trin.Clustered.bak.clstr"), emit: cluster
    tuple val(file_id), path("Trin.Clustered"), emit: fasta 

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
cd-hit-est \
  -i ${fasta} \
  -o Trin.Clustered \
  ${params.clustering} \
  -T ${task.cpus} \
  -M 0 
"""
}

process cluster_member {
  container = "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.clustering_out != "") {
    publishDir "results/${params.clustering_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(clstr), path(bak_clstr)

  output:
    path "cluster_*.name", emit: names
    path "mono_cluster.name", emit: mono_cluster 

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
sort -k1n,1 ${bak_clstr} \
  | sed 's/\\.\\.\\.//' \
  | sed 's/>//' \
  | awk 'BEGIN{group = \$1}; {if(\$1 == group){print \$3  >> "cluster_"\$1".name"}; group = \$1}'

# get mono sequence clusters
wc -l cluster_* \
  | grep " 1 " \
  | awk '{system("cat "\$2)}' \
  | sed -E 's/((TRINITY_DN[0-9]+_c[0-9]+_g[0-9]+)_i[0-9]+)/\1 \2/' \
  > mono_cluster_tmp.name

# remove mono sequence clusters
wc -l cluster_* \
  | grep " 1 " \
  | awk '{system("rm "\$2)}'

# reassign mono cluster to cluster with the same gene name
# mono cluster with no match are added to mono_cluster.name
awk '{system("\
  echo \\\$(\
    grep -c \\""\$2"\\" cluster_*.name \
      | sort -t: -k2nr \
      | head -n 1) "\$1" "\$2\
  )}' mono_cluster_tmp.name \
  | sort -t: -k2nr \
  | sed -E 's/name:/name /' \
  | awk '{if(\$2 != 0){print \$3  >> \$1}; if(\$2 == 0){print \$3  >> "mono_cluster.name"}}'  
rm mono_cluster_tmp.name
"""
}

process extract_cluster_member {
  container = "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
  label "small_mem_mono_cpus"
  tag "$file_id"
  if (params.clustering_out != "") {
    publishDir "results/${params.clustering_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)
    tuple val(cluster_id), path(cluster_sequence_name)

  output:
    tuple val(cluster_id), path("*.fasta"), emit: fasta 

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
seqtk subseq ${fasta} ${cluster_sequence_name} \
  | sed -E 's/(>[^ ]+) .*/\\1/' \
  > ${cluster_sequence_name.simpleName}.fasta
"""
}

process align_cluster {
  container = "quay.io/biocontainers/mafft:7.490--h779adbc_0"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.clustering_out != "") {
    publishDir "results/${params.clustering_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)

  output:
    tuple val(file_id), path("*.align"), emit: align

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
mafft --thread ${task.cpus} \
  --quiet \
  --clustalout \
  ${fasta} \
  > ${fasta.simpleName}_cons.align
"""
}

process consensus_cluster {
  container = "quay.io/biocontainers/emboss:6.6.0--h8719169_4"
  label "big_mem_mono_cpus"
  tag "$file_id"
  if (params.clustering_out != "") {
    publishDir "results/${params.clustering_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(align)

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

"""
cons -sequence ${align} \
  -name ${file_id} \
  -outseq ${align.simpleName}_cons.fasta
sed -i 's/>multiple/>${align.simpleName}/' ${align.simpleName}_cons.fasta
"""
}

process merge_cluster_sequences {
  container = "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.clustering_out != "") {
    publishDir "results/${params.clustering_out}", mode: 'copy'
  }

  input:
    tuple val(file_id), path(fasta)
    path cluster_fasta

  output:
    tuple val(file_id), path("concensus_*.fasta"), emit: fasta

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
cat ${cluster_fasta} \
  > concensus_${file_id}.fasta
"""
}