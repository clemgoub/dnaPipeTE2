version = "2.8.2"
container_url = "lbmc/sratoolkit:${version}"

params.sra_dump = ""
params.sra_dump_out = ""
process sra_dump {
  container = "${container_url}"
  label "big_mem_mono_cpus"
  tag "$sra"
  if (params.sra_dump_out != "") {
    publishDir "results/${params.sra_dump_out}", mode: 'copy'
  }

  input:
    val sra

  output:
    tuple val(sra), path("*.fastq"), emit: fastq

  script:
"""
fastq-dump ${params.sra_dump} --split-files --gzip ${sra}
if [ -f ${sra}_1.fastq.gz ]
then
  mv ${sra}_1.fastq.gz ${sra}_R1.fastq.gz
fi
if [ -f ${sra}_2.fastq.gz ]
then
  mv ${sra}_2.fastq.gz ${sra}_R2.fastq.gz
fi
"""
}
