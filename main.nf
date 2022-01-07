nextflow.enable.dsl=2
/* ========================= params =================================*/
params.sra = "SRR14470610"
params.coverage = ""
params.size = ""
println params.sra;

/* ========================= modules import =================================*/
include {
  sra_dump 
} from './src/sra_dump.nf'

include {
  sampling 
} from './src/sampling.nf' addParams(
  coverage: params.coverage,
  genome_size: params.size
)

include {
  assembly 
} from './src/assembly.nf' addParams(
  min_glue: params.min_glue,
  min_contig_length: params.min_contig_length
)

/* ========================= channel creation =================================*/
channel.from( params.sra.split(" ") ).set{ sra }

workflow {
  sra_dump(sra)
  sampling(sra_dump.out.fastq)
}