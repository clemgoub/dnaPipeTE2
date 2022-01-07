nextflow.enable.dsl=2
/* ========================= modules import =================================*/
include {
  sra_dump 
} from './src/sra_dump.nf'

/* ========================= params =================================*/
params.sra = "SRR14470610"
println params.sra;


/* ========================= channel creation =================================*/
channel.from( params.sra.split(" ") ).set{ sra }

workflow {
  sra_dump(sra)
}