nextflow.enable.dsl=2
/* ========================= params =================================*/
params.sra = "SRR14470610";
params.coverage = 0.1;
params.size = "170mb";
params.min_glue = 1;
params.min_contig_length = 200;
params.seed = 1234; //Math.round(Math.random() * 1e9);
params.sample = "3";

println "==========================================================================================="
println "SRA number:                  " + params.sra;
println "genome size:                 " + params.size;
println "genome coverage:             " + params.coverage;
println "sample number:               " + params.sample;
println "trinity min glue:            " + params.min_glue;
println "trinity min contig length:   " + params.min_contig_length;
println "==========================================================================================="

/* ========================= modules import =================================*/
include {
  sra_dump 
} from './src/sra_dump.nf';

include {
  sampling 
} from './src/sampling.nf' addParams(
  coverage: params.coverage,
  genome_size: params.size,
  seed: params.seed
);

include {
  assembly 
} from './src/assembly.nf' addParams(
  min_glue: params.min_glue,
  min_contig_length: params.min_contig_length
);

/* ========================= channel creation =================================*/
def str_to_int_list(str) {
  return (1..Integer.parseInt(str)).toList();
}

channel.from( params.sra.split(" ") ).set{ sra };
channel.fromList( str_to_int_list(params.sample) ).set{ samples };

workflow {
  sra_dump(sra)
  sampling(
    samples
      .combine(sra_dump.out.fastq)
      .map{ it -> it[1..2] }
  )
  assembly(sampling.out.fastq)
};