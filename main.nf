nextflow.enable.dsl=2
nextflow.preview.recursion=true

/* ========================= params =================================*/
params.sra = "SRR14470610";
params.coverage = 0.1;
params.size = "170mb";
params.min_glue = 1;
params.min_contig_length = 200;
params.sample = 3;

println "==========================================================================================="
println "SRA number (--sra):                                " + params.sra;
println "genome size (--size):                              " + params.size;
println "genome coverage (--coverage):                      " + params.coverage;
println "sample number (--sample):                          " + params.sample;
println "trinity min glue (--min_glue):                     " + params.min_glue;
println "trinity min contig length (--min_contig_length):   " + params.min_contig_length;
println "==========================================================================================="

/* ========================= modules import =================================*/
include {
  sra_dump 
} from './src/sra_dump.nf';

include {
  assembly
} from './src/assembly.nf' addParams(
  coverage: params.coverage,
  genome_size: params.size,
  sample: params.sample,
  min_glue: params.min_glue,
  min_contig_length: params.min_contig_length
);

include {
  clustering
} from './src/clustering.nf' addParams(
  min_glue: params.min_glue,
  min_contig_length: params.min_contig_length
);

/* ========================= channel creation =================================*/
channel.from( params.sra.split(" ") ).set{ sra };

workflow {
  sra_dump(sra)
  assembly(sra_dump.out.fastq)
  clustering(assembly.out.fasta)
};