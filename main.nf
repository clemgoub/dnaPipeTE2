nextflow.enable.dsl=2
/* ========================= params =================================*/
params.sra = "SRR14470610";
params.coverage = 0.1;
params.size = "170mb";
params.min_glue = 1;
params.min_contig_length = 200;
params.sample = "3";

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
  sampling 
} from './src/sampling.nf' addParams(
  coverage: params.coverage,
  genome_size: params.size
);

include {
  assembly;
  merge;
  final_assembly
} from './src/assembly.nf' addParams(
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
def str_to_int_list(str) {
  return (1..Integer.parseInt(str)).toList();
}

def SRA_to_sample_name(str, sample) {
  return(str + "_sample_" + sample);
}

channel.from( params.sra.split(" ") ).set{ sra };
channel.fromList( str_to_int_list(params.sample) ).set{ samples };

workflow {
  sra_dump(sra)
  sampling(
    samples
      .combine(sra_dump.out.fastq)
      .map{ it -> [it[1], SRA_to_sample_name(it[1], it[0]), it[2]] }
  )
  assembly(sampling.out.fastq)
  merge(assembly.out.folder.groupTuple())
  final_assembly(sra_dump.out.fastq, merge.out.folder)
  clustering(final_assembly.out.fasta)
};