nextflow.enable.dsl=2
nextflow.preview.recursion=true

/* ========================= params =================================*/
params.sra = "SRR14470610";
params.fastq = "";
params.coverage = 0.1;
params.size = "170mb";
params.min_glue = 1;
params.min_contig_length = 200;
params.sample = 3;
params.dfam_db = false
params.custom_db = ""
params.repeatmasker_threshold = 0.0

println "==========================================================================================="
if (params.fastq == "") {
  println "SRA number (--sra):                                    " + params.sra;
} else {
  println "fastq file (--fastq):                                  " + params.fastq;
}
println "genome size (--size):                                  " + params.size;
println "genome coverage (--coverage):                          " + params.coverage;
println "sample number (--sample):                              " + params.sample;
println "trinity min glue (--min_glue):                         " + params.min_glue;
println "trinity min contig length (--min_contig_length):       " + params.min_contig_length;
println "use dfam database (--dfam_db):                         " + params.dfam_db;
println "use custom database (--custom_db):                     " + (params.custom_db == "" ? "no" : params.custom_db)
println "use repeatmakser threshold (--repeatmakser_threshold): " + params.repeatmasker_threshold;
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

include {
  annotation
} from './src/annotation.nf' addParams(
  dfam_db: params.dfam_db,
  custom_db: params.custom_db,
  repeatmasker_threshold: params.repeatmasker_threshold,
);

include {
  quantification
} from './src/quantification.nf'

include {
  analysis 
} from './src/analysis.nf'

/* ========================= channel creation =================================*/
if (params.fastq == "") {
  channel.from( params.sra.split(" ") ).set{ sra };
  channel.empty().set{ fastq };
} else {
  channel.empty().set{ sra };
  channel.fromFilePairs( params.fastq ).set{ fastq };
}

if (params.dfam_db) {
  channel.fromPath( "https://www.dfam.org/releases/current/families/Dfam.h5.gz" ).map{it -> [it.simpleName, it]}.set{ dfam_db };
} else {
  channel.empty().set{ dfam_db };
}
if (params.custom_db == "") {
  channel.empty().set{ custom_db };
} else {
  channel.fromPath( params.custom_db ).set{ custom_db };
}

workflow {
  sra_dump(sra)
  assembly(sra_dump.out.fastq.mix(fastq))
  clustering(assembly.out.fasta)
  annotation(
    clustering.out.fasta,
    dfam_db,
    custom_db
  )
  quantification(assembly.out.super_transcript, fastq)
  analysis(
    quantification.out.tsv,
    annotation.out.annotation,
    clustering.out.cluster,
    assembly.out.super_transcript,
    assembly.out.super_transcript_gtf
  )
}