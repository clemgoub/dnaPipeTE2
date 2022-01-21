# dnaPipeTE2

The version 2 of the [dnaPipeTE pipeline](https://github.com/clemgoub/dnaPipeTE) aims to improve the performance of the original pipeline and to simplify its installation.

## Installation

Make sure that you have java 8 or later is installed
on your computer by using the command:

```bash
java -version 
```

Then, you need to install `nextflow`.
The following commands creates a file `nextflow` in the current dir.

```bash
export NXF_EDGE=1
curl -s https://get.nextflow.io | bash
./nextflow self-update
```

Then you can run the pipeline with the following command:

```bash
./nextflow run clemgoub/dnaPipeTE2
```

## Usage

To run the pipeline on your data you can happend the following command line options :

- `--sra` : an SRA accession number of your fastq file (default: `"SRR144706"`)
- `--fastq` : the path of your fastq files (example: `"SRR144706.fastq.gz"` for single-end data or `"SRR144706_R{1,2}.fastq.gz"` for paired-end data)
- `--size` : the estimated genome size, in `kb`, `mb`, `gb` or `tb` (default: `"170mb"` to match the `"SRR144706"` example)
- `--coverage` : the genome coverage that you want to sample, it must be < 1.0 (default: `0.1`)
- `--sample` : the number of sample to build to enrich the assembly in transposable element reads (default: `3`)
- `--dfam_db` : use dfam database for the annotation, instead of RepeatMasker default on (default: `false`)
- `--custom_db` : use custom database for the annotation, instead of RepeatMasker default on (default: `""`)
- `--repeatmakser_threshold` : threshold to filter RepeatMasker results (default: `0.0`)

For example to run the pipeline on the following fastq files: `SRR144706_R1.fastq.gz` and `SRR144706_R2.fastq.gz`

```bash
./nextflow run clemgoub/dnaPipeTE2 --fastq "SRR144706_R{1,2}.fastq.gz" --size "170mb" --coverage 0.1
```

For the same run using dfam database:

```bash
./nextflow run clemgoub/dnaPipeTE2 --fastq "SRR144706_R{1,2}.fastq.gz" --size "170mb" --coverage 0.1 --dfam_db true
```

# [Containerization](https://en.wikipedia.org/wiki/OS-level_virtualization)

By default `dnaPipeTE2` is run using *Docker* container, but you can use *singularity* container instead with the `-profile singularity` option.




