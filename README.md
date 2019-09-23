# ACGTN

ACGTN (Archival Compression for hiGh-Throughput sequeNcing data) is an efficient storage system with the ability to lossless compress >1 billion various reads. Test on kinds of large real datasets, the final storage sizes are reduced to <10% of the original data. Currently, ACGTN works with FASTQ/FASTA files and supports only DNA stream. For quality stream, [QVZ]( https://github.com/mikelhernaez/qvz) can be directly used.


## Installation

ACGTN currently provides Makefile for building on Linux platform. All required libraries are available in directory [deps](https://github.com/compbiolab/acgtn/tree/master/deps). By default, binaries are compiled using _g++_.

First, obtain the repo and its dependencies:

	git clone --recursive https://github.com/compbiolab/acgtn.git

	cd acgtn

Then, install ACGTN's dependencies:

	make deps

When you are ready, build with `make`, and run with `./bin/acgtn`.

You can also produce a static binary with `make static`.


## Usage

ACGTN has three main stages: `index`, `encode` and `decode`. The `index` stage takes input of reference genome and optional SNP database. Once the species-level index has been created, it will be used in the `encode` stage for both DNA-seq and RNA-seq datasets. The `decode` stage will work with the single compressed file and doesn't need the index any more.

The command line of ACGTN is run as:

	acgtn <index|encode|decode> [options] <files|file list>

Run `./bin/acgtn` for available options.

## Examples

Download required genome and SNP files to build index:

	./script/download_grch38_snp.sh
	
	./script/extract_snps_UCSC.py genome.fa snp151Common.txt genome


Build index for ACGTN:

	./bin/acgtn index -o genome -s genome.snp genome.fa

