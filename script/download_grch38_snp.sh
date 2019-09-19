#!/bin/sh

ENSEMBL_RELEASE=95
GENOME_FILE=Homo_sapiens.GRCh38.dna.primary_assembly.fa
ENSEMBL_GRCh38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna

DBSNP_RELEASE=151
SNP_FILE=snp${DBSNP_RELEASE}Common.txt
UCSC_COMMON_SNP=http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/${SNP_FILE}

get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget $1
		return $?
	fi
}

if [ ! -f $GENOME_FILE ] ; then
	get ${ENSEMBL_GRCh38_BASE}/$GENOME_FILE.gz || (echo "Error getting $GENOME_FILE" && exit 1)
	gunzip $GENOME_FILE.gz || (echo "Error unzipping $GENOME_FILE" && exit 1)
	mv $GENOME_FILE genome.fa
fi


if [ ! -f $SNP_FILE ] ; then
       get ${UCSC_COMMON_SNP}.gz || (echo "Error getting ${UCSC_COMMON_SNP}" && exit 1)
       gunzip ${SNP_FILE}.gz || (echo "Error unzipping ${SNP_FILE}" && exit 1)
       awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}' ${SNP_FILE} > ${SNP_FILE}.tmp
       mv ${SNP_FILE}.tmp ${SNP_FILE}
fi
