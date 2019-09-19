#!/usr/bin/python

import re
from collections import deque
from collections import defaultdict
from argparse import ArgumentParser, FileType


def read_genome(genome_fname):
    chr_dict = defaultdict(list)
    prev_chr, sequence, size = "", "", 0
    for line in open(genome_fname):
        if line.startswith(">"):
            if prev_chr:
                chr_dict[prev_chr].append(size)
            chr_name = line.strip().split()[0][1:]
            chr_dict[chr_name].append(size)
            prev_chr = chr_name
        else:
            line = line.strip().upper()
            size += len(line)
            sequence += line
    chr_dict[prev_chr].append(size)
    return chr_dict, sequence


def generate_snps(snp_file, snp_que):
    used_snp = [snp for snp in snp_que if snp[-1]]
    assert len(used_snp) < 2
    if used_snp:
        return snp_que

    snp = snp_que[-1]
    line = "%s\t%s\t%d\t%s\n" % (snp[3], snp[1], snp[2], snp[4])
    snp_file.write(line)
    snp_que[-1][-1] = True
    return snp_que


def main(genome_fname, snp_fname, base_fname, dist):
    rc_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    chr_dict, sequence = read_genome(genome_fname)
    snp_out_file = open(base_fname + ".snp", 'wt')
    snp_que = deque(maxlen=dist)

    pos_seen = set()
    for line in open(snp_fname):
        if not line or line.startswith('#'):
            continue

        fields = line.rstrip('\r\n').split('\t')
        """
        id, chr, start, end, rs_id, score, strand, refNCBI, refUCSC, observed, molType, classType, valid, \
            avHet, avHetSE, func, locType, weight, exceptions, submitterCount, submitters, \
            alleleFreqCount, alleles, alleleNs, alleleFreqs, bitfields = fields
        """
        _, chrom, start, end, rsid, _, strand, _, _, _, molType, classType = fields[:12]
        if molType != "genomic" or classType != "single" or chrom not in chr_dict:
            continue

        start, end = int(start), int(end)
        assert start < chr_dict[chrom][1] - chr_dict[chrom][0]
        pos = chr_dict[chrom][0] + start
        if start + 1 != end or pos in pos_seen:
            continue
        pos_seen.add(pos)

        freqs, alleles, supports = [], [], []
        if re.match("^([01]\.\d+?,){2,}$", fields[-1]) is None:
            freqs = fields[-2].rstrip(',').split(',')
            alleles = fields[-4].rstrip(',').split(',')
            supports = fields[-3].rstrip(',').split(',')
        else:
            freqs = fields[-1].rstrip(',').split(',')
            alleles = fields[-3].rstrip(',').split(',')
            supports = fields[-2].rstrip(',').split(',')

        assert len(alleles) == len(freqs)
        assert len(alleles) == len(supports)
        freqs = [float(f) for f in freqs]
        supports = [float(f) for f in supports]

        if strand == "-":
            alleles = [rc_dict[x] if x in "ACGT" else x for x in alleles]

        ref_base = sequence[pos]
        if ref_base == 'N' or ref_base not in alleles:
            continue

        snp_info = [
            [alleles[x], freqs[x]]
            for x in range(len(alleles))
            if alleles[x] in "ACGT" and alleles[x] != ref_base
            and freqs[x] >= 0.01 and supports[x] >= 100.0
        ]

        if snp_info:
            snp_info = sorted(snp_info, key=lambda x: x[1], reverse=True)
            curr_snp = [pos, chrom, start, rsid, snp_info[0][0], False]

            while snp_que:
                if snp_que[0][0] + dist <=  pos:
                    snp_que.popleft()
                else:
                    break

            snp_que.append(curr_snp)
            snp_que = generate_snps(snp_out_file, snp_que)

    snp_out_file.close()


if __name__ == '__main__':
    parser = ArgumentParser(
        description=
        'Extract SNPs and haplotypes from a SNP file downloaded from UCSC (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz)'
    )
    parser.add_argument(
        'genome_fname',
        nargs='?',
        type=str,
        help='input genome file (e.g. genome.fa)')
    parser.add_argument(
        'snp_fname',
        nargs='?',
        type=str,
        help=
        'input snp file downloaded from UCSC (plain text or gzipped file is accepted: snp144Common.txt or snp144Common.txt.gz)'
    )
    parser.add_argument(
        "base_fname",
        nargs='?',
        type=str,
        help="base filename for SNPs and haplotypes")

    args = parser.parse_args()
    if not args.genome_fname or not args.snp_fname or not args.base_fname:
        parser.print_help()
        exit(1)
    main(args.genome_fname, args.snp_fname, args.base_fname, 54)
