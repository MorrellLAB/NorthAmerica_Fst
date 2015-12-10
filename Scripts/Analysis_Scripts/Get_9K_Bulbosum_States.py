#!/usr/bin/env python
"""A script that uses a tab-delimited file and a reference assembly to extract
homologous base calls for pre-defined positions. The tab-delimited file
describes the positions to extract in BED-like format, and the reference
assemblies contain the sequence. Requires Biopython."""

import sys
from Bio import SeqIO

snps_file = sys.argv[1]
morex_seq = sys.argv[2]
bulbosum_seq = sys.argv[3]

barley_snps = {}
#   We also want to save the list of contigs, since this will make it faster to
#   isolate the contigs of interest in the future
snp_contigs = []
#   Read in the SNP data and store it in a dictionary
with open(snps_file, 'r') as f:
    for line in f:
        tmp = line.strip().split()
        #   The final column is the SNP ID
        snp_id = tmp[-1]
        #   The contig is the first column
        contig = tmp[0]
        #   And the position is the second column. It's in 1-based coords, so
        #   we need to subtract 1
        snp_pos = int(tmp[1]) - 1
        #   Then, put it all into the dictionary
        barley_snps[snp_id] = (contig, snp_pos)
        snp_contigs.append(contig)

#   We start a dictionary to hold the ancestral state information
ancestral = {}
#   First, we read through the morex reference assembly and get the morex bases
morex_parser = SeqIO.parse(morex_seq, 'fasta')
for seq in morex_parser:
    #   Check if the seq id is in the list of SNP-containing contigs
    if seq.id in snp_contigs:
        #   Get the right SNP position
        #   But this is a little tricky, since we will iterate through the
        #   dictionary of SNP positions
        for s, position in barley_snps.iteritems():
            if position[0] == seq.id:
                #   Save the morex base
                sys.stderr.write('Finding Morex base of ' + s + '\n')
                #   Sometimes the length isn't right!
                if position[1] >= len(seq.seq):
                    morex_base = 'N'
                else:
                    morex_base = str(seq.seq)[position[1]].upper()
                ancestral[s] = {'morex': morex_base, 'bulbosum': ''}

#   Then we do the same with the bulbosum assembly
bulbosum_parser = SeqIO.parse(bulbosum_seq, 'fasta')
for seq in bulbosum_parser:
    if seq.id in snp_contigs:
        for s, position in barley_snps.iteritems():
            if position[0] == seq.id:
                sys.stderr.write('Finding Bulbosum base of ' + s + '\n')
                if position[1] >= len(seq.seq):
                    bulbosum_base = 'N'
                else:
                    bulbosum_base = str(seq.seq)[position[1]].upper()
                ancestral[s]['bulbosum'] = bulbosum_base

#   Print out a header
print '\t'.join(['SNP_ID', 'MorexBase', 'BulbosumBase'])
#   Next, we just dump it all out
for snp, bases in ancestral.iteritems():
    print '\t'.join([snp, bases['morex'], bases['bulbosum']])
