#!/usr/bin/env python
"""Calculates derived allele frequecy, using the H. bulbosum states of the BOPA
SNPs listed in Fang et al. 2014 in G3 (Table S3). Sites that are trans-specific
between barley and H. bulbosum are ignored. Heterozygous calls are treated as
missing data."""

import sys

#   Take two input files - Table S3 from Fang et al. 2014 and the genotypes
bulbosum_states = sys.argv[1]
genotyping_matrix = sys.argv[2]

#   Read in and store the data. We will store the ancestral data as a
#   dictionary, and the genotyping data as a list of lists
ancestral = {}
with open(bulbosum_states, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            snp_id = tmp[0]
            morex = tmp[1]
            bulbosum = tmp[2]
            ancestral[snp_id] = bulbosum

genotypes = []
with open(genotyping_matrix, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            #   In the genotpying matrix, each SNP name has an X in front
            snpnames = line.strip().split()
            #   Remove that X
            snpnames = [s[1:] for s in snpnames]
        else:
            #   Put the genotyping data into the matrix. We pull off the first
            #   column, since it is just individual name, and we do not need
            #   to save that.
            genotypes.append(line.strip().split()[1:])
#   Next we zip the genotype matrix, since we want to iterate over markers and
#   not individuals.
genotypes = zip(*genotypes)

#   Now, we can start to calculate up the derived allele frequencies
for snpid, calls in zip(snpnames, genotypes):
    #   If, for some reason, the SNP does not have an inferred ancestral sate
    #   we skip it.
    if snpid not in ancestral:
        continue
    #   Let's get the bulbosum state
    anc = ancestral[snpid]
    #   First, check if bulbosum is segregating at that SNP
    if len(anc) > 1:
        continue
    #   There are also a couple isolated cases where bulbosum has bases that
    #   are completely different than barley. For instance, if bulbosum has
    #   a T, but barley is segregating for A/G
    if anc not in ''.join(calls):
        continue
    #   Then, we go through and remove missing calls and heterozygous calls
    hom = [a if (len(set(a)) == 1 and a != 'NN') else None for a in calls]
    #   Next, count up the number of times we see the bulbosum state.
    #   Since we have already checked the conditions that would violate our
    #   assumptions, we can just do a straight count
    #   Start counters for the number of non-missing individuals, and the
    #   number of times we see the ancestral state
    n_inds = 0
    n_anc = 0
    #   Iterate through the list of genotype calls and spit out the derived
    #   allele frequency.
    for snp in hom:
        #   Check if the SNP isn't NoneType
        if snp:
            n_inds += 1
            if anc in snp:
                n_anc += 1
        else:
            continue
    #   Then, calculate the *ancestral* allele frequency
    anc_freq = float(n_anc)/n_inds
    #   And the derived allele frequency is 1-ancestral
    derived_freq = 1 - anc_freq
    print derived_freq
