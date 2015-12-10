#!/usr/bin/env python
"""Calculates derived allele frequecy, using the H. bulbosum states of the BOPA
SNPs listed in Fang et al. 2014 in G3 (Table S3). Sites that are trans-specific
between barley and H. bulbosum are ignored."""

import sys

#   Define the missing value here
missing = 'NN'

#   Take two input files - Table S3 from Fang et al. 2014 and the genotypes
bulbosum_states = sys.argv[1]
genotyping_matrix = sys.argv[2]


def derived_allele_frequency(genotypes, derived):
    """Uses a supplied derived state to calculate the derived allele
    frequency in a list of genotype calls. Expects the calls to represent
    a vector of calls for an individual marker across samples, and the calls
    to be diploid (2 characters long). Returns a single floating point number
    for the derived allele frequency."""
    #   Count up the number of derived alleles in each call, and sum them
    n_derived = 0
    n_notmissing = 0
    for call in genotypes:
        #   check if the call is missing data, skip if it is missing
        if call == missing:
            continue
        else:
            #   Count up the number of derived alleles
            derived_count = call.count(derived)
            n_derived += derived_count
            #   Add two to the notmissing value, since we are dealing with
            #   diploid calls. For haploid calls, we would use 1.
            n_notmissing += 2
    #   Then calculate the frequency of the derived allele
    derived_freq = float(n_derived)/float(n_notmissing)
    return(derived_freq)


def ancestral_derived(anc):
    """Read in the file that contains ancestral and derived allele information
    and return two dictionaries that contain it. Remove markers with missing
    information and trans-specific polymorphisms. Returns two dictionaries:
    the ancestral and derived states, respectively."""
    #   Read in and store the data. We will store the ancestral data as a
    #   dictionary, and the genotyping data as a list of lists
    ancestral = {}
    derived = {}
    with open(anc, 'r') as f:
        for index, line in enumerate(f):
            #   Skip the header line
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                snp_id = tmp[0]
                morex = tmp[1]
                bulbosum = tmp[2]
                #   If H. bulbosum is segregating, missing, or identical to
                #   morex, then we skip it
                if bulbosum == 'N' or len(bulbosum) != 1 or morex == bulbosum:
                    continue
                else:
                    ancestral[snp_id] = bulbosum
                    derived[snp_id] = morex
    return(ancestral, derived)


def read_genotypes(g_mat):
    """Read in the genotyping matrix and return it as a transposed matrix.
    Assumes that rows are individuals and columns are markers. Returns a list
    of SNP names and a list of tuples for the genotyping data."""
    genotypes = []
    with open(g_mat, 'r') as f:
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
    return(snpnames, genotypes)


def main():
    """Main function. Dispatches to other functions to calculate the derived
    allele frequency for each SNP that it can. Prints a newline-delimited list
    of floating point numbers."""
    #   Parse the SNP names and genotyping data out of the genotype matrix
    snp_ids, genotypes = read_genotypes(genotyping_matrix)
    #   And the ancestral and derived states
    anc, deriv = ancestral_derived(bulbosum_states)
    #   Then, we iterate through the genotypes, and get the derived allele
    #   frequencies
    for snp, calls in zip(snp_ids, genotypes):
        #   What is the derived allele for the SNP?
        if snp not in deriv:
            continue
        else:
            derived_allele = deriv[snp]
            derived_freq = derived_allele_frequency(calls, derived_allele)
            print derived_freq


#   Now, we run the whole thing
main()
