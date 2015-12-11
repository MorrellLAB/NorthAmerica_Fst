#!/usr/bin/env python

"""
Title           :Calculate_FIS.py  
Description     :Calculates FIS as 1 - (Ho/He)
Author		 	 :T. Kono
Notes		     :Takes one argument: a matrix of genotype data
  		 	 	 Individuals are rows, markers are columns
========================================================================================

"""

import sys

#   a little function to calculate the average
def average(x):
    return sum(x)/float(len(x))

#   Read in the matrix as a whole
matrix = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        matrix.append(tmp)

#   Transpose, so we iterate over markers instead of individuals
t_matrix = zip(*matrix)

fis = []
#   Iterate over markers, and...
for marker in t_matrix:
    #   The number of alleles at the locus
    #   we use this as a float so that we can do fp arithmetic
    t = float(2 * (marker.count('AA') + marker.count('CC') + marker.count('AC') + marker.count('CA')))
    #   there are some cases where there is all missing data? We skip these
    if t == 0.0:
        continue
    #   The frequency of the 'A'
    #   (2*homozygotes + hets) / total alleles
    p = (2 * marker.count('AA') + marker.count('AC') + marker.count('CA')) / t
    #   The frequency of the 'C'
    q = (2 * marker.count('CC') + marker.count('AC') + marker.count('CA')) / t
    #   The expected heterozygosity is 2pq
    He = 2 * p * q
    #   And the observed is just the proportion that are
    Ho = (marker.count('AC') + marker.count('CA')) / t
    #   Sometimes, He is 0, and we have to skip these, since there can't be
    #   "inbreeding" in this case
    if He == 0.0:
        continue
    #   FIS is 1 - (Ho/He)
    f = 1 - (Ho / He)
    #   tack it on to the end
    fis.append(f)

#   Now we just print the average
print average(fis)
