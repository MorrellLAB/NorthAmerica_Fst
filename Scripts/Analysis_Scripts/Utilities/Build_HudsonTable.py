#!/usr/bin/env python

"""
Title           :Build_HudsonTable.py
Description     :Takes the tab-delimited text file (samples in rows, markers in columns) without headers,
				 and builds a Hudson-like polymorphism table
Author		 	 :T. Kono

========================================================================================
"""

import sys

#   Dictionary for our data
#   We do it this way because we also need to keep track of the
#   breeding program
g_data = {}
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        tmp = line.strip().split('\t')
        #   For the first line, we want the marker names
        #   which start at field 6 (5 in 0-based counting)
        if index == 0:
            markers = tmp[5:]
        # Othwise, start parsing our data
        else:
            #   The name of our individual (or sample)
            ind = tmp[1]
            g_data[ind] = [x[0] if '0' not in x else 'N' for x in tmp[5:]]

#   First, we are iterating through breeding programs
#   We open a handle to a new file
fname = sys.argv[1].replace('.txt', '_Hudson.txt')
handle = open(fname, 'w')
#   We get he number of individuals in the breeding program
ninds = len(g_data)
#   And the number of markers
nmarkers = len(markers)
#   We write these to our table
handle.write(str(ninds)+'\t'+str(nmarkers)+'\n')
#   Then, we write the names of the markers
handle.write('\t'+'\t'.join(markers)+'\n')
#   Write the line about ancestral state. We don't know it, so we fill with ?
handle.write('anc\t'+'?\t'*nmarkers+'\n')
#   Now we have to print out a file for each breeding program
#   iterate through our dictionary
for i in g_data:
    #   Now we iterate through the individuals within each breeding program
    #   and write the individual name and the marker information
    handle.write(i+'\t'+'\t'.join(g_data[i])+'\n')

#   Close our handle
handle.close()
