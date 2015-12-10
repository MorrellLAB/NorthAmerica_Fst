#!/usr/bin/env python
import sys
Usage = """
Change the name of genotype files for CAP names.
Read starndard input. Two files are required:
1)The first one is the genotypes with SNPs as header, and samples as rows.
2)The second file contain a list of names, 1st column the new names (short),
2nd has the long names from T3

%./Step2_Change_names.py Barley_genotypes_QC_AB.txt Alias_CAPname_Xrefnames.txt 
"""


if len(sys.argv) < 3:
	print Usage
else:

	InputFile = sys.argv[1]
	NamesFile = sys.argv[2]
	# InputFile=('genotype_test.dat')
	INPUT = open(InputFile,'r')
	# NamesFile=('name_test.txt')
	NAMES = open(NamesFile, 'r')

	OutputName = ("Barley_NorthAm_QC_AB.txt")
	OutFile = open(OutputName, 'w')

# Create a dictionary
	FileDic = {}
	for line in NAMES:
		# Remove if empty key
		linenew = line.strip('\n').split('\t')
		NAME_ID = linenew[0]
		NAME_long = linenew[1]
		if NAME_long not in FileDic: FileDic[NAME_long] = []
		FileDic[NAME_long] = NAME_ID

	for index, line in enumerate(INPUT):
		if index == 0:
			OutFile.write(line)
		else:
			# Modify names in genotype file
			linenew = line.strip('\n').split('\t')
			linenew_names = linenew[0]
			if linenew_names in FileDic:
				New_name = FileDic[linenew_names]
			else:
				New_name = linenew_names
			OUTPUT = line.replace(linenew_names, New_name)
			
			OutFile.write(OUTPUT)

INPUT.close()
NAMES.close()
OutFile.close()
