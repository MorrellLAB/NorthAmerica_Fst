
Datasets information
=======================
Original Genotyping files can be downloaded from www.triticeaetoolbox.org.

Alias_CAPname_Xrefnames.txt: Alias for each sample in the data set. Alias represent the program to which each accession belong, and a numeric
count of samples within a breeding program

samples_information.txt: Has the basic information for each sample, easy to use for analysis

Sample_information_FULL.csv: All samples' information downloaded from T3 website. I modified it by adding my new sample assigination (Alias) and changing the 29 BushAg Internation accessions' program assignation from "BA" to "BAI"


Barley_NorthAm_QC_AB_no_duplicates_or.txt  (and Barley_NorthAm_QC_ACTG_no_duplicates_or.txt): 
Genotype file (for genotypes and nucleotides) after QC controling for >25% missingness, duplicated 
individuals and missing information (row type, growth habit) and removal of single samples representing a population.
from this population. SNPs are sorted by their genetic map position according to Munoz et al 2011. Total individuals 3,613 x 2,542 SNPs


Barley_Annotations.txt : SNP annotations using SNP contextual sequences compared to a barley sequence to infer their position. We used SNPmeta (Kono et al., 2014) to obtain this information.
Customized R script (SortingSNPs_genMapOrder.R) was used to order the SNPs into a genetic map position. This file will be use in Fst and PHS analysis.

Fang_et_al_TableS3.txt: From Fang et al 2014, got H.Bulbosum states to be used as ancestral states.

GeneticMap_T3_020315.txt: BOPA SNPs genetic map (Munoz et. al., 2011) downloaded from T3 for 2,994 SNPs.

List_NIL_duplicates.txt: Accessions (using alias names) that are Near Isogenic Lines. 

Lines_missing_sampleInfo.txt: Accessions missing row-type or growth habit information that were removed during QC, as indicated in the manuscript.
