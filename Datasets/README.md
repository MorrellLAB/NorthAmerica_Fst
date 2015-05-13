
Datasets information
=======================

ANA_names_Barley_NorthAm_QC_AB_no_duplicates_or.txt  (and ANA_names_Barley_NorthAm_QC_ACTG_no_duplicates_or.txt): 
Genotype file (for genotyeps and nucleotides) after QC controling for >25% missingness, duplicated 
individuals and missing information (row type, growth habit) and removal of one MT sample that was the only Spring 6-row 
from this population.SNPs are sorted by their genetic map position according to Munoz et al 2011. Total individuals 3,637 x 2,542 SNPs

ANAnames_samples_information.txt: Has the basic information for each sample

Sample_information_addNAMES.csv: All samples' information downloaded from T3 website. I modified it by adding my new sample assigination and changing the 29 BushAg Internation accessions' program assignation from "BA" to "BAI"

Customized R script (SortingSNPs_genMapOrder.R) was used to order the SNPs into a genetic map position. This file will be use in Fst and PHS analysis.

Fang_et_al_TableS3.txt: From Fang et al 2014, got H.Bulbosum states to be used as ancestral states.

GeneticMap_T3_020315.txt: BOPA SNPs genetic map (Munoz et. al., 2011) downloaded from T3 for 2,994 SNPs.
 
