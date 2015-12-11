#!/bin/bash

#Title           :Step3_Processing_output.pl 
#Description     :Add sample names to inputed chromosomes.
#Author		 	 :A. Poets 
#Note		     :Requires output from fastPHASE
#========================================================================================

DIR=~/output/output_processing

for CHR in 1 2 3 4 5 6 7
do
sed '1,21d' MYoutput_${CHR}_hapguess_switch.out |sed '$d'|sed '/#/d'  >./temp.txt

#Samples_NorthAm_prephase_Duplicates.txt is a list of samples. Each sample is present 
#twice therefore Sample1 and Sample1_2 are two chromosomes of same individual.
paste $DIR/Samples_NorthAm_prephase_Duplicates.txt $DIR/temp.txt  >~/output/output_processing/MyOutput_genotypes/MYoutput_${CHR}.txt

rm ./temp.txt
done

