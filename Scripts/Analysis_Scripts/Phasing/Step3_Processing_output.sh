#!/bin/bash


DIR=~/output/output_processing

for CHR in 1 2 3 4 5 6 7
do
sed '1,21d' MYoutput_${CHR}_hapguess_switch.out |sed '$d'|sed '/#/d'  >./temp.txt

#Samples_NorthAm_prephase_Duplicates.txt is a list of samples names doubled, for example Sample1 and Sample1_2
paste $DIR/Samples_NorthAm_prephase_Duplicates.txt $DIR/temp.txt  >~/output/output_processing/MyOutput_genotypes/MYoutput_${CHR}.txt

rm ./temp.txt
done

