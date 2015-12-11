#!/bin/bash
#Title           :Step3_Plink_IBS_OUT.sh
#Description     :Select genomic segments with IBD 
#Author		 	 :A. Poets
#Date			 :Jun 26, 2015
#Note		     :Requires output files from Step2_Plink_IBS_input.R
#========================================================================================

DIR=~/50_SNPs

#Create directory where to put only IBS segments with maximum 10% mistmatch
mkdir $DIR/Analysis

cd $DIR/output
cat $DIR/List_gral.txt | while read line
do

SEGMENT='echo $line'

gunzip -c ${SEGMENT}out.genome.gz | awk '$12 == "1.0000" {print $0}' >$DIR/Analysis/${SEGMENT}_perfectMatch.txt

done

