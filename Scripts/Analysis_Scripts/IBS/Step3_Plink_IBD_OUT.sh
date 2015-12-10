#!/bin/bash
#PBS -l mem=1000mb,nodes=1:ppn=1,walltime=72:00:00
#PBS -m abe

DIR=~/50_SNPs

#Create directory where to put only perfect matches files
mkdir $DIR/Analysis

cd $DIR/output
cat $DIR/List_gral.txt | while read line
do

SEGMENT=`echo $line`

gunzip -c ${SEGMENT}out.genome.gz | awk '$12 == "1.0000" {print $0}' >$DIR/Analysis/${SEGMENT}_perfectMatch.txt

done

