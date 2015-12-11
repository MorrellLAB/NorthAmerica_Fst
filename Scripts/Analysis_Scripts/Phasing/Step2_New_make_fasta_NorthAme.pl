#!/usr/bin/perl -wpl

#Title           :Step2_New_make_fasta_NorthAme.pl 
#Description     :Format files for fastPhase (Scheet and Stephens 2006) step2. Represent 
#				  each individual by two chromosomes.
#Author		 	 :A. Poets 
#Note		     :Requires output from Step1_PHASE_input_NorthAm.R
#Usage			 : Step2_New_make_fasta_NorthAme.pl  NorthAm_prephase_chr_1.txt >./INPUT_NorthAm_chr1.txt
#========================================================================================


s/^(AB_\d+\_2|BA_\d+\_2|MN_\d+\_2|N2_\d+\_2|N6_\d+\_2|OR_\d+\_2|UT_\d+\_2|VT_\d+\_2|WA_\d+\_2|MT_\d+\_2)//ig;

s/^(AB_\d+|BA_\d+|MN_\d+|N2_\d+|N6_\d+|OR_\d+|UT_\d+|VT_\d+|WA_\d+|MT_\d+)/\#$1\n/ig;
s/\t//g;
s/\ /\t/g;
