#!/usr/bin/perl -wpl

#s/(\w+_2)//ig;
s/^(AB_\d+\_2|BA_\d+\_2|MN_\d+\_2|N2_\d+\_2|N6_\d+\_2|OR_\d+\_2|UT_\d+\_2|VT_\d+\_2|WA_\d+\_2|MT_\d+\_2)//ig;

s/^(AB_\d+|BA_\d+|MN_\d+|N2_\d+|N6_\d+|OR_\d+|UT_\d+|VT_\d+|WA_\d+|MT_\d+)/\#$1\n/ig;
s/\t//g;

