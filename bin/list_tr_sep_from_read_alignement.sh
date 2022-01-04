#! /bin/bash

#Developped by Nicolas Guillaudeux
#Dyliss team (INRIA Rennes Bretagne Atlantique/IRISA/Université de Rennes 1)
#02/2020

# Description

echo "Number of junctions:"
awk -F\JUNC '{print $1}' nonull_list.txt | wc -l

echo "Number of transcripts with junctions with read alignment:"
awk -F\JUNC '{print $1}' nonull_list.txt | sort -u | wc -l

echo "Extract this list in nonull_list_tr.txt"
awk -F\JUNC '{print $1}' nonull_list.txt | sort -u > nonull_list_tr.txt

echo "Number of junctions:"
awk -F\JUNC '{print $1}' null_list.txt | wc -l

echo "Number of transcripts with junctions with no read alignment:"
awk -F\JUNC '{print $1}' null_list.txt | sort -u | wc -l

echo "Extract this list in null_list_tr.txt"
awk -F\JUNC '{print $1}' null_list.txt | sort -u > null_list_tr.txt
