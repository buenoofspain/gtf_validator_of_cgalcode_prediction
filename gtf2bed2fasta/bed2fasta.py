#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:40:00 2019

@author: Nicolas Guillaudeux
@team: Dyliss
"""

import argparse

def extract_list_fasta_from_tr_list(fasta_file, tr_list):
    """
    Function description:
        
    """
    list_tr = list()
    
    with open(tr_list, 'r') as f_i:
        lecture = f_i.readlines()
        for tr in lecture:
            tr = tr.replace('\n', '')
            if not tr in list_tr:
                list_tr.append(tr)
            
    print("Number of transcripts:", len(list_tr))
    
    list_tr_uniq = list()
    with open(fasta_file, 'r') as f_i:
        output_file_name = tr_list.split(".")[0]+".fasta"
        with open(output_file_name, 'w') as f_o:
            lecture = f_i.readlines()
            for line in lecture:
                line = line.replace("\n", "")
                if line.startswith(">"):
                    line_test = 0
                    tr_id = line.split(" ")[1]
                    if tr_id in list_tr and not tr_id in list_tr_uniq:
                        f_o.write(line+"\n")
                        list_tr_uniq.append(tr_id)
                        line_test = 1
                else:
                    if line_test == 1:
                        f_o.write(line+"\n")
    
def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="A SCRIPT TO OBTAIN SPECIFIC EXON JUNCTION NOT KNOWN IN ENSEMBL.")
    # Create command-line parser for all options and arguments to give
    
    parser.add_argument("-f", "--fasta",
                        dest = "fasta_ref_file",
                        metavar = "FASTA FILE WITH ALL PREDICTED SEQUENCES",
                        help = "Give a fasta file with all coding sequences from transcripts.",
                        required = True)
    
    parser.add_argument("-t", "--tr",
                        dest = "tr_set_file",
                        metavar = "TRANSCRIPTS LIST",
                        help = "Give a file with a transcript list.",
                        required = True)
    
    
    return parser.parse_args()       
    
if __name__ == "__main__":
    
    OPTIONS = parse()
    
    fasta_file = OPTIONS.fasta_ref_file
    list_tr = OPTIONS.tr_set_file
    
#    fasta_file = "dog.fasta"
#    list_tr = "dog_tr_622.txt"
    
    extract_list_fasta_from_tr_list(fasta_file, list_tr)
    