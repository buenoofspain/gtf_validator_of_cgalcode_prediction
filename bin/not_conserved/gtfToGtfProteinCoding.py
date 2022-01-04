#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 09:19:25 2019

@author: niguilla
"""

import argparse

def extract_spe_junction_via_ref(file):
    """
    Function description:
        This function creates a dictionary with reference exon junctions 
        position from a gtf file.
    
    Input:
        - ref_gtf_file: a file who contain all exon information.
    
    Return:
        - junction_genes_dict: a dictionary which contains all junction exons 
                               positions.
    """
    strand = ""
    gene = ""
    count_exon = 0
    nb_tr_tot = 0
    
    junction_list = list()
    
    with open(file, 'r') as ref_file_i:
        lecture = ref_file_i.readlines()
        for line in lecture:
            line = line.replace('\n', '').replace(';', '').replace('"', '').split()
#            print(line)
            junction_list.append((line[0], line[1], line[2], line[3], line[5], line[10], line[11]))
            
    with open("specific_junctions_final.bed", 'w') as output:
        with open(file, 'r') as file_i:
            lecture = file_i.readlines()
            for line in lecture:
                line = line.replace('\n', '').split('\t')
                read_chr_nb = line[0]
                read_start = line[1]
                read_stop = line[2]
                read_id = line[3]
                read_strand = line[5]
                read_block_number = line[9]
                read_block_length = line[10]
                read_block_length_split = read_block_length.split(',')
                read_block_position = line[11]
                read_block_position_split = read_block_position.split(',')
                junction_chr_nb = line[12]
                junction_start = line[13]
                junction_stop = line[14]
                junction_id = line[15]
                junction_strand = line[17]
                junction_block_number = line[21]
                junction_block_length = line[22]
                junction_block_length_split = junction_block_length.split(',')
                junction_block_position = line[23]
                
                
                if read_block_number == "2" and read_chr_nb == junction_chr_nb:
#                    if junction_id == "CGACAFTENST00000488267JUNC00020":
#                        print(line)
                    start_first_alignment_read = int(read_start)+1
                    end_first_alignment_read = int(read_start) + int(read_block_length_split[0])
                    
                    start_second_alignment_read = start_first_alignment_read + int(read_block_position_split[1])
                    end_second_alignment_read = int(read_stop)
                    
                    start_junction = int(junction_start) + int(junction_block_length_split[0])
                    end_junction = int(junction_stop) - int(junction_block_length_split[1]) + 1
                    
#                        print(end_first_alignment_read, start_junction, start_second_alignment_read, end_junction)
                    if end_first_alignment_read == start_junction and start_second_alignment_read == end_junction: 
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                        
                    elif end_first_alignment_read == start_junction and start_second_alignment_read+1 == end_junction: 
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                
    nb_junctions = len(junction_genes_list)
    
    print(">>>EXTRACT SPECIFIC EXON COUPLES")
#    print("Number of genes:", len(junction_genes_list))
    print("Number of transcripts:", nb_tr_tot)
    print("Number of exon junctions:", nb_junctions)
    
#    print(junction_genes_list)
    return junction_genes_list


def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="A SCRIPT TO OBTAIN SPECIFIC EXON JUNCTION NOT KNOWN IN ENSEMBL.")
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-f", "--file",
                        dest = "gtf_file", 
                        metavar = "GTF FILE", 
                        help = "Give a gtf file to obtain a gtf file with protein_coding tag.",
                        required = True)
    
    return parser.parse_args()       
    
if __name__ == "__main__":
    
    OPTIONS = parse()
    
    gtf_file_i = OPTIONS.gtf_file
    
    extract_spe_junction_via_ref(gtf_file_i)