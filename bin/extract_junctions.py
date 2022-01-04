#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:40:00 2019

@author: Nicolas Guillaudeux
@team: Dyliss
"""

import argparse

def extract_reference_exons_chaining(ref_gtf_file, delta_extension):
    """
    Function description:
        This function creates a dictionary with reference exon junctions 
        position from a gtf file and extract all new exon predicted to compare 
        with a reference information.
    
    Input:
        - ref_gtf_file: a file which contains all exon information from a
                        reference.
        - delta_extension: length of exon-exon junction.
        
    Output:
        - junction_reference.bed: a bed file which contains all exon junction
                                  positions.
    
    Return:
        - junction_genes_list: a list which contains all exon junction 
                               positions.
    """
    strand = ""
    gene = ""
    
    count_junction = 0
    count_junction_traited = 0
    
    count_exon = 0
    nb_tr_tot = 0
    
    junction_genes_list = list()
    len_exon_all = []
    len_exon_part1 = []
    len_exon_part2 = []
    
    delta = delta_extension
    junc_nb = 0
    
    
    with open(ref_gtf_file, 'r') as ref_file_i:
        with open("junction_reference.bed", "w") as junction_bed_file_o:
            lecture = ref_file_i.readlines()
            for line in lecture:
                line = line.replace('\n', '').replace(';', '').replace('"', '').split()
                
#                print(line)
                if line[2] == "gene":
                    gene = line[9]
                
                if line[2] == "transcript":
                    nb_tr_tot += 1
                    gene = line[9]
                    transcript = line[11]
                    chr_nb = line[0]
                    strand = line[6]
                    exon_number = line[13]
                    beg1 = 0
                    end1 = 0
                    beg2 = 0
                    end2 = 0
                    count_exon = 0
                    
                if line[2] == "exon":
                    count_exon += 1
                    
                    if count_exon == 1:
                        beg1 = int(line[3])
                        end1 = int(line[4])
                    if count_exon > 1:
                        beg2 = int(line[3])
                        end2 = int(line[4])
                        
                        if strand == "+":
                            junction = (end1, beg2, strand)
                            
                        elif strand == "-":
                            junction = (beg1, end2, strand)
                        
                        
                        if not junction in junction_genes_list:
                            junction_genes_list.append(junction)
                            count_junction += 1
                            junc_nb +=1
                            if len(str(junc_nb)) == 1:
                                name_zero = "0000"
                            elif len(str(junc_nb)) == 2:
                                name_zero = "000"
                            elif len(str(junc_nb)) == 3:
                                name_zero = "00"
                            elif len(str(junc_nb)) == 4:
                                name_zero = "0"
                                
                            if strand == "+":
                                strand_nb = "1"
                                extension1 = 0
                                extension2 = 0
                                
                                if end1-beg1+1 >= delta:
                                    exon_junction_prev = end1-delta
                                    extension1 = delta
                                elif end1-beg1+1 < delta:
                                    exon_junction_prev = beg1-1
                                    extension1 = end1-beg1+1
                                
                                if end2-beg2+1 >= delta:
                                    exon_junction_next = beg2+delta
                                    extension2 = delta
                                else:
                                    exon_junction_next = end2
                                    extension2 = end2-beg2+1
                                    
                                junction_bed_file_o.write(str(chr_nb) +'\t'+ str(exon_junction_prev) +'\t'+ str(exon_junction_next) +'\t'+ transcript+"JUNC"+name_zero+str(junc_nb) +'\t.\t'+ strand +'\t'+ str(exon_junction_prev) +'\t'+ str(exon_junction_next) +'\t.\t'+ str(2) +'\t'+ str(extension1)+','+str(extension2) +'\t'+ str(0)+','+str(beg2-end1+extension1-1) +'\n')
#                                            junction_bed_file_o.write(transcript+'\n')
                                count_junction_traited += 1
                                
                                len_exon_all.append(end1-beg1)
                                len_exon_all.append(end2-beg2)
                                len_exon_part1.append(end1-beg1)
                                len_exon_part2.append(end2-beg2)
                                    
                            elif strand == "-":
                                strand_nb = "-1"
                                extension1 = 0
                                extension2 = 0
                                end_prev = end2
                                beg_next = beg1
                                
                                if end2-beg2+1 >= delta:
                                    exon_junction_prev = end2-delta
                                    extension1 = delta
                                elif end2-beg2+1 < delta:
                                    exon_junction_prev = beg2-1
                                    extension1 = end2-beg2+1
                                    
                                if end1-beg1+1 >= delta:
                                    exon_junction_next = beg1+delta
                                    extension2 = delta
                                else:#if end1-beg1+1 < delta:
                                    exon_junction_next = end1
                                    extension2 = end1-beg1+1
                                    
                                junction_bed_file_o.write(str(chr_nb) +'\t'+ str(exon_junction_prev) +'\t'+ str(exon_junction_next) +'\t'+ transcript+"JUNC"+name_zero+str(junc_nb) +'\t.\t'+ strand +'\t'+ str(exon_junction_prev) +'\t'+ str(exon_junction_next) +'\t.\t'+ str(2) +'\t'+ str(extension1)+','+str(extension2) +'\t'+ str(0)+','+str(beg1-end2+extension1-1) +'\n')
#                                            junction_bed_file_o.write(transcript+'\n')
                                count_junction_traited += 1
                                
                                len_exon_all.append(end2-beg2)
                                len_exon_all.append(end1-beg1)
                                len_exon_part1.append(end1-beg1)
                                len_exon_part2.append(end2-beg2)
                                
                            
                        beg1 = beg2
                        end1 = end2
                
    nb_junctions = len(junction_genes_list)
    
    print(">>>EXTRACT SPECIFIC EXON COUPLES")
#    print("Number of genes:", len(junction_genes_list))
    print("Number of transcripts:", nb_tr_tot)
    print("Number of exon junctions:", nb_junctions)
    print("Number of specific exon junctions traited:", count_junction_traited, "/", count_junction)
    
#    print(junction_genes_list)
    return junction_genes_list

        
def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="A SCRIPT TO OBTAIN SPECIFIC EXON JUNCTION NOT KNOWN IN ENSEMBL.")
    # Create command-line parser for all options and arguments to give
    
    parser.add_argument("-r", "--ref",
                        dest = "gtf_ref_file",
                        metavar = "EXON GTF FILE FROM REFERENCE",
                        help = "Give an exon gtf file with exon from reference to search specific exon from prediction.",
                        required = True)
    
    parser.add_argument("-d", "--delta",
                        dest = "extension",
                        metavar = "DELTA CHOSE FOR EXON EXTENSION",
                        help = "Give a delta score for the exon extension.",
                        default = 25,
                        required = False,
                        type = int)
    
    
    return parser.parse_args()       
    
if __name__ == "__main__":
    
    OPTIONS = parse()
    
    ref_file = OPTIONS.gtf_ref_file
    delta = OPTIONS.extension
    
    extract_reference_exons_chaining(ref_file, delta)
