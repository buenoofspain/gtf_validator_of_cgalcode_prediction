#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:40:00 2019

@author: Nicolas Guillaudeux
@team: Dyliss
"""

import argparse

def extract_reference_exons_chaining(ref_gtf_file):
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
    
    junction_genes_list = list()
    
    with open(ref_gtf_file, 'r') as ref_file_i:
        lecture = ref_file_i.readlines()
        for line in lecture:
            line = line.replace('\n', '').replace(';', '').replace('"', '').split()
            
#            if line[2] == "gene":
#                gene = line[9]
            
            if line[2] == "transcript":
                nb_tr_tot += 1
                gene = line[9]
#                transcript = line[11]
                strand = line[6]
                beg1 = 0
                end1 = 0
                beg2 = 0
                end2 = 0
                count_exon = 0
                
            if line[2] == "exon":
#                if "protein_coding" in line:
                if line[6] == "+":
                    count_exon += 1
                    if count_exon == 1:
                        beg1 = int(line[3])
                        end1 = int(line[4])
                    if count_exon > 1:
                        beg2 = int(line[3])
                        end2 = int(line[4])
                        
                        junction = (end1, beg2, strand)
                        
                        if not junction in junction_genes_list:
                            junction_genes_list.append(junction)
                        
                        beg1 = beg2
                        end1 = end2
                        
                elif line[6] == "-":
                    count_exon += 1
                    if count_exon == 1:
                        beg2 = int(line[3])
                        end2 = int(line[4])
                    if count_exon > 1:
                        beg1 = int(line[3])
                        end1 = int(line[4])
                        
                        junction = (beg2, end1, strand)
                        
                        if not junction in junction_genes_list:
                            junction_genes_list.append(junction)
                        
                        beg2 = beg1
                        end2 = end1
                
    nb_junctions = len(junction_genes_list)
    
    print(">>>EXTRACT SPECIFIC EXON COUPLES")
#    print("Number of genes:", len(junction_genes_list))
    print("Number of transcripts:", nb_tr_tot)
    print("Number of exon junctions:", nb_junctions)
    
#    print(junction_genes_list)
    return junction_genes_list

def read_files(pred_gtf_file, ref_gtf_file, delta_extension):
    """
    Function description:
        This function extract all new exon predicted to compare with a
        reference information.
    
    Input:
        - pred_gtf_file; a file which contains all exon information from 
                         prediction.
        - ref_gtf_file: a file which contains all exon information from a
                        reference.
                        
    Output:
        - output.exonsSpecificPrediction.gtf: a file which contains all new
                                              predicted exons.
    """
    # Extraire les junctions de chaque transcrit dans un dictionnaire
    junction_ref_list = extract_reference_exons_chaining(ref_gtf_file)
    
    junction_pred_list = list()
    
    strand = ""
    gene = ""
    count_specific_junction = 0
    count_specific_junction_traited = 0
    count_exon = 0
    
    len_exon_all = []
    len_exon_part1 = []
    len_exon_part2 = []
    
    
    delta = delta_extension
    junc_nb = 0
    
    with open(pred_gtf_file, 'r') as pred_file_i:
        with open("specific_junction_cga.bed", 'w') as junction_bed_file_o:
            lecture = pred_file_i.readlines()
            for line in lecture:
                line_read = line.replace('\n', '').replace('"', '').replace(';', '').split()
                
                if line_read[2] == "gene":
                    gene = line_read[9]
                        
                elif line_read[2] == "transcript":
                    gene = line_read[9]
                    transcript = line_read[11]
                    chr_nb = line_read[0]
                    strand = line_read[6]
                    exon_nb = line_read[13]
                    beg1 = 0
                    end1 = 0
                    beg2 = 0
                    end2 = 0
                    count_exon = 0
                    
                elif line_read[2] == "exon":
                    count_exon += 1
                    
                    if count_exon == 1:
                        beg1 = int(line_read[3])
                        end1 = int(line_read[4])
                        
                    elif count_exon > 1:
                        beg2 = int(line_read[3])
                        end2 = int(line_read[4])
                        
                        if strand == "+":
                            junction = (end1, beg2, strand)
                            
                        elif strand == "-":
                            junction = (beg1, end2, strand)
                            
                        
                        if not junction in junction_ref_list and not junction in junction_pred_list:
                            junction_pred_list.append(junction)
                            count_specific_junction += 1
                            junc_nb += 1
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
                                count_specific_junction_traited += 1
                                
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
                                count_specific_junction_traited += 1
                                
                                len_exon_all.append(end2-beg2)
                                len_exon_all.append(end1-beg1)
                                len_exon_part1.append(end1-beg1)
                                len_exon_part2.append(end2-beg2)

                                
                            
                        beg1 = beg2
                        end1 = end2
                    
    print("Number of specific exon junctions traited:", count_specific_junction_traited, "/", count_specific_junction)
    
    return
        
def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="A SCRIPT TO OBTAIN SPECIFIC EXON JUNCTION NOT KNOWN IN ENSEMBL.")
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-p", "--pred",
                        dest = "gtf_pred_file", 
                        metavar = "EXON GTF FILE FROM PREDICTION", 
                        help = "Give an exon gtf file with exon from prediction to search specific exon from prediction.",
                        required = True)
    
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
    
    pred_file = OPTIONS.gtf_pred_file
    ref_file = OPTIONS.gtf_ref_file
    delta = OPTIONS.extension
    
    read_files(pred_file, ref_file, delta)
