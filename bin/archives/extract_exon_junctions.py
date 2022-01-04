#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:40:00 2019

@author: Nicolas Guillaudeux
@team: Dyliss
"""

import argparse
import requests, sys
import matplotlib.pyplot as plt

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
    
    junction_genes_dict = dict()
    
    with open(ref_gtf_file, 'r') as ref_file_i:
        lecture = ref_file_i.readlines()
        for line in lecture:
            line = line.replace('\n', '').replace(';', '').replace('"', '').split()
            
            if line[2] == "gene":
                gene = line[9]
                if gene not in junction_genes_dict:
                    junction_genes_dict[gene] = []
            
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
                if line[6] == "+":
                    count_exon += 1
                    if count_exon == 1:
                        beg1 = int(line[3])
                        end1 = int(line[4])
                    if count_exon > 1:
                        beg2 = int(line[3])
                        end2 = int(line[4])
                        
                        junction = (end1, beg2, strand)
                        
                        try:
                            if not junction in junction_genes_dict[gene]:
                                junction_genes_dict[gene].append(junction)
                        except KeyError:
                            junction_genes_dict[gene] = [junction]
                        
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
                        
                        try:
                            if not junction in junction_genes_dict[gene]:
                                junction_genes_dict[gene].append(junction)
                        except KeyError:
                            junction_genes_dict[gene] = [junction]
                        
                        beg2 = beg1
                        end2 = end1
                
    nb_junctions = 0
    for list_of_junctions in junction_genes_dict.values():
        nb_junctions += len(list_of_junctions)
    
    print(">>>EXTRACT SPECIFIC EXON COUPLES")
    print("Number of genes:", len(junction_genes_dict))
    print("Number of transcripts:", nb_tr_tot)
    print("Number of exon junctions:", nb_junctions)
    
    
    return junction_genes_dict
    
def tr_validated(tr_file):
    """
    Function description:
        This function transform a list of transcript within a file to a list.
        
    Input:
        - tr_file: the file with the list of transcripts.
        
    Return:
        - tr_list: a list of transcripts.
    """
    tr_list = list()
    
    with open(tr_file, 'r') as file_i:
        lecture = file_i.readlines()
        for line in lecture:
            tr = line.replace('\n', '')
            print(line)
            if tr not in tr_list:
                tr_list.append(tr)
            
    return tr_list

def obtain_sequence_from_ensembl(chrNb, seqBeg, seqEnd, strand):
    """
    Function description:
        From http://rest.ensembl.org/documentation/info/sequence_region
        This function request Ensembl from genomics positions and return the
        string sequence.
        
    Input:
        - chrNb: The chromosome reference number (1, 2..., X)
        - seqBeg: The first position of the genomics position (< seqEnd).
        - seqEnd: The last position of the genomics position (> seqBeg).
        - strand: The strand of the reading frame (1 or -1).
        
    Return:
        - The sequence in string format.
        
    """
    server = "http://rest.ensembl.org"
    ext = "/sequence/region/dog/"+chrNb+":"+seqBeg+".."+seqEnd+":"+strand+"?"
    #To return the output in fasta format:
#    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
    #To return the sequence without description:
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
 
    if not r.ok:
        r.raise_for_status()
        sys.exit()
 
 
    return r.text

def read_files(pred_gtf_file, ref_gtf_file, tr_validated, delta_extension):
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
    junction_ref_dict = extract_reference_exons_chaining(ref_gtf_file)
    
    tr_list = list()
    
    with open(tr_validated, 'r') as file_i:
        lecture = file_i.readlines()
        for line in lecture:
            tr = line.replace('\n', '')
            if tr not in tr_list:
                tr_list.append(tr)
    
    junction_pred_dict = dict()
    
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
        with open("specific_exon_junctions.fasta", 'w') as junction_file_o:
            with open("specific_junction.bed", 'w') as junction_bed_file_o:
                lecture = pred_file_i.readlines()
                for line in lecture:
                    line_read = line.replace('\n', '').replace('"', '').replace(';', '').split()
                    
                    if line_read[2] == "gene":
                        gene = line_read[9]
                        if gene not in junction_pred_dict:
                            junction_pred_dict[gene] = []
                            
                    if line_read[2] == "transcript":
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
                        
                    if line_read[2] == "exon":
                        count_exon += 1
                        
                        if count_exon == 1:
                            beg1 = int(line_read[3])
                            end1 = int(line_read[4])
                            
                        if count_exon > 1:# and count_exon != nb_exon_tot:
                            beg2 = int(line_read[3])
                            end2 = int(line_read[4])
                            
                            if strand == "+":
                                junction = (end1, beg2, strand)
                                
                            elif strand == "-":
                                junction = (beg1, end2, strand)
                                
                            
#                            if transcript == "CGACAFTENST00000380138":
#                                print(junction in junction_ref_dict[gene])
#                                print((end1, beg2, strand) in junction_ref_dict[gene])
#                                print(end1, beg2, strand)
#                                print(junction_ref_dict[gene])
#                            if transcript == "CGACAFTENST00000318445":
#                            if transcript == "CGACAFTENST00000318445":
#                                print(beg1, end1, beg2, end2)
#                                if not junction in junction_ref_dict[gene]: 
#                                    print("No")
#                                else:
#                                    print("Yes")
#                                print(junction_ref_dict[gene])
                                
                            
                            if not junction in junction_ref_dict[gene]:
                                try:
#                                    if not junction in junction_ref_dict[gene]:
                                    junction_ref_dict[gene].append(junction)
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
#                                        if end2-beg2+1 < 10 or end1-beg1+1 < 10:
#                                            print(transcript, end2-beg2+1, end1-beg1+1)
                                        if transcript in tr_list:
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
# =============================================================================
#                                            part_of_first_exon = obtain_sequence_from_ensembl(str(chr_nb), str(exon_junction_prev), str(end1), str(strand_nb))
#                                            part_of_second_exon = obtain_sequence_from_ensembl(str(chr_nb), str(beg2), str(exon_junction_next), str(strand_nb))
#                                             
#                                            junction_file_o.write("> "+transcript+" "+gene+" "+str(chr_nb)+" "+str(exon_junction_prev)+"-"+str(end1)+" "+str(beg2)+"-"+str(exon_junction_next)+" "+strand_nb+'\n')
#                                            junction_file_o.write(part_of_first_exon+'|'+part_of_second_exon+'\n')
# =============================================================================
                                            
                                            junction_bed_file_o.write(str(chr_nb) +'\t'+ str(exon_junction_prev) +'\t'+ str(exon_junction_next) +'\t'+ transcript+"JUNC"+name_zero+str(junc_nb) +'\t.\t'+ strand +'\t'+ str(exon_junction_prev) +'\t'+ str(exon_junction_next) +'\t.\t'+ str(2) +'\t'+ str(extension1)+','+str(extension2) +'\t'+ str(0)+','+str(beg2-end1+extension1-1) +'\n')
#                                            junction_bed_file_o.write(transcript+'\n')
                                            count_specific_junction_traited += 1
                                            
                                            len_exon_all.append(end1-beg1)
                                            len_exon_all.append(end2-beg2)
                                            len_exon_part1.append(end1-beg1)
                                            len_exon_part2.append(end2-beg2)
                                            
                                    if strand == "-":
#                                        if end2-beg2+1 < 10 or end1-beg1+1 < 10:
#                                            print(transcript, end2-beg2+1, end1-beg1+1)
                                        if transcript in tr_list:
                                            strand_nb = "-1"
                                            extension1 = 0
                                            extension2 = 0
                                            end_prev = end2
                                            beg_next = beg1
    #                                        print(obtain_sequence_from_ensembl(str(chr_nb), str(end_prev), str(beg_next), str(strand_nb)))
                                            
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
                                                
# =============================================================================
#                                            part_of_first_exon = obtain_sequence_from_ensembl(str(chr_nb), str(exon_junction_prev), str(end2), str(strand_nb))
#                                            part_of_second_exon = obtain_sequence_from_ensembl(str(chr_nb), str(beg1), str(exon_junction_next), str(strand_nb))
#                                        
#                                            junction_file_o.write("> "+transcript+" "+gene+" "+str(chr_nb)+" "+str(exon_junction_prev)+"-"+str(end2)+" "+str(beg1)+"-"+str(exon_junction_next)+" "+strand_nb+'\n')
#                                            junction_file_o.write(part_of_first_exon+'|'+part_of_second_exon+'\n')
# =============================================================================

                                            junction_bed_file_o.write(str(chr_nb) +'\t'+ str(exon_junction_prev) +'\t'+ str(exon_junction_next) +'\t'+ transcript+"JUNC"+name_zero+str(junc_nb) +'\t.\t'+ strand +'\t'+ str(exon_junction_prev) +'\t'+ str(exon_junction_next) +'\t.\t'+ str(2) +'\t'+ str(extension1)+','+str(extension2) +'\t'+ str(0)+','+str(beg1-end2+extension1-1) +'\n')
#                                            junction_bed_file_o.write(transcript+'\n')
                                            count_specific_junction_traited += 1
                                            
                                            len_exon_all.append(end2-beg2)
                                            len_exon_all.append(end1-beg1)
                                            len_exon_part1.append(end1-beg1)
                                            len_exon_part2.append(end2-beg2)

                                except KeyError:
                                    junction_ref_dict[gene] = [junction]
                                    
                                
                            beg1 = beg2
                            end1 = end2
                    
    print("Number of specific exon junctions traited:", count_specific_junction_traited, "/", count_specific_junction)
    
#    plt.hist(len_exon_all)
##    print(len_exon_all)
#    plt.title('Histogram of length of part left in adjacence couple:', fontsize=10)
##    plt.savefig("len_exon_all.png")
##    plt.show()
#    
#    plt.hist(len_exon_part1)
##    print(len_exon_part1)
#    plt.title('Histogram of length of part left in adjacence couple:', fontsize=10)
##    plt.savefig("len_exon_part1.png")
##    plt.show()
#    
#    plt.hist(len_exon_part2)
##    print(len_exon_part2)
#    plt.title('Histogram of length of part right in adjacence couple:', fontsize=10)
##    plt.savefig("len_exon_part2.png")
##    plt.show()
    
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
    
    parser.add_argument("-t", "--tr",
                        dest = "tr_validated_list",
                        metavar = "TRANSCRIPTS ALREADY VALIDATED",
                        help = "Give a transcripts validated list to extract specific exons of this transcripts.",
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
    tr_validated = OPTIONS.tr_validated_list
    delta = OPTIONS.extension
    
    read_files(pred_file, ref_file, tr_validated, delta)
