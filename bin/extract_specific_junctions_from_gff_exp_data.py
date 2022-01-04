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

def read_files(pred_gtf_file, ref_gtf_file, delta_extension, junction_list=False):
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
    
    #list_specific_view = [(26177716, 26137222, '-'), (24738323, 24843812, '+'), (10046519, 10051081, '+'), (21007221, 21003337, '-'), (21003329, 21000088, '-'), (62387155, 62391338, '+'), (62391442, 62394597, '+'), (59514109, 59513429, '-'), (59513395, 59512910, '-'), (54890871, 54899934, '+'), (25885932, 25900372, '+'), (25819513, 25852271, '+'), (18291817, 18398022, '+'), (5427455, 5493860, '+'), (45278942, 45329424, '+'), (45334916, 45339015, '+'), (45329534, 45334125, '+'), (21118198, 21117695, '-'), (21117623, 21117114, '-'), (55551588, 55547030, '-'), (55547016, 55544121, '-'), (20053916, 20064007, '+'), (20064201, 20071367, '+'), (20057215, 20071367, '+'), (20064201, 20071373, '+'), (20057215, 20069472, '+'), (42586586, 42586331, '-'), (42586263, 42585902, '-'), (31320555, 31325566, '+'), (76336242, 76328952, '-'), (11558980, 11555828, '-'), (24817165, 24821247, '+'), (24767772, 24821247, '+'), (54589076, 54587986, '-'), (55406304, 55479124, '+'), (55406304, 55504653, '+'), (55451203, 55504653, '+'), (12863021, 12858832, '-'), (72012524, 72005031, '-'), (42966739, 42983409, '+'), (301204, 302101, '+'), (303076, 304849, '+'), (62331450, 62319733, '-'), (28104968, 28129025, '+'), (28134156, 28138842, '+'), (110090437, 110089147, '-'), (34983480, 34977872, '-'), (35000163, 34982292, '-'), (35000163, 34977872, '-'), (39935844, 39927803, '-'), (23226995, 23227731, '+'), (102078695, 102078486, '-'), (25844835, 25854240, '+'), (25847909, 25854240, '+'), (25844835, 25849240, '+'), (38461117, 38461775, '+'), (33856270, 33862981, '+'), (12093973, 12078121, '-'), (12078113, 12077320, '-'), (20590937, 20595290, '+'), (62528763, 62546438, '+'), (21977720, 21982442, '+'), (8881709, 8881570, '-'), (11464269, 11462477, '-'), (54890154, 54894977, '+'), (54892663, 54894977, '+'), (44765160, 44776069, '+'), (52447076, 52475540, '+'), (42375534, 42374243, '-'), (43792467, 43766355, '-'), (18023399, 18023739, '+'), (43053829, 43055150, '+'), (43055490, 43064085, '+'), (22742645, 22723368, '-'), (49247003, 49245099, '-'), (42401433, 42402450, '+'), (106853767, 106851780, '-'), (106853857, 106851777, '-'), (28529606, 28517187, '-')]
    #list_specific_view2 = [(26177716, 26137222, '-'), (24738323, 24843812, '+'), (10046519, 10051081, '+'), (21007221, 21003337, '-'), (21003329, 21000088, '-'), (62391442, 62394597, '+'), (59514109, 59513429, '-'), (59513395, 59512910, '-'), (54890871, 54899934, '+'), (25885932, 25900372, '+'), (25819513, 25852271, '+'), (18291817, 18398022, '+'), (5427455, 5493860, '+'), (45334916, 45339015, '+'), (45329534, 45334125, '+'), (21118198, 21117695, '-'), (55551588, 55547030, '-'), (55547016, 55544121, '-'), (20053916, 20064007, '+'), (20064201, 20071367, '+'), (20057215, 20071367, '+'), (20064201, 20071373, '+'), (20057215, 20069472, '+'), (42586586, 42586331, '-'), (42586263, 42585902, '-'), (31320555, 31325566, '+'), (76336242, 76328952, '-'), (11558980, 11555828, '-'), (24817165, 24821247, '+'), (24767772, 24821247, '+'), (54589076, 54587986, '-'), (55406304, 55479124, '+'), (55406304, 55504653, '+'), (55451203, 55504653, '+'), (12863021, 12858832, '-'), (72012524, 72005031, '-'), (42966739, 42983409, '+'), (303076, 304849, '+'), (62331450, 62319733, '-'), (28104968, 28129025, '+'), (28134156, 28138842, '+'), (35000163, 34982292, '-'), (23226995, 23227731, '+'), (102078695, 102078486, '-'), (25847909, 25854240, '+'), (38461117, 38461775, '+'), (33856270, 33862981, '+'), (12093973, 12078121, '-'), (12078113, 12077320, '-'), (20590937, 20595290, '+'), (62528763, 62546438, '+'), (21977720, 21982442, '+'), (8881709, 8881570, '-'), (11464269, 11462477, '-'), (54890154, 54894977, '+'), (54892663, 54894977, '+'), (44765160, 44776069, '+'), (52447076, 52475540, '+'), (42375534, 42374243, '-'), (43792467, 43766355, '-'), (18023399, 18023739, '+'), (43053829, 43055150, '+'), (43055490, 43064085, '+'), (22742645, 22723368, '-'), (42401433, 42402450, '+'), (106853767, 106851780, '-'), (106853857, 106851777, '-')]

    
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
    
    if junction_list:
        output = "specific_junction_cga_not_in_exp.bed"
    else:
        output = "specific_junction_cga.bed"
    
    with open(pred_gtf_file, 'r') as pred_file_i:
        with open(output, 'w') as junction_bed_file_o:
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
                            
                        try:
                            if not junction in junction_ref_list and not junction in junction_pred_list and junction in junction_list:# and junction in list_specific_view:# and junction in list_specific_view2:
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
                                        
                                    junction_bed_file_o.write(str(chr_nb) +'\t'+ str(exon_junction_prev) +'\t'+ str(exon_junction_next) +'\t'+ transcript+"JUNC"+name_zero+str(junc_nb) +'\t.\t'+ strand +'\t'+ str(exon_junction_prev) +'\t'+ str(exon_junction_next) +'\t.\t'+ str(2) +'\t'+ str(extension1)+','+str(extension2) +'\t'+ str(0)+','+str(beg1-end2+extension1-1) +'\n')
    #                                            junction_bed_file_o.write(transcript+'\n')
                                    count_specific_junction_traited += 1
                                    
                                    len_exon_all.append(end2-beg2)
                                    len_exon_all.append(end1-beg1)
                                    len_exon_part1.append(end1-beg1)
                                    len_exon_part2.append(end2-beg2)
    
                                    
                                
                            beg1 = beg2
                            end1 = end2
                        
                        except TypeError:
                            if not junction in junction_ref_list and not junction in junction_pred_list:# and junction in list_specific_view:# and junction in list_specific_view2:
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

    return junction_pred_list

def extract_specific_junction_not_in_experimental(pred_junction_list, exp_gtf_file, delta_extension, pred_gtf_file, ref_gtf_file):
    """
    Function description:
        
    """
#    junction_exp_list = extract_reference_exons_chaining(exp_gtf_file)
#    print(pred_junction_list)
    junction_exp_list = list()
    
    strand = ""
    gene = ""
    count_junction = 0
    count_specific_junction_traited = 0
    count_exon = 0
    
    len_exon_all = []
    len_exon_part1 = []
    len_exon_part2 = []
    
    
    delta = delta_extension
    junc_nb = 0
    
    with open(exp_gtf_file, 'r') as exp_file_i:
        with open("junction_exp.txt", 'w') as junction_exp_o:
            lecture = exp_file_i.readlines()
            for line in lecture:
                line_read = line.replace('\n', '').replace('"', '').replace(';', '').split()
                if line.startswith("#"):
                    pass
#                if line_read[2] == "gene":
#                    gene = line_read[8].split('Name=Gene:')[1]
                        
                elif line_read[2] == "transcript":
#                    gene = line_read[8].split('Parent=Gene:')[1]
#                    transcript = line_read[11]
#                    chr_nb = line_read[0]
                    strand = line_read[6]
#                    exon_nb = line_read[13]
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
#                        print(beg1, end1, strand, "1er")
                        
                    
                    elif count_exon > 1:
                        beg2 = int(line_read[3])
                        end2 = int(line_read[4])
#                        print(beg2, end2, strand)    
#                        if strand == "+":
                        junction = (end1, beg2, strand)
#                            print(junction)
                            
#                        elif strand == "-":
#                            junction = (beg2, end1, strand)
#                            print(junction)
#
#                            print(junction[2])
                        for j in pred_junction_list:
                            if junction[2] == j[2] and junction[0] == j[0] and junction[1] == j[1]:
#                                print(junction)
#                                if strand == "+":
#                                    print(end1, beg2, strand, j)
#                                else:
#                                    print(beg2, end1, strand, j)
                                pred_junction_list.remove(j)
                            elif junction[2] == j[2] and junction[0] == j[1] and junction[1] == j[0]:
#                                print(junction)
#                                if strand == "+":
#                                    print("1", end1, beg2, strand, j)
#                                else:
#                                    print("2", beg2, end1, strand, j)
                                pred_junction_list.remove(j)
                                
                        junction_exp_o.write(str(junction)+'\n')
                        count_junction += 1
                        
                        beg1 = beg2
                        end1 = end2           
                        
    print(len(pred_junction_list))
    #print(pred_junction_list)
    print("Nb junction in Feelnc", count_junction)
    
    read_files(pred_gtf_file, ref_gtf_file, delta_extension, pred_junction_list)
    
    return #junction_pred_list

        
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
    
    parser.add_argument("-e", "--exp",
                        dest = "gtf_exp_file",
                        metavar = "EXON GTF FILE FROM EXPERIMENTAL",
                        help = "Give an exon gtf file with exon from an experimental procedure to search specific exon from prediction.",
                        required = True)
    
    parser.add_argument("-d", "--delta",
                        dest = "extension",
                        metavar = "DELTA CHOSE FOR EXON EXTENSION",
                        help = "Give a delta score for the exon extension.",
                        default = 65,
                        required = False,
                        type = int)
    
    
    return parser.parse_args()       
    
if __name__ == "__main__":
    
    OPTIONS = parse()
    
    pred_file = OPTIONS.gtf_pred_file
    ref_file = OPTIONS.gtf_ref_file
    exp_file = OPTIONS.gtf_exp_file
    delta = OPTIONS.extension
    
    junction_pred_list = read_files(pred_file, ref_file, delta)
    
    print("-----------------------------------------------")
    print("EXTRACT SPECIFIC JUNCTION NOT IN FEELNC RESULTS")
    print("-----------------------------------------------")
    extract_specific_junction_not_in_experimental(junction_pred_list, exp_file, delta, pred_file, ref_file)
    
