#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 15:43:01 2019

@author: Nicolas Guillaudeux
@team: Dyliss
"""
import argparse

def add_correct_reads(dict_reads_in_junction, junction_id, read_id):
    '''
    '''
    try:
        if read_id not in dict_reads_in_junction[junction_id]:
            dict_reads_in_junction[junction_id].append(read_id)
    except KeyError:
        dict_reads_in_junction[junction_id] = [read_id]
        
    return dict_reads_in_junction

def read_file(alignment_file):#, junction_file):
    '''
    FUNCTION DESCRIPTION:
        This function use an alignment file and check if reads are correctly 
        mapped on specific junctions desired.
    
    INPUT:
        alignment_file: a file obtain with BedTools intersect in BED12 format 
                        and with -wa, -wb, -split and -bed options.
                        Example: bedtools intersect \
                                -abam file.bam \
                                -b junction.bed \
                                -wa \
                                -wb \
                                -split \ 
                                -bed \
                                > output_BedToolsIntersect_specific_junction.bed
                                
    OUTPUT:
        output_reads_mapped_on_specific_junctions.bed: a bed file with all
                                                       desired reads
                                                       
    RETURN:
        junction_count_read: a dictionary with all reads align correctly for
                             each junction.
    '''
    
    junction_count_read = dict()
    
    nb_read_ok = 0
    with open("output_reads_mapped_on_specific_junctions.bed", 'w') as output:
        with open(alignment_file, 'r') as alignment_i:
            lecture_alignment = alignment_i.readlines()
            for line_alignment in lecture_alignment:
                line = line_alignment.replace('\n', '').split('\t')
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
                        
                        
#                    if junction_strand == "+" or junction_strand == "-":# and "JUNC00101" in junction_id:
#                        if read_strand == "+" or read_strand == "-":
#                            
#                            if int(read_start)+int(read_block_length_split[0]) == int(junction_start)+int(junction_block_length_split[0]):
#                                #TO DO : si +1 ?
#                                if int(read_stop)-int(read_block_length_split[1])+1 == int(junction_stop)-int(junction_block_length_split[1]):
#                                    try:
#                                        junction_count_read[junction_id].append(read_id)
#                                    except KeyError:
#                                        junction_count_read[junction_id] = [read_id]
#                                    nb_read_ok += 1
#                                    output.write(line_alignment)
#                                elif int(read_stop)-int(read_block_length_split[1]) == int(junction_stop)-int(junction_block_length_split[1]):
##                                    print(junction_chr_nb, read_id, read_start, read_stop, read_block_length_split, junction_id, junction_start, junction_stop, junction_block_length_split)
##                                    print(int(read_start)+int(read_block_length_split[0]), int(junction_start)+int(junction_block_length_split[0]))
##                                    print(int(read_stop)-int(read_block_length_split[1]), int(junction_stop)-int(junction_block_length_split[1]))
#                                    try:
#                                        junction_count_read[junction_id].append(read_id)
#                                    except KeyError:
#                                        junction_count_read[junction_id] = [read_id]
#                                    nb_read_ok += 1
#                                    output.write(line_alignment)
#
##                    elif junction_strand == "-":
##                        if read_strand == "+" or read_strand == "-":
##                             if int(read_start)+int(read_block_length_split[0]) == int(junction_start)+int(junction_block_length_split[0]):
##                                 if int(read_stop)-int(read_block_length_split[1]) == int(junction_stop)-int(junction_block_length_split[1]):
##                                     try:
##                                         junction_count_read[junction_id].append(read_id)
##                                     except KeyError:
##                                         junction_count_read[junction_id] = [read_id]
##                                     nb_read_ok += 1
##                                     output.write(line_alignment)
##                                 elif int(read_stop)-int(read_block_length_split[1])+1 == int(junction_stop)-int(junction_block_length_split[1]):
##                                    print(junction_chr_nb, read_id, read_start, read_stop, read_block_length_split, junction_id, junction_start, junction_stop, junction_block_length_split)
##                                    print(int(read_start)+int(read_block_length_split[0]), int(junction_start)+int(junction_block_length_split[0]))
##                                    print(int(read_stop)-int(read_block_length_split[1])+1, int(junction_stop)-int(junction_block_length_split[1]))
##                                    try:
##                                        junction_count_read[junction_id].append(read_id)
##                                    except KeyError:
##                                        junction_count_read[junction_id] = [read_id]
##                                    nb_read_ok += 1
##                                    output.write(line_alignment)
                elif read_block_number == "3":
#                    if junction_id == "CGACAFTENSMUST00000096766JUNC00135":# and ("CGACAFTENSMUST00000096766JUNC00135" in junction_id or "CGACAFTENSMUST00000096766JUNC00136" in junction_id):# and read_chr_nb == junction_chr_nb:
    #                    print('yes', read_id, junction_id, read_block_number)
    #                    print(junction_chr_nb, read_id, read_start, read_stop, read_block_number, read_block_length_split, junction_id, junction_start, junction_stop, junction_block_length_split)
    #                    print(line_alignment)
                    start_first_alignment_read = int(read_start)+1
                    end_first_alignment_read = int(read_start) + int(read_block_length_split[0])
                    
                    start_second_alignment_read = start_first_alignment_read + int(read_block_position_split[1])
                    end_second_alignment_read = int(read_start) + int(read_block_position_split[1]) + int(read_block_length_split[1])
                    
                    start_third_alignment_read = start_first_alignment_read + int(read_block_position_split[2])
                    end_third_alignment_read = int(read_stop)
#                    print(str(start_first_alignment_read)+"-"+str(end_first_alignment_read), str(start_second_alignment_read)+"-"+str(end_second_alignment_read), str(start_third_alignment_read)+"-"+str(end_third_alignment_read))
                    
                    start_junction = int(junction_start) + int(junction_block_length_split[0])
                    end_junction = int(junction_stop) - int(junction_block_length_split[1]) + 1
                    
                    # If the junction is in the first alignment (1 and 2)
                    if end_first_alignment_read == start_junction and start_second_alignment_read == end_junction:
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                        
                    elif end_first_alignment_read == start_junction and start_second_alignment_read+1 == end_junction:
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                        
                    
                    # If the junction is in the second alignment (2 and 3)
                    elif end_second_alignment_read == start_junction and start_third_alignment_read == end_junction:
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                        
                    elif end_second_alignment_read == start_junction and start_third_alignment_read+1 == end_junction:
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)


                elif read_block_number == "4":
                    start_first_alignment_read = int(read_start)+1
                    end_first_alignment_read = int(read_start) + int(read_block_length_split[0])
                    
                    start_second_alignment_read = start_first_alignment_read + int(read_block_position_split[1])
                    end_second_alignment_read = int(read_start) + int(read_block_position_split[1]) + int(read_block_length_split[1])
                    
                    start_third_alignment_read = start_first_alignment_read + int(read_block_position_split[2])
                    end_third_alignment_read = int(read_start) + int(read_block_position_split[2]) + int(read_block_length_split[2])
                    
                    start_fourth_alignment_read = start_first_alignment_read + int(read_block_position_split[3])
                    end_fourth_alignment_read = int(read_stop)
#                    print(str(start_first_alignment_read)+"-"+str(end_first_alignment_read), str(start_second_alignment_read)+"-"+str(end_second_alignment_read), str(start_third_alignment_read)+"-"+str(end_third_alignment_read))
                    
                    start_junction = int(junction_start) + int(junction_block_length_split[0])
                    end_junction = int(junction_stop) - int(junction_block_length_split[1]) + 1
                    
                    # If the junction is in the first alignment (1 and 2)
                    if end_first_alignment_read == start_junction and start_second_alignment_read == end_junction:
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                        
                    elif end_first_alignment_read == start_junction and start_second_alignment_read+1 == end_junction:
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                    
                    # If the junction is in the second alignment (2 and 3)
                    elif end_second_alignment_read == start_junction and start_third_alignment_read == end_junction:
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                        
                    elif end_second_alignment_read == start_junction and start_third_alignment_read+1 == end_junction:
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                        
                    #If the junction is in the third alignment (3 and 4)
                    elif end_third_alignment_read == start_junction and start_fourth_alignment_read == end_junction:
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                        
                    elif end_third_alignment_read == start_junction and start_fourth_alignment_read+1 == end_junction:
                        add_correct_reads(junction_count_read, junction_id, read_id)
                        nb_read_ok += 1
                        output.write(line_alignment)
                    else:
                        pass
    print(">>> Number of reads correctly mapped:", nb_read_ok)
#    print(junction_count_read.keys())
    return junction_count_read
    
def compare_read_align_against_list_junction_reference(dict_junction_reads, junction_file):
    '''
    FUNCTION DESCRIPTION:
        This function count the number of reads align for each junction.
        
    INPUT:
        - dict_junction_reads: a dictionary with all reads mapped correctly 
                               for each junction.
        - junction_file: the junctions file in bed format.
    '''
    with open(junction_file, 'r') as junction_file_i:
        lecture = junction_file_i.readlines()
        print(">>> READS aligned:")
        for junction in lecture:
            junction = junction.split('\t')
            try:
                print(junction[3], len(dict_junction_reads[junction[3]]))
            except KeyError:
                print(junction[3], "0")

def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="A SCRIPT TO OBTAIN SPECIFIC EXON JUNCTION NOT KNOWN IN ENSEMBL.")
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-a", "--align",
                        dest = "align_file", 
                        metavar = "ALIGNMENT FILE OBTAIN WITH BEDTOOLS FOR EXAMPLE", 
                        help = "Give the file with the output alignment obtain with BedTools intersect for example, in BED12 format.",
                        required = True)
    
    parser.add_argument("-j", "--junction",
                        dest = "junction_file", 
                        metavar = "EXON-EXON JUNCTIONS FILE", 
                        help = "Give the file with the exon-exon junctions in BED12 format.",
                        required = False)
    
    
    return parser.parse_args()       
    
if __name__ == "__main__":
    
    OPTIONS = parse()
    
    alignment_file_input = OPTIONS.align_file
    junction_file_input = OPTIONS.junction_file
    
    
    junction_reads = read_file(alignment_file_input)
    
    if junction_file_input:    
        print(">>> EXON JUNCTIONS FILE: YES")
        compare_read_align_against_list_junction_reference(junction_reads, junction_file_input)
    else:
        print(">>> EXON JUNCTIONS FILE: NO")
    
    
    