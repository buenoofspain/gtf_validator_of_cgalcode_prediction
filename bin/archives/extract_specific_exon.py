#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 16:39:51 2019

@author: Nicolas Guillaudeux
"""

import argparse

def extract_reference_position(ref_gtf_file):
    """
    Function description:
        This function creates a list of reference positions from a gtf file.
    
    Input:
        - ref_gtf_file: a file who contain all exon information.
    
    Return:
        - positions_ref_list: a list which contains all exon positions (start
                              and stop) confused.
    """
    positions_ref_list = list()
    
    with open(ref_gtf_file, 'r') as ref_gtf_f:
        lecture = ref_gtf_f.readlines()
        for line in lecture:
            line = line.replace('\n', '').split('\t')
            pos_exon_start = line[1]
            pos_exon_stop = line[2]
            if not pos_exon_start in positions_ref_list:
                positions_ref_list.append(pos_exon_start)
            if not pos_exon_stop in positions_ref_list:
                positions_ref_list.append(pos_exon_stop)
    
    return positions_ref_list
    

def read_files(pred_gtf_file, ref_gtf_file):
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
    # Extraire les positions de chaque exons dans une liste
    positions_ref_list = extract_reference_position(ref_gtf_file)
#    print(positions_ref_list)
    # vérifier que chaque position n'est pas dans la liste de ref
    with open("output.exonsSpecificPrediction.bed", 'w') as output_f:
        with open(pred_gtf_file, 'r') as pred_gtf_f:
            lecture = pred_gtf_f.readlines()
            for line in lecture:
                line_exon = line.replace('\n', '').split("\t")
                pos_exon_start = line_exon[1]
                pos_exon_stop = line_exon[2]
#                print(pos_exon_start, pos_exon_stop)
#                print(not pos_exon_start in positions_ref_list, not pos_exon_stop in positions_ref_list)
                # Si aucune des deux, réécrire la ligne dans la sortie
                if not pos_exon_start in positions_ref_list and not pos_exon_stop in positions_ref_list:
                    print("SORTIE DIFFERENTE A 100%")
                    print(pos_exon_start, pos_exon_stop)
#                    output_f.write(line)
                elif not pos_exon_start in positions_ref_list or not pos_exon_stop in positions_ref_list:
                    print("SORTIE DIFFERENTE A 1 EXTREMITE")
                    print(pos_exon_start, pos_exon_stop)
                
                # Si un ou les deux sont connues, ne pas conserver pour le moment
    return
        
def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="A SCRIPT TO OBTAIN SPECIFIC EXON JUNCTION NOT KNOWN IN ENSEMBLE.")
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-a",
                        dest = "gtf_pred_file", 
                        metavar = "EXON GTF FILE FROM PREDICTION", 
                        help = "Give an exon gtf file with exon from prediction to search specific exon from prediction.")
    
    parser.add_argument("-b",
                        dest = "gtf_ref_file",
                        metavar = "EXON GTF FILE FROM REFERENCE",
                        help = "Give an exon gtf file with exon from reference to search specific exon from prediction.")
    
    
    return parser.parse_args()       
    
if __name__ == "__main__":
    
    OPTIONS = parse()
    
    pred_file = OPTIONS.gtf_pred_file
    ref_file = OPTIONS.gtf_ref_file
    
    read_files(pred_file, ref_file)
