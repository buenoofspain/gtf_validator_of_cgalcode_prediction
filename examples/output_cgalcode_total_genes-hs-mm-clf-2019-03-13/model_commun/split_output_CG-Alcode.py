#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 09:43:14 2018

@author: Nicolas Guillaudeux
@team: DYLISS
@Supervisor: O. Dameron, C. Belleannée, S. Blanquart
"""

import argparse

def extract_CGAlcode_output(output_CGA_file):
    '''
    Function description: 
        This function split CGAlcode programm output to return three output 
        file : 
            - hs-mm comparison
            - hs-clf comparison
            - mm-clf comparison
            
    @INPUT:
            - output_CGA_file : The CGAlcode programm output within three 
              compararison
    
    @OUTPUT: 
            - EXPORT_HSMM.txt: comparison between human and mouse
            - EXPORT_HSCLF.txt: comparison between human and dog
            - EXPORT_MMCLF.txt: comparison between mouse and dog
    '''
    with open(output_CGA_file, 'r') as input_file:
        with open("EXPORT_HSMM.txt", 'w') as hsmm_o:
            with open("EXPORT_HSCLF.txt", 'w') as hsclf_o:
                with open("EXPORT_MMCLF.txt", 'w') as mmclf_o:
                    
                    lecture_file_= input_file.readlines()
                    count = 0
                    # Count: each comparison begin by OP   0, with this symbol 
                    # it is possible to separate comparison with a counter
                    for line in lecture_file_:
                        if line.startswith("OP     0 ENS") :
                            count += 1
                        # HS-MM comparison:
                        if count == 1:
                            hsmm_o.write(line)
                        # HS-CLF comparison:
                        if count == 2:
                            hsclf_o.write(line)
                        # MM-CLF comparison:
                        if count == 3:
                            mmclf_o.write(line)
                        
if __name__ == "__main__":
#    path = "/home/niguilla/Documents/test_donnees_stat_cliques_doublons_singleton_181011/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/test-hs-mm-clf-2018-10-11/"
#    file_CGA = path+"output_run2.txt"
    
    # Parser
    parser = argparse.ArgumentParser(prog = "extract_corresponding_triplets_with_same_letter_number.py")
    
    # Output required parameters
    parser.add_argument("-o", "--outputcga", dest = "outputcga", metavar = "CGA-Alcode output", help = "CG-Alcode output within the three comparison: HS-MM, HS-CLF & MM-CLF", required = True)

    # Args definition
    args = parser.parse_args()
    
    CGA_output = args.outputcga

    # Function execution    
    extract_CGAlcode_output(CGA_output)
