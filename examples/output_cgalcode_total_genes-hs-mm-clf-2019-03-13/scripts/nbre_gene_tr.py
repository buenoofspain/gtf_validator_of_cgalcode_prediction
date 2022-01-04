#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 14:34:58 2018

@author: niguilla
"""

import os

def test_number_genes_transcrips(species):
    '''
    This function count the number of transcripts found in genes/transcripts
    files correspondance for a species tested.
    '''
    
    species_tested = species
    sp_tested = ""
    
    if species_tested == "human":
        sp_tested = "hs"
    elif species_tested == "mouse":
        sp_tested = "mm"
    elif species_tested == "dog":
        sp_tested = "clf"
    
    rep_sp = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/"+species_tested+"/"
    genes_species_list = list()
    tr_species_list = list()
    
    with open("/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/total_genes_with_error_to_ytbfcn-hs-mm-clf-2018-12-18/transcripts_with_prediction/transcripts_"+sp_tested+"_withprediction.csv", 'r') as species_file_i:
        lecture_species = species_file_i.readlines()
        for lines_species in lecture_species:
            lines_species = lines_species.replace('\n', '').split(",")
            if lines_species[0] not in genes_species_list and lines_species[0] != "gene_ID":
                genes_species_list.append(lines_species[0])
            if lines_species[1] not in tr_species_list and lines_species[1] != "transcript_ID":
                tr_species_list.append(lines_species[1])
    
    nb_tr = 0
    nb_tr_abs = 0
    
    nb_g = 0
    nb_g_abs = 0
    
    for g in genes_species_list:
        if os.path.exists(rep_sp+g+".fasta"):
            nb_g += 1
        else:
            nb_g_abs += 1
    
    for tr in tr_species_list:
        if os.path.exists(rep_sp+tr+".fasta"):
            nb_tr += 1
        else:
            nb_tr_abs += 1
            
    print("Nombre de gènes dans le filtre :", len(genes_species_list))
    print("Nombre de transcrits dans le filtre :", len(tr_species_list))
    
    print("Nombre de gènes :", nb_g)
    print("Nombre de gènes absents :", nb_g_abs)
    
    print("Nombre de transcrits :", nb_tr)
    print("Nombre de transcrits absents :", nb_tr_abs) 
    
if __name__ == "__main__":
        print("HUMAN")
        test_number_genes_transcrips("human")
        print("******\nMOUSE")
        test_number_genes_transcrips("mouse")
        print("******\nDOG")
        test_number_genes_transcrips("dog")