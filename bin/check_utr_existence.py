#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:00:21 2020

@author: niguilla
"""

import os
import networkx as nx

def list_tr(file_tr_list):
    """
    This function uses a file with all transcript list and return the 
    transcripts in a list format data.
    """
    tr_list_known = list()
    with open(file_tr_list, 'r') as tr_list_i:
        lecture = tr_list_i.readlines()
        for tr in lecture:
            tr = tr.replace('\n', '')
            if not tr in tr_list_known:
                tr_list_known.append(tr)
    
    return tr_list_known
        

def check_files_utr(repository, tr_list):
    """
    This function lists all transcripts without utr at least in one side.
    """
    tr_without_utr = list()
    
    for file in os.listdir(repository):
        line_nb = 0
        tr_info = ""
        utr_left = True
        utr_right = True
        tr_id = file[:-6]

        if tr_id in tr_list:
            with open(repository+file, 'r') as f_i:
                lecture = f_i.readlines()
                for line in lecture:
                    line_nb += 1
                    line = line.replace('\n', '')
                    if line_nb == 1:
                        line = line.split(' ')
                        transcript_id = line[1]
                        gene_id = line[2]
                        transcript_name = line[7]
                        
                    if line_nb == 3 and line == '[]':
                        utr_left = False
                    if line_nb == 4 and line == '[]':
                        utr_right = False
                    
                ############################################
                # Remplacer or par and selon ce qu'on veut #
                if utr_left == False or utr_right == False:
                    tr_without_utr.append(transcript_id)#+","+str(utr_left)+","+str(utr_right))

    return tr_without_utr

def check_files_utr_side(repository, tr_list):
    """
    This function lists all transcripts without utr at least in one side.
    """
    tr_without_utr_left_right = list()
    tr_without_utr_left = list()
    tr_without_utr_right = list()
    
    for file in os.listdir(repository):
        line_nb = 0
        tr_info = ""
        utr_left = True
        utr_right = True
        tr_id = file[:-6]

        if tr_id in tr_list:
            with open(repository+file, 'r') as f_i:
                lecture = f_i.readlines()
                for line in lecture:
                    line_nb += 1
                    line = line.replace('\n', '')
                    if line_nb == 1:
                        line = line.split(' ')
                        transcript_id = line[1]
                        gene_id = line[2]
                        transcript_name = line[7]
                        
                    if line_nb == 3 and line == '[]':
                        utr_left = False
                    if line_nb == 4 and line == '[]':
                        utr_right = False
                    
                ############################################
                # Remplacer or par and selon ce qu'on veut #
                if utr_left == False and utr_right == False:
                    tr_without_utr_left_right.append(transcript_id)
                else:
                    if utr_left == False :
                        tr_without_utr_left.append(transcript_id)
                    if utr_right == False:
                        tr_without_utr_right.append(transcript_id)#+","+str(utr_left)+","+str(utr_right))

    return tr_without_utr_left_right, tr_without_utr_left, tr_without_utr_right

def check_set_gene(list_tr_without_utr, list_utr_interest):
    """
    """
    tr_set = list()
    with open(list_utr_interest, 'r') as f_i:
        lecture = f_i.readlines()
        for tr in lecture:
            tr = tr.replace('\n', '')
            tr_set.append(tr)

    print("Number of transcripts known in the set of 118:", len(tr_set))    
    
    tr_without_utr_nb = 0
    for tr in tr_set:
        if tr in list_tr_without_utr:
            tr_without_utr_nb += 1
#            print(tr)
    print("Number of transcripts without utr:", tr_without_utr_nb)
    
def gene_tr_relationship(gene_tr_file):
    """
    This function uses a file with the gene - transcript correspondance and
    return a dictionary with this relation.
    """
    gene_tr_dict = dict()
    
    with open(gene_tr_file, 'r') as f_i:
        lecture = f_i.readlines()
        for line in lecture:
            line = line.replace('\n', '').split(',')
            gene, tr = line[0], line[1]
            if not gene in gene_tr_dict:
                gene_tr_dict[gene] = [tr]
            else:
                gene_tr_dict[gene].append(tr)
            
    return gene_tr_dict
    
def check_set_gene_tr_if_utr(gene_tr_relation, tr_list_without_utr, set_gene_list, species):
    """
    This function tests if a gene has all of these transcripts with UTR.
    """
    gene_conserved = list()
    
    with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/135/"+species+"_gene_with_all_utr.txt", 'w') as f_o:
        for gene in set_gene_list:
            all_tr = True
            for transcript in gene_tr_relation[gene]:
                if not transcript.startswith("CGA"):
                    if transcript in tr_list_without_utr:
                        all_tr = False
                        print(transcript, "not UTR in", gene)
                        
            if all_tr == False:
                print(gene+" not conserved")
            else:
                gene_conserved.append(gene)
                f_o.write(gene+'\n')
                
#    print("ENSG00000001631" in gene_conserved)
        
    print("We conserved:", len(gene_conserved), "genes in", species)
    
    return gene_conserved
   
def check_set_gene_tr_if_utr_test(gene_tr_relation, tr_list_without_utr, set_gene_list, species):
    """
    This function tests if a gene has all of these transcripts with UTR.
    """
    gene_conserved = list()
    
    with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/135/"+species+"_gene_with_all_utr.txt", 'w') as f_o:
        for gene in set_gene_list:
#            all_tr = True
            all_tr_count = len(gene_tr_relation[gene]) 
            all_tr = len(gene_tr_relation[gene])
            for transcript in gene_tr_relation[gene]:
                if not transcript.startswith("CGA"):
                    if transcript in tr_list_without_utr:
                        all_tr -= 1
                        print(transcript, "not UTR in", gene)
                else:
                    all_tr -= 1
                        
            if all_tr < all_tr_count:
                print(gene+" not conserved")
            else:
                gene_conserved.append(gene)
                f_o.write(gene+'\n')
        
    print("We conserved:", len(gene_conserved), "genes in", species)
    
    return gene_conserved
    
if __name__ == "__main__":
    
    ###########################################################################
    ### >>> HUMAN TEST <<<
    print("Test in Human species:")
    repository = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/"
    list_tr_known = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/tr_list_hs_tot.sort.known.txt"
    
    # Create a list with all transcript known:
    list_tr_hs = list_tr(list_tr_known)
#    print(list_tr)
    
    # Check if transcript are utr:
    tr_without_utr_list_hs = check_files_utr(repository, list_tr_hs)
#    tr_without_utr_list_hs_l_r, tr_without_utr_list_hs_l, tr_without_utr_list_hs_r = check_files_utr_side(repository, list_tr_hs)
#    print(tr_without_utr_list_hs_l_r)
#    print(tr_without_utr_list_hs_l)
#    print(tr_without_utr_list_hs_r)
    
    
    # Check if set of transcripts are utr:
    #tr_set = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/transcript_known_set_253_hs.txt"
    tr_set = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/135/transcript_list_118_hs.txt"
    check_set_gene(tr_without_utr_list_hs, tr_set)
    
    # Give gene/transcript relation:
    gene_tr_hs_file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/transcripts_hs_withprediction_final.csv"
    gene_tr_dict_hs = gene_tr_relationship(gene_tr_hs_file)
    
    #list 118:
    hs_118 = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/135/118_gene_list_hs.txt"
    list_hs_118 = list_tr(hs_118)
    
    # Test if UTR in all transcript of 118:
    gene_hs_list_conserved = check_set_gene_tr_if_utr(gene_tr_dict_hs, tr_without_utr_list_hs, list_hs_118, "hs")
    
    ###########################################################################
    ### >>> MOUSE TEST <<<
    print("")
    print("Test in Mouse species:")
    repository = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/mouse/"
    list_tr_known = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/tr_list_mm_tot.sort.known.txt"
    
    # Create a list with all transcript known:
    list_tr_mm = list_tr(list_tr_known)
#    print(list_tr)
    
    # Check if transcript are utr:
    tr_without_utr_list_mm = check_files_utr(repository, list_tr_mm)
#    tr_without_utr_list_mm_l_r, tr_without_utr_list_mm_l, tr_without_utr_list_mm_r = check_files_utr_side(repository, list_tr_mm)
#    print(tr_without_utr_list_mm_l_r)
#    print(tr_without_utr_list_mm_l)
#    print(tr_without_utr_list_mm_r)
    
    # Check if set of transcripts are utr:
    #tr_set = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/transcript_known_set_253_mm.txt"
    tr_set = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/135/transcript_list_118_mm.txt"
    check_set_gene(tr_without_utr_list_mm, tr_set)
    
    # Give gene/transcript relation:
    gene_tr_mm_file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/transcripts_mm_withprediction_final.csv"
    gene_tr_dict_mm = gene_tr_relationship(gene_tr_mm_file)
    
    #list 118:
    mm_118 = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/135/118_gene_list_mm.txt"
    list_mm_118 = list_tr(mm_118)
    
    # Test if UTR in all transcript of 118:
    gene_mm_list_conserved = check_set_gene_tr_if_utr(gene_tr_dict_mm, tr_without_utr_list_mm, list_mm_118, "mm")
    
    ###########################################################################
    ### >>> DOG TEST <<<
    print("")
    print("Test in Dog species:")
    repository = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/dog/"
    list_tr_known = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/tr_list_clf_tot.sort.known.txt"
    
    # Create a list with all transcript known:
    list_tr_clf = list_tr(list_tr_known)
#    print(list_tr)
    
    # Check if transcript are utr:
    tr_without_utr_list_clf = check_files_utr(repository, list_tr_clf)
    tr_without_utr_list_clf_l_r, tr_without_utr_list_clf_l, tr_without_utr_list_clf_r = check_files_utr_side(repository, list_tr_clf)
#    print(tr_without_utr_list_clf_l_r)
#    print(tr_without_utr_list_clf_l)
#    print(tr_without_utr_list_clf_r)
    
    # Check if set of transcripts are utr:
    #tr_set = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/transcript_known_set_253_clf.txt"
    tr_set = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/135/transcript_list_118_clf.txt"
    check_set_gene(tr_without_utr_list_clf, tr_set)
    
    # Give gene/transcript relation:
    gene_tr_clf_file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/transcripts_clf_withprediction_final.csv"
    gene_tr_dict_clf = gene_tr_relationship(gene_tr_clf_file)
    
    #list 118:
    clf_118 = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/135/118_gene_list_clf.txt"
    list_clf_118 = list_tr(clf_118)
    
    # Test if UTR in all transcript of 118:
    gene_clf_list_conserved = check_set_gene_tr_if_utr(gene_tr_dict_clf, tr_without_utr_list_clf, list_clf_118, "clf")
    
    
    orthologous_gene = list()
    utr_ok_gene = 0
    genes_conserved = list()
    with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/135/resume_gene_information_human_mouse_dog.csv", 'r') as f_i:
        lecture = f_i.readlines()
        for line in lecture:
            if not line.startswith("human_ensembl"):
                line = line.replace('\n', '').split('\t')
                hs_gene, mm_gene, clf_gene = line[0], line[3], line[6]
                if not [hs_gene, mm_gene, clf_gene] in orthologous_gene:
                    orthologous_gene.append([hs_gene, mm_gene, clf_gene])
#                    print(orthologous_gene)
#                    print(hs_gene, mm_gene, clf_gene)
                    if hs_gene in gene_hs_list_conserved and mm_gene in gene_mm_list_conserved and clf_gene in gene_clf_list_conserved:
                        utr_ok_gene += 1
                        print(hs_gene, mm_gene, clf_gene)
                        genes_conserved.append((hs_gene, mm_gene, clf_gene))
    print("We conserved", utr_ok_gene, "genes")
    

    ###########################################################################
    # Check transcript graph to test alternative transcription:
    
    #list_hs_118 : liste des 118 gènes
    #gene_tr_dict_hs : liste des relations gènes/tr pour tous les gèness (2 167)
    
#    for gene_conserved in genes_conserved:
#        
#        hs_gene = gene_conserved[0]
#        mm_gene = gene_conserved[1]
#        clf_gene = gene_conserved[2]
#        print(hs_gene, mm_gene, clf_gene)
#        print("hs_tr", len(gene_tr_dict_hs[hs_gene]), "mm_tr:", len(gene_tr_dict_mm[mm_gene]), "clf_tr:", len(gene_tr_dict_clf[clf_gene]))
#        
#        print(gene_tr_dict_hs[hs_gene])
#        print(gene_tr_dict_mm[mm_gene])
#        print(gene_tr_dict_clf[clf_gene])
#        
#        for tr in gene_tr_dict_hs[hs_gene]:
#            if tr in tr_without_utr_list_hs:
#                print("TEST", tr in tr_without_utr_list_hs)
#           
#        for tr in gene_tr_dict_mm[mm_gene]:
#            if tr in tr_without_utr_list_mm:
#                print("TEST", tr in tr_without_utr_list_mm)
#        
#        for tr in gene_tr_dict_clf[clf_gene]:
#            if tr in tr_without_utr_list_clf:
#                print("TEST", tr in tr_without_utr_list_clf)