#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 09:01:20 2018

@author: niguilla
"""

import os

def genes_tr_relation(gene_list, graph_analysis, tr_otr_relationship, sp1, sp2, sp3):
    """
    """
    # To define gene and its transcripts
    trGeneRelation = trGeneSeparator()
    
    # To analyse the genes triplet list
        
    with open("graph_tr_analysis_detail.csv", 'w') as output:
        
        with open("relation_gene_known_predict.csv", 'w') as output_analysis :
            output_analysis.write("hs_gene,hs_tr_known,hs_tr_pred,mm_gene,mm_tr_known,mm_tr_pred,clf_gene,clf_tr_known,clf_tr_pred\n")
    
            with open(graph_analysis, 'r') as graph_analyse:
                lecture = graph_analyse.readlines()
                for line in lecture:
                    if not line.startswith("triplet_number"):
                        line = line.replace("\n", "").split(',')
#                        print(line)
                        hs_gene = line[1]
                        mm_gene = line[2]
                        clf_gene = line[3]
                        
                        hs_tr = 0
                        hs_pr = 0
                        hs_mm_P = 0
                        hs_clf_P = 0
                        
                        mm_tr = 0
                        mm_pr = 0
                        mm_hs_P = 0
                        mm_clf_P = 0
                        
                        clf_tr = 0
                        clf_pr = 0
                        clf_hs_P = 0
                        clf_mm_P = 0
                        
                        for tr_hs in trGeneRelation[hs_gene]:
                            if tr_hs.startswith("CGATENSMUST"):
#                                print(os.path.isfile(human+tr_hs+".fasta"))
                                if os.path.isfile(human+tr_hs+".fasta"):
                                    hs_mm_P += 1    
                                    hs_pr += 1
                            if tr_hs.startswith("CGATENSCAFT"):
                                if os.path.isfile(human+tr_hs+".fasta"):
                                    hs_clf_P += 1
                                    hs_pr += 1
                            if not tr_hs.startswith("CGA"):
                                hs_tr += 1
                                
                        for tr_mm in trGeneRelation[mm_gene]:
                            if tr_mm.startswith("CGAMUSTENST"):
                                if os.path.isfile(mouse+tr_mm+".fasta"):
                                    mm_hs_P += 1
                                    mm_pr += 1
                            if tr_mm.startswith("CGAMUSTENSCAFT"):
                                if os.path.isfile(mouse+tr_mm+".fasta"):
                                    mm_clf_P += 1
                                    mm_pr += 1
                            if not tr_mm.startswith("CGA"):
                                mm_tr += 1
                                
                        for tr_clf in trGeneRelation[clf_gene]:
                            if tr_clf.startswith("CGACAFTENST"):
                                if os.path.isfile(dog+tr_clf+".fasta"):
                                    clf_hs_P += 1
                                    clf_pr += 1
                            if tr_clf.startswith("CGACAFTENSMUST"):
                                if os.path.isfile(dog+tr_clf+".fasta"):
                                    clf_mm_P += 1
                                    clf_pr += 1
                            if not tr_clf.startswith("CGA"):
                                clf_tr += 1
                                
                        wr = hs_gene+","+str(hs_tr)+" + "+str(hs_mm_P)+" mm + "+str(hs_clf_P)+" clf,"+mm_gene+","+str(mm_tr)+" + "+str(mm_hs_P)+" hs + "+str(mm_clf_P)+" clf,"+clf_gene+","+str(clf_tr)+" + "+str(clf_hs_P)+" hs + "+str(clf_mm_P)+" mm"+"\n"#","+clique_nb+","+tot_nb+","+model+"\n"
#                        print(hs_gene+","+str(hs_tr)+" + "+str(hs_mm_P)+" mm + "+str(hs_clf_P)+" clf,"+mm_gene+","+str(mm_tr)+" + "+str(mm_hs_P)+" hs + "+str(mm_clf_P)+" clf,"+clf_gene+","+str(clf_tr)+" + "+str(clf_hs_P)+" hs + "+str(clf_mm_P)+" mm")
                        output.write(str(wr))
                        
                        wr2 = hs_gene +','+ str(hs_tr) +','+ str(hs_pr) +','+ mm_gene +','+ str(mm_tr) +','+ str(mm_pr) +','+ clf_gene +','+ str(clf_tr) +','+ str(clf_pr) +'\n'
#                        print(wr2)
                        output_analysis.write(wr2)
    
    sort_tr_otr_relationship = sortTranscriptsDir(tr_otr_relationship, human, mouse, dog)
    
#    with open(tr_otr_relationship, 'r') as tr_otr_i:
#        lecture = tr_otr_i.readlines()
#        for line in lecture:
#            print(line.replace("\n", ""))
            
    
    
def sortTranscriptsDir(triplet_prediction_file, dataX, dataY, dataZ):
    """
    Function description:
        This function sort transcripts to remove redundant transcripts.
        
    Parameters:
    - triplet_prediction_file: a file with transcrips and prediction with a 
      redundant possibility. 
      File content: "GeneSp1, GeneSp2, TrSp1, TrSp2, Relation"
    - dataX, dataY, dataZ: Directories of species data.
    
    Output:
    - Remove predicted transcripts files with redondance proved.
    """
    
    # Dictionary of relation correspondance
    descr_tr_relation = dict()
    # List of discovered redundancies
    comparison_tr = list()
    
    descr_tr_relation["Found"] = []
    descr_tr_relation["Found-CN"] = []
    descr_tr_relation["Yet-to-be-found"] = []
    descr_tr_relation["Yet-to-be-found-CN"] = []
    descr_tr_relation["No-executable"] = []
    
    # Sorting data with a dictionary
    with open(triplet_prediction_file, 'r') as verif_file:
        lecture = verif_file.readlines()
        for line in lecture:
            line = line.replace('\n', '').split(',')
            descr_tr_relation[line[4]].append([line[0], line[1], line[2], line[3]])
        
#    # Number of redundancies
##    count_tr = 0
#    
    # Searching ortholog correspondance between the three species
    for pairwise1 in descr_tr_relation['Found']:
        tr1 = pairwise1[2]
        tr2 = pairwise1[3]
        g1 = pairwise1[0]
        g2 = pairwise1[1]
        g3 = ""
        trCGA1, trCGA2 = '', ''
        for pairwise2 in descr_tr_relation["Yet-to-be-found"]:
            if tr1 == pairwise2[2] and pairwise2[3].startswith('CGA') and trCGA1 == '':
                trCGA1 = pairwise2[3]
                g3 = pairwise2[1]
            if tr2 == pairwise2[2] and pairwise2[3].startswith('CGA') and trCGA2 == '':
                trCGA2 = pairwise2[3]
                g3 = pairwise2[1]
        for pairwise3 in descr_tr_relation["Yet-to-be-found-CN"]:
            if tr1 == pairwise3[2] and pairwise3[3].startswith('CGA') and trCGA1 == '':
                trCGA1 = pairwise3[3]
                g3 = pairwise3[1]
            if tr2 == pairwise3[2] and pairwise3[3].startswith('CGA') and trCGA2 == '':
                trCGA2 = pairwise3[3]
                g3 = pairwise3[1]
        if (g1 and g2 and g3 and tr1 and tr2 and trCGA1 and trCGA2) != "" and (g1, g2, g3, tr1, tr2, trCGA1, trCGA2) not in comparison_tr and (g2, g1, g3, tr2, tr1, trCGA2, trCGA1) not in comparison_tr:
            comparison_tr.append((g1, g2, g3, tr1, tr2, trCGA1, trCGA2))
#            count_tr += 1
        
    for pairwise1 in descr_tr_relation['Found-CN']:
        tr1 = pairwise1[2]
        tr2 = pairwise1[3]
        g1 = pairwise1[0]
        g2 = pairwise1[1]
        trCGA1, trCGA2 = '', ''
        for pairwise2 in descr_tr_relation["Yet-to-be-found"]:
            if tr1 == pairwise2[2] and pairwise2[3].startswith('CGA') and trCGA1 == '':
                trCGA1 = pairwise2[3]
                g3 = pairwise2[1]
            if tr2 == pairwise2[2] and pairwise2[3].startswith('CGA') and trCGA2 == '':
                trCGA2 = pairwise2[3]
                g3 = pairwise2[1]
        for pairwise3 in descr_tr_relation["Yet-to-be-found-CN"]:
            if tr1 == pairwise3[2] and pairwise3[3].startswith('CGA') and trCGA1 == '':
                trCGA1 = pairwise3[3]
                g3 = pairwise3[1]
            if tr2 == pairwise3[2] and pairwise3[3].startswith('CGA') and trCGA2 == '':
                trCGA2 = pairwise3[3]
                g3 = pairwise3[1]
        if (g1 and g2 and g3 and tr1 and tr2 and trCGA1 and trCGA2) != "" and (g1, g2, g3, tr1, tr2, trCGA1, trCGA2) not in comparison_tr and (g2, g1, g3, tr2, tr1, trCGA2, trCGA1) not in comparison_tr:
            comparison_tr.append((g1, g2, g3, tr1, tr2, trCGA1, trCGA2))
#            count_tr += 1
    
#    print("Number of redundancies discovered", count_tr)
            
    # Deleting redundancies
#    removeRedundantTranscripts(comparison_tr, dataX, dataY, dataZ)
    return descr_tr_relation
    
def trGeneSeparator():
    """
    """
    benchmark = os.getcwd() + "/transcripts_with_prediction"
    dog = benchmark+"/transcripts_clf_withprediction.csv"
    human = benchmark+"/transcripts_hs_withprediction.csv"
    mouse = benchmark+"/transcripts_mm_withprediction.csv"
    
    dict_gene_tr_relation = dict()
    
    with open(dog, 'r') as dog_input:
        lecture = dog_input.readlines()
        for line in lecture:
            if not line.startswith("gene_ID"):
                line = line.replace("\n", "").split(",")
                try:
                    if line[1] not in dict_gene_tr_relation[line[0]]:
                        dict_gene_tr_relation[line[0]].append(line[1])
                except KeyError:
                    dict_gene_tr_relation[line[0]] = [line[1]]
                    
    with open(human, 'r') as dog_input:
        lecture = dog_input.readlines()
        for line in lecture:
            if not line.startswith("gene_ID"):
                line = line.replace("\n", "").split(",")
                try:
                    if line[1] not in dict_gene_tr_relation[line[0]]:
                        dict_gene_tr_relation[line[0]].append(line[1])
                except KeyError:
                    dict_gene_tr_relation[line[0]] = [line[1]]
                    
    with open(mouse, 'r') as dog_input:
        lecture = dog_input.readlines()
        for line in lecture:
            if not line.startswith("gene_ID"):
                line = line.replace("\n", "").split(",")
                try:
                    if line[1] not in dict_gene_tr_relation[line[0]]:
                        dict_gene_tr_relation[line[0]].append(line[1])
                except KeyError:
                    dict_gene_tr_relation[line[0]] = [line[1]]

    return dict_gene_tr_relation
    
if __name__ == "__main__":
    
    # Files:
    benchmark = os.getcwd()
    genes_list = benchmark + "/data_set_interessant/112_2.txt"
    graph_analysis = benchmark + "/data_set_interessant/graph_tr_analysis_112.csv"
    tr_otr_relationship = benchmark + "/corresponding_tr_otr.txt"
    dog = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/dog/"
    mouse = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/mouse/"
    human = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/"
    
    genes_tr_relation(genes_list, graph_analysis, tr_otr_relationship, human, mouse, dog)
