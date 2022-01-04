#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 11:14:24 2018

@author: niguilla

This script analyses the gene triplet with a common model and if they are only 
contains clique3 in graph delimitor
"""
import argparse
import ast

from extract_corresponding_genes_triplets import extractGOWithinOutputFile, dictCompareModels
from caracteristics import caracteristicsInCSV

def genesTripletTable(cliquesOutput, commonComparisonOutput):
    """
    Functiopn description:
        This function give a table in csv format with genes triplets with same 
        gene structural model GO.
        
    Input:
    - cliquesOutput: this file contains the detail in intron delimitor in 
        subgraph format (clique3 != 0, doublon = 0, singleton = 0) for each 
        genes triplet.
    - commonComparisonOutput: this file contains the genes triplet with the
        same gene structural model GO.
    
    Output:
    - cliques_common.csv: a file with all the information of the two files 
        combined.
    """
    with open("cliques_common.csv", 'w') as output:
        title = "HUMAN GENE,HUMAN TRANSCRIPTS,MOUSE GENE,MOUSE TRANSCRIPTS,DOG GENE,DOG TRANSCRIPTS,CLIQUES NUMBER,TOTAL SIGNAL,COMMON MODEL\n"
        output.write(title)
        with open(commonComparisonOutput, 'r') as common:
            with open(cliquesOutput, 'r') as clique:
                lecture_clique = clique.readlines()
                lecture_common = common.readlines()
                for line_clique in lecture_clique:
                    if not line_clique.startswith("triplet_number"):
                        line_clique = line_clique.replace("\n", "").split(",")
                        triplet = line_clique[1]
                        print(triplet)
                        clique_nb = line_clique[2]
                        tot_nb = line_clique[9]
                        for line_common in lecture_common:
                            if not line_common.startswith("GENE NUMBER"):
                                line_common = line_common.replace("\n", "").split(",")
                                g1 = line_common[1]
                                tr1 = line_common[2]
                                g2 = line_common[3]
                                tr2 = line_common[4]
                                g3 = line_common[5]
                                tr3 = line_common[6]
                                model = line_common[7]
                                if g1 in triplet:
                                    hs_mm_P = 0
                                    hs_clf_P = 0
                                    mm_hs_P = 0
                                    mm_clf_P = 0
                                    clf_hs_P = 0
                                    clf_mm_P = 0
                                    trGeneRelation = trGeneSeparator()
                                    for tr_hs in trGeneRelation[g1]:
#                                        print(trGeneRelation[g1])
#                                        print(trGeneRelation[g2])
#                                        print(trGeneRelation[g3])
                                        if tr_hs.startswith("CGATENSMUST"):
                                            hs_mm_P += 1
                                        if tr_hs.startswith("CGATENSCAFT"):
                                            hs_clf_P += 1
                                    for tr_mm in trGeneRelation[g2]:
                                        if tr_mm.startswith("CGAMUSTENST"):
                                            mm_hs_P += 1
                                        if tr_mm.startswith("CGAMUSTENSCAFT"):
                                            mm_clf_P += 1
                                    for tr_clf in trGeneRelation[g3]:
                                        if tr_clf.startswith("CGACAFTENST"):
                                            clf_hs_P += 1
                                        if tr_clf.startswith("CGACAFTENSMUST"):
                                            clf_mm_P += 1
                                    print(True)
                                    wr = g1+","+tr1+" + "+str(hs_mm_P)+" mm + "+str(hs_clf_P)+" clf,"+g2+","+tr2+" + "+str(mm_hs_P)+" hs + "+str(mm_clf_P)+" clf,"+g3+","+tr3+" + "+str(clf_hs_P)+" hs + "+str(clf_mm_P)+" mm"+","+clique_nb+","+tot_nb+","+model+"\n"
                                    print(g1+","+tr1+" + "+str(hs_mm_P)+" mm + "+str(hs_clf_P)+" clf,"+g2+","+tr2+" + "+str(mm_hs_P)+" hs + "+str(mm_clf_P)+" clf,"+g3+","+tr3+" + "+str(clf_hs_P)+" hs + "+str(clf_mm_P)+" mm")
                                    output.write(str(wr))
                                    break
                            
                            
def trGeneSeparator():
    benchmark = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/model_commun/transcripts_with_prediction"
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
                            
                            
def genesTripletsTableDelimitor(delimitorDict, transcripts):
    """
    """
    with open(delimitorDict) as delimitors:
        lecture_delimitor = delimitors.readlines()
        for line_delimitor in lecture_delimitor:
            line_delimitor = line_delimitor.split("#####")
                              
            first_genes_comparison_str = line_delimitor[0]
            first_genes_comparison = str(first_genes_comparison_str[1:-1])
            first_genes_comparison = ast.literal_eval(first_genes_comparison)
            
            second_genes_comparison_str = line_delimitor[1]
            second_genes_comparison = str(second_genes_comparison_str[1:-1])
            second_genes_comparison = ast.literal_eval(second_genes_comparison)
            
            third_genes_comparison_str = line_delimitor[2]
            third_genes_comparison = str(third_genes_comparison_str[1:-1])
            third_genes_comparison = ast.literal_eval(third_genes_comparison)
            
            dict1 = first_genes_comparison[0]
            dict2 = first_genes_comparison[1]
            dict3 = second_genes_comparison[0]
            dict4 = second_genes_comparison[1]
            dict5 = third_genes_comparison[0]
            dict6 = third_genes_comparison[1]
            
            print(first_genes_comparison[2])

if __name__ == "__main__":
    
    # Parser
    parser = argparse.ArgumentParser(prog = "extract_corresponding_triplets_with_same_letter_number.py")
    
    # Output required parameters
    parser.add_argument("-o1", "--output_1", dest = "output_1", metavar = "CG_ALCODE_OUTPUT_1", help = "Output of CG-alcode programm", required = True)
    parser.add_argument("-o2", "--output_2", dest = "output_2", metavar = "CG_ALCODE_OUTPUT_2", help = "Output of CG-alcode programm", required = True)
    parser.add_argument("-o3", "--output_3", dest = "output_3", metavar = "CG_ALCODE_OUTPUT_3", help = "Output of CG-alcode programm", required = True)
    parser.add_argument("-g", "--graphsAnalysis", dest = "graphAnalysis", metavar = "GRAPH_ANALYSIS.CSV", help = "Output graph_analysis.csv of CG-alcode multi-species comparison programm", required = True)
    parser.add_argument("-c", "--cliques", dest = "cliques", metavar = "CLIQUES.TXT", help = "Output cliques.txt of CG-alcode multi-species comparison programm", required = True)
    parser.add_argument("-d", "--doublons", dest = "doublons", metavar = "DOUBLONS.TXT", help = "Output doublons.txt of CG-alcode multi-species comparison programm", required = True)
    parser.add_argument("-s", "--singletons", dest = "singletons", metavar = "SINGLETONS.TXT", help = "Output singletons.txt of CG-alcode multi-species comparison programm", required = True)
    parser.add_argument("-dd", "--delimitorDict", dest = "delimitor", metavar = "OUTPUT_DICTIONARIS.TXT", help = "Output output_dictionaries.txt of CG-alcode multi-species comparison programm", required = False)
    parser.add_argument("-tr", "--transcripts", dest = "transcripts", metavar = "CORRESPONDING_TR_OTR.TXT", help = "Output corresponding_tr_otr.txt of CG-alcode multi-species comparison programm", required = False)

    # Args definition
    args = parser.parse_args()
    
    output_1 = args.output_1
    output_2 = args.output_2
    output_3 = args.output_3
    
    graphAnalysis = args.graphAnalysis
    cliques = args.cliques
    doublons = args.doublons
    singletons = args.singletons
    
    delimitorDict = args.delimitor
    transcripts = args.transcripts
    
    # Step 1: to generate output_common_comparison.csv
    dict_1 = extractGOWithinOutputFile(output_1)
    dict_2 = extractGOWithinOutputFile(output_2)
    dict_3 = extractGOWithinOutputFile(output_3)
    
    genes_models_dict, interest_genes_models = dictCompareModels(dict_1, dict_2, dict_3)
    
    with open('output_common_comparison.csv', 'w') as filout:
        title="GENE NUMBER,HUMAN GENE,NUMBER OF HUMAN TRANSCRIPTS,MOUSE GENE,NUMBER OF MOUSE TRANSCRIPTS,DOG GENE,NUMBER OF DOG TRANSCRIPTS,COMMON MODEL\n"
        filout.write(title)
        for genes, model_list in genes_models_dict.items():
            for model in model_list:
                if model.startswith("."):
                    model = model[1:]
                line = genes+','+model
                filout.write(line+'\n')
        print("")
    
    # Step 2: to generate cliques_output.csv
    caracteristicsInCSV(graphAnalysis, cliques, doublons, singletons)
    
    # Step 3: to generate clique_common.csv
    genesTripletTable("cliques_output.csv", "output_common_comparison.csv")
    
    # Step 4: to generate....
#    genesTripletsTableDelimitor(delimitorDict, transcripts)
