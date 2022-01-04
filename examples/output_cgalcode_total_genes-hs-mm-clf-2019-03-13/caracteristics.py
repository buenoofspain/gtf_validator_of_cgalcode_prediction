#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 17:36:03 2018

@author: niguilla
"""

import argparse
import os

def caracteristicsInCSV(graph_analysis, cliques, doublons, singletons):

    #graph_analysis = "/home/Stage_M2/DATA/stastical_analysis-hs-mm-clf-2018-05-22/graph_analysis.csv"
    #graph_analysis = "/home/Stage_M2/DATA/export_files_to_statistical_analysis-hs-mm-clf-2018-05-29/graph_analysis.csv"
#    graph_analysis = "/home/Stage_M2/DATA/To conserve/final_run2-hs-mm-clf-2018-06-08/graph_analysis.csv"
    
    #link = "/home/Stage_M2/DATA/stastical_analysis-hs-mm-clf-2018-05-22/"
    #link = "/home/Stage_M2/DATA/export_files_to_statistical_analysis-hs-mm-clf-2018-05-29/"
#    link = "/home/Stage_M2/DATA/To conserve/final_run2-hs-mm-clf-2018-06-08/"
#    cliques = link+"clique.txt"
#    doublons = link+"doublons.txt"
#    singletons = link+"singletons.txt"
    
    cliques_list, doublons_list, singletons_list = list(), list(), list()
    singletons_list, singletons_list = list(), list()
    
    with open(cliques, 'r') as cliques_file:
        lecture = cliques_file.readlines()
        for line in lecture:
            line = line.replace("\n", "")
            cliques_list.append(line)
    
    with open(doublons, 'r') as doublons_file:
        lecture = doublons_file.readlines()
        for line in lecture:
            line = line.replace("\n", "")
            doublons_list.append(line)
    
    with open(singletons, 'r') as singletons_file:
        lecture = singletons_file.readlines()
        for line in lecture:
            line = line.replace("\n", "")
            singletons_list.append(line)
        
    test_d1, test_d2, test_d3, test_s1, test_s2, test_s3, test_c = 0, 0, 0, 0, 0, 0, 0
    with open("cliques_output.csv", 'w') as cliques_output:
        with open("doublons_non_hs_output.csv", 'w') as doublons_non_hs_output:
            with open("doublons_non_mm_output.csv", 'w') as doublons_non_mm_output:
                with open("doublons_non_clf_output.csv", 'w') as doublons_non_clf_output:
                    with open("singletons_hs_output.csv", 'w') as singletons_hs_output: 
                        with open("singletons_mm_output.csv", 'w') as singletons_mm_output: 
                            with open("singletons_clf_output.csv", 'w') as singletons_clf_output: 
                                with open(graph_analysis, 'r') as graph_file:
                                    lecture = graph_file.readlines()
                                    for line in lecture:
                                        if line.startswith("triplet_number"):
                                            cliques_output.write(line)
                                            doublons_non_hs_output.write(line)
                                            doublons_non_mm_output.write(line)
                                            doublons_non_clf_output.write(line)
                                            singletons_hs_output.write(line)
                                            singletons_mm_output.write(line)
                                            singletons_clf_output.write(line)
                                        else:
                                            line_s = line.split(',')
                                            if line_s[1][:15] in cliques_list:
                                                if line_s[3] == '0' and line_s[4] == '0' and line_s[5] == '0' and line_s[6] == '0' and line_s[7] == '0' and line_s[8] == '0':
                                                    cliques_output.write(line)
                                                    test_c += 1
                                            if line_s[1][:15] in doublons_list:
                                                if line_s[3] != '0' and line_s[4] == '0' and line_s[5] == '0':
                                                    #print line_s[3]
                                                    test_d1 += 1
                                                    doublons_non_hs_output.write(line)
                                                if line_s[3] == '0' and line_s[4] != '0' and line_s[5] == '0':
                                                    test_d2 += 1
                                                    doublons_non_mm_output.write(line)
                                                if line_s[3] == '0' and line_s[4] == '0' and line_s[5] != '0':
                                                    test_d3 += 1
                                                    doublons_non_clf_output.write(line)
                                            if line_s[1][:15] in singletons_list:
                                                if line_s[6] != '0' and line_s[7] == '0' and line_s[8] == '0':
                                                    #print line_s[3]
                                                    test_s1 += 1
                                                    singletons_hs_output.write(line)
                                                if line_s[6] == '0' and line_s[7] != '0' and line_s[8] == '0':
                                                    test_s2 += 1
                                                    singletons_mm_output.write(line)
                                                if line_s[6] == '0' and line_s[7] == '0' and line_s[8] != '0':
                                                    test_s3 += 1
                                                    singletons_clf_output.write(line)
                                                #singletons_output.write(line)
                            
    print "CARACTERISTICS", test_c, test_d1, test_d2, test_d3, test_s1, test_s2, test_s3
    '''
        - csv_file: the file with all caracteristics og genes triplet with 10 columns: 
            #triplet number -> 0
            #triplet genes name -> 1
            #number of clique -> 2
            #number of signal abs in human gene -> 3
            #number of signal abs in mouse gene -> 4
            #number of signal abs in dog gene -> 5
            #number of signal specific to human gene -> 6
            #number of signal specific to mouse gene -> 7
            #number of signal specific to dog gene -> 8
            #total number of signal -> 9
    '''
    with open(singletons, 'r') as singletons_file:
        lecture = singletons_file.readlines()
        for line in lecture:
            line = line.replace("\n", "")
            singletons_list.append(line)
                
    #with open(graph_analysis,'r') as compare_file:
    #        #lecture = compare_file.read().splitlines("#####")
    #        lines = compare_file.readlines()
    #        for line in lines:
    #            if not line.startswith("triplet_number"):
    #                print line
    #                line = line.split(',')
    #                triplet = line[1]
    #                genes_triplet = triplet.split("-")
    #                if genes_triplet[0] not in genes_list:
    #                    genes_list.append(genes_triplet[0])
    #                if genes_triplet[1] not in genes_list:
    #                    genes_list.append(genes_triplet[1])
    #                if genes_triplet[2] not in genes_list:
    #                    genes_list.append(genes_triplet[2])
    
if __name__ == "__main__":
    
    # Parser
    parser = argparse.ArgumentParser(prog = "extract_corresponding_triplets_with_same_letter_number.py")
    
    # Output required parameters
    parser.add_argument("-g", "--graphsAnalysis", dest = "graphAnalysis", metavar = "GRAPH ANALYSIS", help = "Output graph_analysis.csv of CG-alcode multi-species comparison programm", required = True)
    parser.add_argument("-c", "--cliques", dest = "cliques", metavar = "CLIQUES TXT", help = "Output cliques.txt of CG-alcode multi-species comparison programm", required = True)
    parser.add_argument("-d", "--doublons", dest = "doublons", metavar = "DOUBLONS TXT", help = "Output doublons.txt of CG-alcode multi-species comparison programm", required = True)
    parser.add_argument("-s", "--singletons", dest = "singletons", metavar = "SINGLETONS TXT", help = "Output singletons.txt of CG-alcode multi-species comparison programm", required = True)

    # Args definition
    args = parser.parse_args()
    
    graphAnalysis = args.graphAnalysis
    cliques = args.cliques
    doublons = args.doublons
    singletons = args.singletons
    
    caracteristicsInCSV(graphAnalysis, cliques, doublons, singletons)
