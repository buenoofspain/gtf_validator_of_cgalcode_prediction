#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 14:38:47 2020

@author: niguilla

This script extract all information from a set of graphs of interest.
"""

import os
import argparse
import networkx
from matplotlib_venn import venn3, venn2, venn2_circles, venn3_circles
import matplotlib.pyplot as plt
import numpy as np

#def graphTranscriptAnalysis(list_conserved_graph_delimitor, graph_dir, name_sp1, name_sp2, name_sp3):
#    """
#    Function description:
#        Analyse phylogenetic relationship between transcripts in graph.
#        
#    Parameters:
#    - list_conserved_graph_delimitor: list of files where graph delimitor are conserved to analysis
#    - Graph: the transcripts graph containing nodes and edges
#    - name_sp1/2/3: species name 1/2/3 ('HS', 'MUS', 'CAF')
#    
#    Output:
#    - graph_tr_analysis.csv: an output table with the detail of graph analysis
#    """
#    #graph_dir = os.path(graph_directory)
#    graph_list_files = os.listdir(graph_dir)
#    nb_file = 0
#    graph_not_conserved = 0
#    triplet_number = 0
#    graph_conserved_list = list()
#    
#    with open("graph_tr_analysis.csv", "w") as output:
#        output.write("triplet_number,Genes_triplet,clique_number,human_abs_number,mouse_abs_number,dog_abs_number,human_specific_number,mouse_specific_number,dog_specific_number,total\n")
#    
#        for file in graph_list_files:
#            tr_file_name_test = file[:15]+file[18:]
#            if tr_file_name_test in list_conserved_graph_delimitor:
#    
#                graph_conserved_list.append(file)
#                nb_file +=1
#                path = graph_dir+"/"+file
#                Graph = nx.read_gexf(path, node_type=None, relabel=False, version='1.1draft')
#                
#                #To analyse transcripts graph
#            
#                problematic_nodes = 0
#                clique_number, nodes_spared, specific_nodes_number = 0, 0, 0
#
#                human_specific_number, mouse_specific_number, dog_specific_number = 0, 0, 0
#                human_abs_number, mouse_abs_number, dog_abs_number = 0, 0, 0
#
#                hs_gene, mm_gene, clf_gene = "", "", ""
#
#                clique_list = list()
#
#                list_nodes_in_graph = list(Graph.nodes())
#                list_nodes = list_nodes_in_graph[:]
#                
#                alternative_transcription = False
#                
#                #nodes_number = Graph.number_of_nodes() #Number of nodes in graph
#                #edges_number = Graph.number_of_edges() #Number if edges in graph
#
#                for node in list_nodes_in_graph:
#                    if Graph.nodes[node]['gene'].startswith('ENSG') and hs_gene == '':
#                        hs_gene = Graph.nodes[node]['gene']
#                    if Graph.nodes[node]['gene'].startswith('ENSMUSG'):
#                        mm_gene = Graph.nodes[node]['gene']
#                    if Graph.nodes[node]['gene'].startswith('ENSCAFG'):
#                        clf_gene = Graph.nodes[node]['gene']
#
#                    if node in list_nodes:
#                        clique = 0
#                        neighbors_nodes = list(Graph.neighbors(node))
#
#                        #To treat specific node in graph
#                        if checkGraphSpecificNode(Graph, node):
#                            if Graph.nodes[node]['species'] == name_sp1:
#                                print "Specific node: "+str(node)+" - "+str(name_sp1)+": AGE PRIMATE OR NEWER"
#                                human_specific_number += 1
#                            elif Graph.nodes[node]['species'] == name_sp2:
#                                print "Specific node: "+str(node)+" - "+str(name_sp2)+": AGE RODENT OR NEWER"
#                                mouse_specific_number += 1
#                            elif Graph.nodes[node]['species'] == name_sp3:
#                                print "Specific node: "+str(node)+" - "+str(name_sp3)+": AGE LAURASIATHERIA OR AGE EUTHERIA OR OLDER"
#                                dog_specific_number += 1
#                            specific_nodes_number += 1
#                            list_nodes.remove(node)
#
#                        #To treat a clique delimitor in graph
#                        elif checkGraphDegreeClique(Graph, node):
#                            for neighbor in neighbors_nodes:
#                                neighbor_neighbors_nodes = list(Graph.neighbors(neighbor))
#                                if len(neighbor_neighbors_nodes) == 2 and (neighbor_neighbors_nodes[0] in neighbors_nodes or neighbor_neighbors_nodes[1] in neighbors_nodes) and node in neighbor_neighbors_nodes:
#                                    if checkGraphDegreeClique(Graph, neighbor) == True:
#                                        # CAS CLIQUE
#                                        clique+=1
#                                    else:
#                                        print "ERROR 1 HERE:"+str(neighbor)+" "+str(node)+str(neighbor_neighbors_nodes)
#                                        problematic_nodes += 1
#                                        list_nodes.remove(node)
##                                        break
#                                else:
#                                    print "ERROR 2 HERE:"+str(neighbor)+" "+str(node)+str(neighbor_neighbors_nodes)
#                                    problematic_nodes += 1
#                                    alternative_transcription = True
#                                    try:
#                                        list_nodes.remove(node)
#                                    except ValueError:
#                                        pass
##                                    break
#                
#                            if clique == 2:
#                                nodes_in_clique = str(node)
#                                clique_list = neighbors_nodes[:]
#                                clique_list.append(node)
#                                clique_number += 1
#                                for neighbor in neighbors_nodes:
#                                    nodes_in_clique += " "+str(neighbor)
#                                    list_nodes.remove(neighbor)
#                                list_nodes.remove(node)
#                                print "Clique: "+str(nodes_in_clique)+" - all species: AGE EUTHERIA OR OLDER"
#
#                        #To treat nodes spared by two species (absent in one species)
#                        elif checkGraphNodeSpared(Graph, node):
#                            for neighbor in neighbors_nodes:
#                                neighbor_neighbors_nodes = list(Graph.neighbors(neighbor))
#                                if len(neighbor_neighbors_nodes) == 1 and node in neighbor_neighbors_nodes:
#                                    if (Graph.nodes[node]['species'] == name_sp2 or Graph.nodes[node]['species'] == name_sp3) and (Graph.nodes[neighbor]['species'] == name_sp2 or Graph.nodes[neighbor]['species'] == name_sp3):
#                                    #(node in positions_sp2 or node in positions_sp3) and (neighbor in positions_sp2 or neighbor in positions_sp3):
#                                        #" ABS chez HS ")
#                                        print "Shared nodes: "+str(node)+" "+str(neighbor)+" - "+str(name_sp2)+" "+str(name_sp3)+" (missing in "+str(name_sp1)+"): AGE EUTHERIA OR OLDER"
#                                        list_nodes.remove(neighbor)
#                                        list_nodes.remove(node)
#                                        nodes_spared += 1
#                                        human_abs_number += 1
#                                        
#                                    elif (Graph.nodes[node]['species'] == name_sp1 or Graph.nodes[node]['species'] == name_sp3) and (Graph.nodes[neighbor]['species'] == name_sp1 or Graph.nodes[neighbor]['species'] == name_sp3):
#                                    #(node in positions_sp1 or node in positions_sp3) and (neighbor in positions_sp1 or neighbor in positions_sp3):
#                                        #" ABS chez MM "
#                                        print "Shared nodes: "+str(node)+" "+str(neighbor)+" - "+str(name_sp1)+" "+str(name_sp3)+" (missing in "+str(name_sp2)+"): AGE EUTHERIA OR OLDER"
#                                        list_nodes.remove(neighbor)
#                                        list_nodes.remove(node)
#                                        nodes_spared += 1
#                                        mouse_abs_number += 1
#                                        
#                                    elif (Graph.nodes[node]['species'] == name_sp1 or Graph.nodes[node]['species'] == name_sp2) and (Graph.nodes[neighbor]['species'] == name_sp1 or Graph.nodes[neighbor]['species'] == name_sp2):
#                                    #(node in positions_sp1 or node in positions_sp2) and (neighbor in positions_sp1 or neighbor in positions_sp2):
#                                        #" ABS chez CLF ")
#                                        print "Shared nodes: "+str(node)+" "+str(neighbor)+" - "+str(name_sp1)+" "+str(name_sp2)+" (missing in "+str(name_sp3)+"): AGE EUARCHONTOGLIRES OR AGE EUTHERIA OR OLDER"
#                                        list_nodes.remove(neighbor)
#                                        list_nodes.remove(node)
#                                        nodes_spared += 1
#                                        dog_abs_number += 1
#                                                    
#                                else:
#                                    print "ERROR 3 HERE:"+str(neighbor)+" "+str(node)+str(neighbor_neighbors_nodes)
#                                    problematic_nodes += 1
#                                    list_nodes.remove(node)
#                                    alternative_transcription = True
##                                    break
#                        else:
#                            print "ERROR 4 HERE:"+str(node)
#                            problematic_nodes += 1
#                            list_nodes.remove(node)
##                            break
#
#                print ">>> Genes: "+str(hs_gene)+" "+str(mm_gene)+" "+str(clf_gene)
#                print ">>> Clique number: "+str(clique_number)
#                print ">>> Specific nodes number: "+str(specific_nodes_number)
#                print ">>> Nodes pairwise shared by two species: "+str(nodes_spared)
#                print ">>> Problematic nodes number: "+str(problematic_nodes)
#                print ">>> Number of nodes not treated: "+str(len(list_nodes))+" / "+str(len(Graph.nodes()))+" "+str(list_nodes)
#                if alternative_transcription == True:
#                    print ">>> Alternative transcription found"
#                
#                if problematic_nodes == 0 and hs_gene != "" and mm_gene != "" and clf_gene != "":
#                    output.write(str(triplet_number)+","+
#                                 str(hs_gene)+"-"+str(mm_gene)+"-"+str(clf_gene)+","+
#                                 str(clique_number)+","+
#                                 str(human_abs_number)+","+str(mouse_abs_number)+","+str(dog_abs_number)+","+
#                                 str(human_specific_number)+","+str(mouse_specific_number)+","+str(dog_specific_number)+","+
#                                 str(len(Graph.nodes()))+"\n")
#                    triplet_number += 1
#                    print "################################################################################"
#                    
#                else:
#                    print ">>> !!!Graph tr not conserved!!!"
#                    print "################################################################################"
#                    graph_not_conserved += 1
#                    try:
#                        graph_conserved_list.remove(file)
#                    except ValueError:
#                        pass
#            
#    print ">>> Number of graph conserved: "+str(nb_file - graph_not_conserved)+"/"+str(nb_file)
#    print ">>> Number of files conserved: "+str(len(graph_conserved_list))
#    print "################################################################################"
#    statisticalAnalysis('graph_tr_analysis.csv', nb_file)
#    
#    return graph_conserved_list

def main_run(rep_graph, file_list_genes):
    """
    Function description:
        ...
        
    Input:
        - 
        
    Output:
        -
        
    """
    list_genes = list()
    
    with open(file_list_genes, 'r') as f_i:
        lecture = f_i.readlines()
        for gene_triplet in lecture:
            gene_triplet = gene_triplet.replace('\n', '')
            print(gene_triplet)

    

#def parse():
#    """
#    Method to create a parser for command-line options
#    """
#    parser = argparse.ArgumentParser(description="A SCRIPT TO OBTAIN SPECIFIC EXON JUNCTION NOT KNOWN IN ENSEMBL.")
#    # Create command-line parser for all options and arguments to give
#    parser.add_argument("-g", "--graph",
#                        dest = "rep_graph", 
#                        metavar = "REPOSITORY OF TRANSCRIPT GRAPHS",
#                        default = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/0_total_genes_exact_without_cn_prediction-hs-mm-clf-2019-03-13/graph_tr_set/",
#                        help = "Give the path of the repository where are stored the transcript graphs.",
#                        required = True)
#    
#    parser.add_argument("-l", "--list",
#                        dest = "list_genes", 
#                        metavar = "LIST OF INTERESTED GENES", 
#                        default = "/home/niguilla/Documents/these_nguillaudeux/analyse_graph_signaux_similaires_tr_differents/list_120.txt",
#                        help = "Give the list of interested genes to analyse transcript graphs.",
#                        required = True)
#    
#    
#    return parser.parse_args()  
    
    
    
if __name__ == "__main__":
#    OPTIONS = parse()
    
#    graph_repository = OPTIONS.rep_graph
    graph_repository = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/0_total_genes_exact_without_cn_prediction-hs-mm-clf-2019-03-13/graph_tr_set/"
#    interested_gene_list = OPTIONS.list_genes
    interested_gene_list = "/home/niguilla/Documents/these_nguillaudeux/analyse_graph_signaux_similaires_tr_differents/list_120.txt"
    
    main_run(graph_repository, interested_gene_list)