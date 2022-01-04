#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 15:02:38 2020

@author: niguilla
"""

import networkx as nx
import os

graph_repository = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/0_total_genes_exact_without_cn_prediction-hs-mm-clf-2019-03-13/graph_tr_set/"

gene_list_253_file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/all_knowledge_ensemble_ccds_2017/135/253_gene_list_hs.txt"
gene_hs_list = list()
with open(gene_list_253_file, 'r') as f_i:
    lecture = f_i.readlines()
    for line in lecture: 
        line = line.replace('\n', '')
        gene_hs_list.append(line)

list_tr_used = list()
for file in os.listdir(graph_repository):
    name = file.replace("obtention_Graph_tr_", '').replace("_for_Gephi.gexf", '').split('-')
    if name[0] in gene_hs_list:
#        print(name)
        g = nx.read_gexf(graph_repository+file)
        for tr in nx.Graph.nodes(g):
            if not tr in list_tr_used:
             list_tr_used.append(tr)
             
             
#print(list_tr_used)

hs_k = list()
hs_p = list()
mm_k = list()
mm_p = list()
clf_k = list()
clf_p = list()
for tr in list_tr_used:
    if tr.startswith("ENST"):
        hs_k.append(tr)
    elif tr.startswith("CGAT"):
        hs_p.append(tr)
    elif tr.startswith("ENSMUST"):
        mm_k.append(tr)
    elif tr.startswith("CGAMUST"):
        mm_p.append(tr)
    elif tr.startswith("ENSCAFT"):
        clf_k.append(tr)
    elif tr.startswith("CGACAFT"):
        clf_p.append(tr)

print("Number of tr known in hs:", len(hs_k))
print("Number of tr known in mm:", len(mm_k))
print("Number of tr known in clf:", len(clf_k))
print("Number of tr predicted in hs:", len(hs_p))
print("Number of tr predicted in mm:", len(mm_p)) 
print("Number of tr predicted in clf:", len(clf_p))