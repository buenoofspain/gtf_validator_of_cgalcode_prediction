#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 16:44:48 2021

@author: niguilla
"""






output_cga_run2 = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/0_total_genes_exact_without_cn_prediction-hs-mm-clf-2019-03-13/output_run2.txt"

with open(output_cga_run2, 'r') as output_cga_i:
    test = 0
    count_known_known = 0
    count_known_predicted = 0
    count_predicted_predicted = 0
    lecture = output_cga_i.readlines()
    for line in lecture:
        if line.startswith("CO "):
            line_co = line.split()#[3:5]
            if not line_co[4] == "No-executable" and not line_co[4] == "No-orthologous-transcript":
                test += 1
                if line_co[3].startswith("CGA"):
                    if line_co[4].startswith("CGA"):
                        count_predicted_predicted += 1
                    else:
                        count_known_predicted += 1
                else:
                    if line_co[4].startswith("CGA"):
                        count_known_predicted += 1
                    else:
                        count_known_known += 1
                
    print("Number of known-known transcript relationship:", count_known_known)
    print("Number of known-predicted transcript relationship:", count_known_predicted)
    print("Number of predicted-predicted transcript relationship:", count_predicted_predicted)
    print("Total number of transcript relationship:", count_known_known+count_known_predicted+count_predicted_predicted)
    print(test)