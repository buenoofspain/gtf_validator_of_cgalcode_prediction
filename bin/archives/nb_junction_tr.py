#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:24:31 2019

@author: niguilla
"""

tr_list = list()
tr_with_specific_junction = list()

with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/list_201_tr_clf.txt", 'r') as t:
    lecture = t.readlines()
    for tr in lecture:
        tr = tr.replace('\n', '')
        tr_list.append(tr)
        
print(len(tr_list))

with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/junction_validation_dog_201_2019-10-01/specific_junction.bed", 'r') as j_f:
    lecture = j_f.readlines()
    for info_j in lecture:
        info_j = info_j.replace('\n', '').split('\t')
        junc_name = info_j[3]
        junc_name = junc_name.split("JUNC")
        if junc_name[0] in tr_list and junc_name[0] not in tr_with_specific_junction:
            tr_with_specific_junction.append(junc_name[0])
        
print(tr_with_specific_junction)
print(len(tr_with_specific_junction))