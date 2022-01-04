#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:05:14 2019

@author: niguilla
"""
tr_validate_list = list()
tr_non_validate_list = list()

with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/list_111_transcript_validated.txt", 'r') as validated_t:
    lecture = validated_t.readlines()
    for tr in lecture:
        tr = tr.replace('\n', '')
        tr_validate_list.append(tr)
   
print(len(tr_validate_list))

with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/list_312_tr_clf.txt", 'r') as non_validated_t:
    lecture = non_validated_t.readlines()
    for tr in lecture:
        tr = tr.replace('\n', '')
        if not tr in tr_validate_list:
            tr_non_validate_list.append(tr)
        
print(len(tr_non_validate_list))

with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/list_201_tr_clf.txt", 'w') as non_validate_t_o:
    for tr in tr_non_validate_list:
        non_validate_t_o.write(tr+'\n')
