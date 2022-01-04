#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 15:13:02 2021

@author: niguilla
"""

list_tr_with_junction = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/list_tr_tot.txt"

dict_human = dict()
dict_mouse = dict()
dict_dog = dict()

with open(list_tr_with_junction, 'r') as f_i:
    lecture = f_i.readlines()
    for line in lecture:
        line = line.replace('\n', '')
        print(line)
        if line.startswith("CGAT"):
            if line in dict_human:
                dict_human[line] += 1
            else:
                dict_human[line] = 1
        elif line.startswith("CGAMUST"):
            if line in dict_mouse:
                dict_mouse[line] += 1
            else:
                dict_mouse[line] = 1
                
        else:
            if line in dict_dog:
                dict_dog[line] += 1
            else:
                dict_dog[line] = 1
                

print(dict_human)
print(dict_mouse)
print(dict_dog)