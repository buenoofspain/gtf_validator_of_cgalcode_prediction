#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 14:24:00 2019

@author: Nicolas Guillaudeux
@team: Dyliss
"""

import pandas

file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/mouse_junction_exons/junction_validation_mouse_2019-11-26/analysed_detail/nb_reads_mapped_tissue.csv"

tissue = list()
junction = list()
tissue_junc = dict()

                
df = pandas.read_table(file, sep ='\t',header = 0)
print(df.head())
print(df.columns)
for i in range(1, len(df.columns)):
    tissue.append(df.columns[i])
    nb_read = 0
    for junc in df[df.columns[i]]:
        junc = int(junc)
        if junc != 0:
            nb_read += 1
    tissue_junc[df.columns[i]] = nb_read
    junction.append(nb_read)

print(tissue_junc)
print(tissue)
print(junction)

print(tissue[30], junction[30])
print(tissue[25], junction[25])
print(tissue[34], junction[34])

with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/mouse_junction_exons/junction_validation_mouse_2019-11-26/analysed_detail/tissue_junction.csv", 'w') as f_o:
    f_o.write("Tissue"+'\t'+"Nb_reads"+'\n')
    for t in tissue_junc:
        f_o.write(str(t)+'\t'+str(tissue_junc[t])+'\n')
        
        
Brain_13': 32, 
'Brain_26': 32, 
'Brain_39': 33, 
'Colon_20': 21, 
'Colon_33': 17, 'Colon_7': 16, 'Duodenum_17': 13, 'Duodenum_30': 15, 'Duodenum_4': 9, 'Esophagus_12': 16, 'Esophagus_25': 21, 'Esophagus_38': 20, 'Heart_24': 15, 'Heart_37': 21, 'Ileum_19': 16, 'Ileum_32': 20, 'Ileum_6': 16, 'Jejunum_18': 21, 'Jejunum_31': 19, 'Jejunum_5': 17, 'Kidney_21': 14, 'Kidney_34': 17, 'Kidney_8': 15, 'Liver_15': 6, 'Liver_28': 0, 'Liver_2': 9, 'Muscle_22': 16, 'Muscle_35': 15, 'Muscle_9': 14, 'Pancreas_14': 2, 'Pancreas_1': 0, 'Pancreas_27': 0, 'Stomach_16': 14, 'Stomach_29': 15, 'Stomach_3': 18, 'Thymus_10': 19, 'Thymus_23': 19, 'Thymus_36': 17}