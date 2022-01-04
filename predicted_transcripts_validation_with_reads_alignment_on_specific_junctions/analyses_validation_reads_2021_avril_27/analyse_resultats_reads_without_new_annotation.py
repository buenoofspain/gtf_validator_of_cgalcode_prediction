#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 14:26:57 2021

@author: niguilla
"""
rep_result_reads = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_reads_alignment_on_specific_junctions/analyses_validation_annotation_2021_avril_23/"
file_read_hs = rep_result_reads+"nb_reads_mapped_tissue_hs.csv"
file_read_mm = rep_result_reads+"nb_reads_mapped_tissue_mm.csv"
file_read_clf = rep_result_reads+"nb_reads_mapped_tissue_clf.csv"

result_new_annot_hs = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/human/analyses_validation_annotation_2021_avril_23/result_validation_with_new_annotation.csv"
result_new_annot_mm = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/mouse/analyses_validation_annotation_2021_avril_23/result_validation_with_new_annotation.csv"
result_new_annot_clf = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/dog/analyses_validation_annotation_2021_avril_23/result_validation_with_new_annotation.csv"

tr_do_not_considered = ["CGACAFTENSMUST00000113589",
                        "CGACAFTENSMUST00000159840",
                        "CGACAFTENSMUST00000090821",
                        "CGACAFTENSMUST00000111287",
                        "CGACAFTENSMUST00000064371",
                        "CGACAFTENSMUST00000070782",
                        "CGACAFTENSMUST00000115205",
                        "CGACAFTENSMUST00000084738",
                        "CGACAFTENSMUST00000063091",
                        "CGACAFTENSMUST00000053540",
                        "CGACAFTENSMUST00000114691",
                        "CGACAFTENSMUST00000107046",
                        "CGACAFTENSMUST00000028377",
                        "CGACAFTENSMUST00000112140",
                        "CGACAFTENSMUST00000025023",
                        "CGACAFTENSMUST00000115454",
                        "CGACAFTENSMUST00000117269",
                        "CGACAFTENSMUST00000209532",
                        "CGACAFTENSMUST00000137769",
                        "CGACAFTENSMUST00000023854",
                        "CGACAFTENSMUST00000115323",
                        "CGACAFTENSMUST00000137750"]


####### HUMAN #######
list_tr_with_new_annot = list()
with open(result_new_annot_hs, 'r') as f_hs_annot_i:
    lecture = f_hs_annot_i.readlines()
    for line in lecture:
        line = line.replace('\n', '').replace("'", "").split(',')
        tr = line[0]
        list_tr_with_new_annot.append(tr)

print("On a validé", len(list_tr_with_new_annot), "avec de nouveaux annotations chez l'humain")

tr_with_read = list()
tr_without_read = list()
tr_total_with_junction_spe = list()

nb_spe = 0

with open(file_read_hs, 'r') as f_hs_read_i:
    lecture = f_hs_read_i.readlines()
    for line in lecture:
        
        if not line.startswith("Junctions"):
            line = line.replace('\n', '').split(',')
            junction = line[0]
            tr = junction.split("JUNC")[0]
            nb_reads = int(line[-1])
            if nb_reads != 0 and not tr in tr_with_read:
                tr_with_read.append(tr)
            elif nb_reads == 0 and not tr in tr_without_read:
                tr_without_read.append(tr)
#                print(tr, junction, nb_reads)
            if not tr in list_tr_with_new_annot and not tr in tr_total_with_junction_spe:
                tr_total_with_junction_spe.append(tr)
                nb_spe += 1
                
tr_with_read_ok = list()
for tr in tr_with_read:
#    print(tr)
#    print(tr in tr_with_read, tr in tr_without_read)
    if not tr in tr_without_read and not tr in tr_with_read_ok and not tr in list_tr_with_new_annot:
        tr_with_read_ok.append(tr)

#print("on retrouve", len(tr_with_read), "avec des reads alignés")
#print("on retrouve", len(tr_without_read), "sans reads alignés")
print("on retrouve", len(tr_total_with_junction_spe), "avec une/des jonctions spécifiques")
print("on valide", len(tr_with_read_ok), "avec les reads chez l'humain")





print()
####### MOUSE #######
list_tr_with_new_annot = list()
with open(result_new_annot_mm, 'r') as f_mm_annot_i:
    lecture = f_mm_annot_i.readlines()
    for line in lecture:
        line = line.replace('\n', '').replace("'", "").split(',')
        tr = line[0]
        list_tr_with_new_annot.append(tr)

print("On a validé", len(list_tr_with_new_annot), "avec de nouveaux annotations chez la souris")

tr_with_read = list()
tr_without_read = list()
tr_total_with_junction_spe = list()

nb_spe = 0

with open(file_read_mm, 'r') as f_mm_read_i:
    lecture = f_mm_read_i.readlines()
    for line in lecture:
        
        if not line.startswith("Junctions"):
            line = line.replace('\n', '').split(',')
            junction = line[0]
            tr = junction.split("JUNC")[0]
            nb_reads = int(line[-1])
            if nb_reads != 0 and not tr in tr_with_read:
                tr_with_read.append(tr)
            elif nb_reads == 0 and not tr in tr_without_read:
                tr_without_read.append(tr)
#                print(tr, junction, nb_reads)
            if not tr in list_tr_with_new_annot and not tr in tr_total_with_junction_spe:
                tr_total_with_junction_spe.append(tr)
                nb_spe += 1
                
tr_with_read_ok = list()
for tr in tr_with_read:
#    print(tr)
#    print(tr in tr_with_read, tr in tr_without_read)
    if not tr in tr_without_read and not tr in tr_with_read_ok and not tr in list_tr_with_new_annot:
        tr_with_read_ok.append(tr)

#print("on retrouve", len(tr_with_read), "avec des reads alignés")
#print("on retrouve", len(tr_without_read), "sans reads alignés")
print("on retrouve", len(tr_total_with_junction_spe), "avec une/des jonctions spécifiques")
print("on valide", len(tr_with_read_ok), "avec les reads chez la souris")




print()
####### DOG #######
list_tr_with_new_annot = list()
with open(result_new_annot_clf, 'r') as f_clf_annot_i:
    lecture = f_clf_annot_i.readlines()
    for line in lecture:
        line = line.replace('\n', '').replace("'", "").split(',')
        tr = line[0]
        if not tr in tr_do_not_considered:
            list_tr_with_new_annot.append(tr)

print("On a validé", len(list_tr_with_new_annot), "avec de nouveaux annotations chez le chien")

tr_with_read = list()
tr_without_read = list()
tr_total_with_junction_spe = list()

nb_spe = 0

with open(file_read_clf, 'r') as f_clf_read_i:
    lecture = f_clf_read_i.readlines()
    for line in lecture:
        if not line.startswith("Junctions"):
            line = line.replace('\n', '').split(',')
            junction = line[0]
            tr = junction.split("JUNC")[0]
            nb_reads = int(line[-1])
            if nb_reads != 0 and not tr in tr_with_read and not tr in tr_do_not_considered:
                tr_with_read.append(tr)
            elif nb_reads == 0 and not tr in tr_without_read and not tr in tr_do_not_considered:
                tr_without_read.append(tr)
#                print(tr, junction, nb_reads)
            if not tr in list_tr_with_new_annot and not tr in tr_total_with_junction_spe and not tr in tr_do_not_considered:
                tr_total_with_junction_spe.append(tr)
                nb_spe += 1
                
tr_with_read_ok = list()
for tr in tr_with_read:
#    print(tr)
#    print(tr in tr_with_read, tr in tr_without_read)
    if not tr in tr_without_read and not tr in tr_with_read_ok and not tr in list_tr_with_new_annot and not tr in tr_do_not_considered:
        tr_with_read_ok.append(tr)

#print("on retrouve", len(tr_with_read), "avec des reads alignés")
#print("on retrouve", len(tr_without_read), "sans reads alignés")
print("on retrouve", len(tr_total_with_junction_spe), "avec une/des jonctions spécifiques")
print("on valide", len(tr_with_read_ok), "avec les reads chez le chien")