#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 15:26:17 2021

@author: niguilla
"""

dict_result_annotation = dict()
nb_tr_validated_by_ens96 = 0
nb_tr_validated_by_ens98 = 0
nb_tr_validated_by_ens102 = 0
nb_tr_validated_by_ens103 = 0
nb_tr_validated_by_ucsc = 0
nb_tr_validated_by_feelnc = 0

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

with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/dog/analyses_validation_annotation_2021_avril_23/result_via_ENS96_clf.txt", 'r') as ens96_f_i:
    lecture = ens96_f_i.readlines()
    db_source = "ENS96"
    for tr in lecture:
#        nb_tr_validated_by_ens96 += 1
        tr = tr.replace('\n', '')
        if not tr in tr_do_not_considered:
            if not tr in dict_result_annotation:
                nb_tr_validated_by_ens96 += 1
                dict_result_annotation[tr] = [db_source]
            elif not db_source in dict_result_annotation[tr]:
                dict_result_annotation[tr].append(db_source)
            else:
                print("ERROR HERE")
        else:
            print(tr, "not considered in", db_source)

with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/dog/analyses_validation_annotation_2021_avril_23/result_via_ENS98_clf.txt", 'r') as ens98_f_i:
    lecture = ens98_f_i.readlines()
    db_source = "ENS98"
    for tr in lecture:
#        nb_tr_validated_by_ens98 += 1
        tr = tr.replace('\n', '')
        if not tr in tr_do_not_considered:
            if not tr in dict_result_annotation:
                nb_tr_validated_by_ens98 += 1
                dict_result_annotation[tr] = [db_source]
            elif not db_source in dict_result_annotation[tr]:
                dict_result_annotation[tr].append(db_source)
            else:
                print("ERROR HERE")
        else:
            print(tr, "not considered in", db_source)

with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/dog/analyses_validation_annotation_2021_avril_23/result_via_ENS102_clf.txt", 'r') as ens102_f_i:
    lecture = ens102_f_i.readlines()
    db_source = "ENS102"
    for tr in lecture:
#        nb_tr_validated_by_ens103 += 1
        tr = tr.replace('\n', '')
        if not tr in tr_do_not_considered:
            if not tr in dict_result_annotation:
                nb_tr_validated_by_ens102 += 1
                dict_result_annotation[tr] = [db_source]
            elif not db_source in dict_result_annotation[tr]:
                dict_result_annotation[tr].append(db_source)
            else:
                print("ERROR HERE")
        else:
            print(tr, "not considered in", db_source)

with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/dog/analyses_validation_annotation_2021_avril_23/result_via_ENS103_clf.txt", 'r') as ens103_f_i:
    lecture = ens103_f_i.readlines()
    db_source = "ENS103"
    for tr in lecture:
#        nb_tr_validated_by_ens103 += 1
        tr = tr.replace('\n', '')
        if not tr in tr_do_not_considered:
            if not tr in dict_result_annotation:
                nb_tr_validated_by_ens103 += 1
                dict_result_annotation[tr] = [db_source]
            elif not db_source in dict_result_annotation[tr]:
                dict_result_annotation[tr].append(db_source)
            else:
                print("ERROR HERE")
        else:
            print(tr, "not considered in", db_source)
     
with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/dog/analyses_validation_annotation_2021_avril_23/result_via_UCSC_clf.txt", 'r') as ucsc_f_i:
    lecture = ucsc_f_i.readlines()
    db_source = "UCSC"
    for tr in lecture:
#        nb_tr_validated_by_ucsc += 1
        tr = tr.replace('\n', '')
        if not tr in tr_do_not_considered:
            if not tr in dict_result_annotation:
                nb_tr_validated_by_ucsc += 1
                dict_result_annotation[tr] = [db_source]
            elif not db_source in dict_result_annotation[tr]:
                dict_result_annotation[tr].append(db_source)
            else:
                print("ERROR HERE")
        else:
            print(tr, "not considered in", db_source)
    
with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/dog/analyses_validation_annotation_2021_avril_23/result_via_FEELnc_clf.txt", 'r') as feelnc_f_i:
    lecture = feelnc_f_i.readlines()
    db_source = "FEELNC"
    for tr in lecture:
#        nb_tr_validated_by_feelnc += 1
        tr = tr.replace('\n', '')
        if not tr in tr_do_not_considered:
            if not tr in dict_result_annotation:
                nb_tr_validated_by_feelnc += 1
                dict_result_annotation[tr] = [db_source]
            elif not db_source in dict_result_annotation[tr]:
                dict_result_annotation[tr].append(db_source)
            else:
                print("ERROR HERE")
        else:
            print(tr, "not considered in", db_source)


print("Number of transcript found in databases:", len(dict_result_annotation))

print("Number of transcript found in ENS96:", nb_tr_validated_by_ens96)
print("Number of transcript found in ENS98:", nb_tr_validated_by_ens98)
print("Number of transcript found in ENS102:", nb_tr_validated_by_ens102)
print("Number of transcript found in ENS103:", nb_tr_validated_by_ens103)
print("Number of transcript found in UCSC:", nb_tr_validated_by_ucsc)
print("Number of transcript found in FEELNC:", nb_tr_validated_by_feelnc)


with open("result_validation_with_new_annotation.csv", 'w') as result_f_o:
    for tr in dict_result_annotation:
        source_val = dict_result_annotation[tr][0]
        result_f_o.write(tr+','+source_val+'\n')
        
for tr in dict_result_annotation:
    if "ENS102" in dict_result_annotation[tr][0] or "ENS103" in dict_result_annotation[tr][0]:
        print(tr)