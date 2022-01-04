#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 15:14:40 2020

@author: niguilla
"""
def number_of_specific_junctions(alignement_specific_junctions_file):
    '''
    '''
    transcript_with_junction_aligned = dict()
    transcript_without_junction_aligned = dict()
    
    with open(alignement_specific_junctions_file, 'r') as alignment_file_i:
        lecture = alignment_file_i.readlines()
        for line in lecture:
            line = line.replace('\n', '')
            if not line.startswith("Junctions"):
                junction = line.split('\t')
                transcript = junction[0].split('JUNC')[0]
                nb_read_aligned = int(junction[-1])
                
                if nb_read_aligned == 0:
                    if not transcript in transcript_without_junction_aligned:
                        transcript_without_junction_aligned[transcript] = 1
                    else:
                        nb_junction = transcript_without_junction_aligned[transcript] +1
                        transcript_without_junction_aligned[transcript] = nb_junction
                    if not transcript in transcript_with_junction_aligned:
                        transcript_with_junction_aligned[transcript] = 0
    
                else:
                    if not transcript in transcript_with_junction_aligned:
                        transcript_with_junction_aligned[transcript] = 1
                    else:
                        nb_junction = transcript_with_junction_aligned[transcript] +1
                        transcript_with_junction_aligned[transcript] = nb_junction
                    if not transcript in transcript_without_junction_aligned:
                        transcript_without_junction_aligned[transcript] = 0
    
    transcript_junction_read_state = dict()
    for transcript in transcript_with_junction_aligned:
        if transcript in transcript_with_junction_aligned and transcript in transcript_without_junction_aligned:
            nb_specific_junction_with_read = transcript_with_junction_aligned[transcript]
            nb_specific_junction_without_read = transcript_without_junction_aligned[transcript]
            nb_specific_junction = nb_specific_junction_with_read + nb_specific_junction_without_read
            transcript_junction_read_state[transcript] = (nb_specific_junction_with_read, nb_specific_junction)
        else:
            print("ERROR")
            
    return transcript_junction_read_state
    
def annotation(annotation_result_file):
    '''
    '''
    dict_annot = dict()
    
    with open(annotation_result_file, 'r') as annot_i:
        lecture = annot_i.readlines()
        for line in lecture:
            line = line.replace('\n', '').split('\t')
            transcript = line[0]
            annot = line[1]
            dict_annot[transcript] = annot
    
    return dict_annot

def sourceId(database):
    if database == "ens96":
        return "1"
    elif database == "ens98":
        return "2"
    elif database == "ucsc":
        return "3"
    elif database == "xbseq":
        return "4"
    elif database == "feelnc":
        return "5"
    else:
        return "NA"

if __name__ == "__main__":
    
    alignment_file_hs = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/add_validation_in_db/nb_reads_mapped_tissue_hs_sum.csv"
    alignment_file_mm = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/add_validation_in_db/nb_reads_mapped_tissue_mm_sum.csv"
    alignment_file_clf = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/add_validation_in_db/nb_reads_mapped_tissue_clf_sum.csv"
    
    annotation_file_hs = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/add_validation_in_db/evidence_hs_annot"
    annotation_file_mm = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/add_validation_in_db/evidence_mm_annot"
    annotation_file_clf = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/add_validation_in_db/evidence_clf_annot"
    
    transcript_details_db = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/add_validation_in_db/idTranscript_sourceTId_state"

    transcript_specific_junction_read_state_hs = number_of_specific_junctions(alignment_file_hs)
    transcript_specific_junction_read_state_mm = number_of_specific_junctions(alignment_file_mm)
    transcript_specific_junction_read_state_clf = number_of_specific_junctions(alignment_file_clf)
    
    annot_hs = annotation(annotation_file_hs)
    annot_mm = annotation(annotation_file_mm)
    annot_clf = annotation(annotation_file_clf)
    
    dict_transcript_info = dict()
    with open(transcript_details_db, 'r') as info_tr_i:
        lecture = info_tr_i.readlines()
        for line in lecture:
            line = line.replace('\n', '').split('\t')
            idTranscript = line[0]
            transcript = line[1]
            state = line[2]
            dict_transcript_info[transcript] = idTranscript
            
    count_known = 0
    with open("resultValidationTcsv", 'w') as validation_o:
        for transcript in dict_transcript_info:
            if transcript.startswith("CGAT"):
                if transcript in annot_hs: 
                    sourceTId = sourceId(annot_hs[transcript])
                    validation_o.write(dict_transcript_info[transcript] +' '+ "1" +' '+ sourceTId +' '+ "NULL" +' '+ "NULL" +'\n')
                    
                elif transcript in transcript_specific_junction_read_state_hs:
                    if transcript_specific_junction_read_state_hs[transcript][0] == transcript_specific_junction_read_state_hs[transcript][1]:
                        validation_o.write(dict_transcript_info[transcript] +' '+ "3" +' '+ "6" +' '+ str(transcript_specific_junction_read_state_hs[transcript][1]) +' '+ str(transcript_specific_junction_read_state_hs[transcript][0]) +'\n')
                    else:
                        validation_o.write(dict_transcript_info[transcript] +' '+ "0" +' '+ "6" +' '+ str(transcript_specific_junction_read_state_hs[transcript][1]) +' '+ str(transcript_specific_junction_read_state_hs[transcript][0]) +'\n')
                        
                else:
                    validation_o.write(dict_transcript_info[transcript] +' '+ "2" +' '+ "0" +' '+ "NULL" +' '+ "NULL" +'\n')
                    
            elif transcript.startswith("CGAMUST"):
                if transcript in annot_mm: 
                    sourceTId = sourceId(annot_mm[transcript])
                    validation_o.write(dict_transcript_info[transcript] +' '+ "1" +' '+ sourceTId +' '+ "NULL" +' '+ "NULL" +'\n')
                    
                elif transcript in transcript_specific_junction_read_state_mm:
                    if transcript_specific_junction_read_state_mm[transcript][0] == transcript_specific_junction_read_state_mm[transcript][1]:
                        validation_o.write(dict_transcript_info[transcript] +' '+ "3" +' '+ "7" +' '+ str(transcript_specific_junction_read_state_mm[transcript][1]) +' '+ str(transcript_specific_junction_read_state_mm[transcript][0]) +'\n')
                    else:
                        validation_o.write(dict_transcript_info[transcript] +' '+ "0" +' '+ "7" +' '+ str(transcript_specific_junction_read_state_mm[transcript][1]) +' '+ str(transcript_specific_junction_read_state_mm[transcript][0]) +'\n')
                        
                else:
                    validation_o.write(dict_transcript_info[transcript] +' '+ "2" +' '+ "0" +' '+ "NULL" +' '+ "NULL" +'\n')
                    
            elif transcript.startswith("CGACAFT"):
                if transcript in annot_clf: 
                    sourceTId = sourceId(annot_clf[transcript])
                    validation_o.write(dict_transcript_info[transcript] +' '+ "1" +' '+ sourceTId +' '+ "NULL" +' '+ "NULL" +'\n')
                    
                elif transcript in transcript_specific_junction_read_state_clf:
                    if transcript_specific_junction_read_state_clf[transcript][0] == transcript_specific_junction_read_state_clf[transcript][1]:
                        validation_o.write(dict_transcript_info[transcript] +' '+ "3" +' '+ "8" +' '+ str(transcript_specific_junction_read_state_clf[transcript][1]) +' '+ str(transcript_specific_junction_read_state_clf[transcript][0]) +'\n')
                    else:
                        validation_o.write(dict_transcript_info[transcript] +' '+ "0" +' '+ "8" +' '+ str(transcript_specific_junction_read_state_clf[transcript][1]) +' '+ str(transcript_specific_junction_read_state_clf[transcript][0]) +'\n')
                        
                else:
                    validation_o.write(dict_transcript_info[transcript] +' '+ "2" +' '+ "0" +' '+ "NULL" +' '+ "NULL" +'\n')
                    
            else:
                if transcript.startswith("ENS"):
                    count_known += 1
    print(count_known)
#    print(dict_transcript_info["CGATENSMUST00000152555"])
            
#    print(transcript_specific_junction_read_state_hs)
#    for tr in dict_transcript_info:
#        if tr.startswith("CGAT"):
#            print(tr, dict_transcript_info[tr], transcript_specific_junction_read_state_hs[tr])