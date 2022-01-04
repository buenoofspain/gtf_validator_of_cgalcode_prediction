#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:44:25 2019

@author: niguilla
"""

def split_gtf_by_chr():
    '''
    Function description:
        This function split a gtf file in different gtf file split by chr number.
    '''
    chr_list = list()
    #with open('/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/dog/igdr/split_gtf/chr.sorted.txt', 'r') as f_i:
    with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/dog/igdr/split_gtf/junction/chr_cga_pred.sorted.txt", 'r') as f_i:
        lecture = f_i.readlines()
        for line in lecture:
            line = line.split('\n')
            if not line[0] in chr_list:
                chr_list.append(line[0])
                
    for chr_nb in chr_list:
        #with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/dog/igdr/split_gtf/canfam3_cons_annot_TritouPuppetMasterChief_23-03-2016_lncClasse_geneBiotype_withEnsCds_withSilicoCds_withGoodGeneSource.gtf", 'r') as i_f:
        #with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/dog/igdr/split_gtf/junction/specific_junction_cga_ens2017.sort.bed", 'r') as i_f:
        #with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/dog/igdr/split_gtf/junction/specific_junction_cga_chr.sort.bed", 'r') as i_f:
        with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/dog/igdr/split_gtf/dog_cgalcode_tot.gtf", 'r') as i_f:
            with open("/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/dog/igdr/split_gtf/gtf_cga/split_chr"+chr_nb+".bed", 'a') as o_f:
                lecture = i_f.readlines()
                for line in lecture:
                    line_split = line.replace('\n', '').split('\t')
                    if line_split[0] == chr_nb:# and line[1] == "\t":
                        o_f.write(line)
        
def test_specific_junctions():
    '''
    Function description:
            this function check if specific junction are present.
    '''
    chr_nb = str(1)
    pred_file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/dog/igdr/split_gtf/cga_pred_split/split_cga_chr"+chr_nb+".gtf"
    ref_file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/dog/igdr/split_gtf/junction_ref_split/split_chr"+chr_nb+".bed"
    
    print(ref_file)
    
    pred_list = dict()
    ref_list = dict()

    with open(pred_file, 'r') as pred_f:
        lecture_pred = pred_f.readlines()
        
        with open(ref_file, 'r') as ref_f:
            lecture_ref = ref_f.readlines()
            
            for line_pred in lecture_pred:
                line_pred = line_pred.replace('\n', '').split('\t')
                tr = line_pred[3][:-9]
                if not tr in pred_list:
                    pred_list[tr] = []
                if (line_pred[0], line_pred[1], line_pred[2], line_pred[4], line_pred[5], line_pred[6], line_pred[7], line_pred[8], line_pred[9], line_pred[10], line_pred[11]) not in pred_list[tr]:
                    pred_list[tr].append((line_pred[0], line_pred[1], line_pred[2], line_pred[4], line_pred[5], line_pred[6], line_pred[7], line_pred[8], line_pred[9], line_pred[10], line_pred[11]))
                
            for line_ref in lecture_ref: 
                line_ref = line_ref.replace('\n', '').split('\t')
                tr = line_ref[3][:-9]
                if not tr in ref_list:
                    ref_list[tr] = []
                if (line_ref[0], line_ref[1], line_ref[2], line_ref[4], line_ref[5], line_ref[6], line_ref[7], line_ref[8], line_ref[9], line_ref[10], line_ref[11]) not in ref_list[tr]:
                    ref_list[tr].append((line_ref[0], line_ref[1], line_ref[2], line_ref[4], line_ref[5], line_ref[6], line_ref[7], line_ref[8], line_ref[9], line_ref[10], line_ref[11]))
            
    
    nb_test = 0
    for tr in pred_list:
        for junc_pred in pred_list[tr]:
            if junc_pred not in ref_list[tr]:
                print(junc_pred)
                nb_test += 1
#            else:
#                nb_test += 1
##                print(junc_pred)
    print(nb_test)
        
if __name__ == "__main__":
    
    split_gtf_by_chr()
#    test_specific_junctions()
    pass
    