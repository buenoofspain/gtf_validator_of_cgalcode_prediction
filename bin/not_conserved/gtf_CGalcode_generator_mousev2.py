#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 11:10:26 2019

@author: Nicolas Guillaudeux

This script create a gtf file adapted to the dog information and analyse 
results.
"""
import os
import re
import argparse
import subprocess

def gene_tr_relation(gene_tr):
    '''
    Function description:
        This function use the file with the relation between genes and their
        transcripts (Ensembl and predicted) to return a dictionary with this 
        relation.
        
    Input:
        gene_tr: a file with the genes/transcripts relation.
    
    Return:
        gene_tr_dict: a dictionary within genes/transcripts relation.        
    '''
    gene_tr_dict = dict()
    with open(gene_tr, 'r') as gene_tr_file:
        lecture = gene_tr_file.readlines()
        for line in lecture:
            if not 'gene_ID' in line:
                line = line.replace('\n', '').split(',')
                if line[0] not in gene_tr_dict:
                    gene_tr_dict[line[0]] = []
                if os.path.exists(rep+line[1]+'.fasta'):
                    gene_tr_dict[line[0]].append(line[1])
                    
    return gene_tr_dict

def gene_tr_ensembl_relation(gene_tr):
    '''
    Function description:
        This function use the file with the relation between genes and their
        Ensembl transcripts to return a dictionary with this relation.
        
    Input:
        gene_tr: a file with the genes/transcripts Ensembl relation.
    
    Return:
        gene_tr_dict: a dictionary within genes/transcripts Ensembl relation.        
    '''
    gene_tr_dict = dict()
    with open(gene_tr, 'r') as gene_tr_file:
        lecture = gene_tr_file.readlines()
        for line in lecture:
            if not 'gene_ID' in line:
                line = line.replace('\n', '').split(',')
                if line[0] not in gene_tr_dict:
                    gene_tr_dict[line[0]] = []
                if os.path.exists(rep+line[1]+'.fasta'):
                    if line[1].startswith("ENS"):
                        gene_tr_dict[line[0]].append(line[1])
                    
    return gene_tr_dict

def gene_tr_pred_relation(gene_tr):
    '''
    Function description:
        This function use the file with the relation between genes and their
        predicted transcripts to return a dictionary with this relation.
        
    Input:
        gene_tr: a file with the genes/transcripts predicted relation.
    
    Return:
        gene_tr_dict: a dictionary within genes/transcripts predicted relation.        
    '''
    gene_tr_dict = dict()
    with open(gene_tr, 'r') as gene_tr_file:
        lecture = gene_tr_file.readlines()
        for line in lecture:
            if not 'gene_ID' in line:
                line = line.replace('\n', '').split(',')
                if line[0] not in gene_tr_dict:
                    gene_tr_dict[line[0]] = []
                if os.path.exists(rep+line[1]+'.fasta'):
                    if line[1].startswith("CGA"):
                        gene_tr_dict[line[0]].append(line[1])
                    
    return gene_tr_dict
        

def gene_to_gtf(repository, gene_tr, gene_triplet_interest):
    '''
    Function description:
        This function use the set of transcripts and their genes used by 
        CG-Alcode programm and return a new file in gtf format.
    
    Input:
        - repository: the repository where are located genes and transcripts
                      files.
        - gene_tr: a file with the correspondance between genes and their 
                   transcripts.
        - gene_triplet_interest: a file with only interested genes.
        
    Output:
        mm_ensembl_cgalcode.gtf: a file in gtf format. For the moment, this fonction is 
        adapted to the species Mus musculus (mouse).
        
    '''
    genes_interested_set = list()
        
    with open(gene_triplet_interest, 'r') as gset_i:
        for gene_sp in gset_i.read().splitlines():
            genes_interested_set.append(gene_sp)
            
    gene_tr_dict = gene_tr_relation(gene_tr)
    
    with open("mm_ensembl_cgalcode.gtf", 'w') as gtf_o:
        for file in os.listdir(repository):
            gene_chr_nb = "" #1
            gene_annot_source = "Ensembl" #2
            gene_entity = "gene" #3
            gene_beg = "" #4
            gene_end = "" #5
            gene_score = "." #6
            gene_strand = "" #7
            gene_phase = "." #8
            #9:
            gene_id = ""
            transcript_nb = ""
            transcript_allbiotypes = "NA"
            ensembl_name = "NA"
            gene_name = "NA"
            gene_biotype = "NA"
                   
            if file.startswith("ENSG"):
                gene = (file)
#                print(os.path.exists(rep+file))
                with open(rep+gene, 'r') as gene_file:
                    lecture_gene = gene_file.readlines()
                    for line_gene in lecture_gene:
                        if line_gene.startswith('>'):
                            line_gene = line_gene.replace('\n', '').split(" ")
                            gene_chr_nb = line_gene[3]
                            gene_beg = line_gene[5]
                            gene_end = line_gene[6]
                            if line_gene[7] == "1":
                                gene_strand = "+"
                            elif line_gene[7] == "-1":
                                gene_strand = "-"
                            gene_id = line_gene[1]
                            if gene_id in genes_interested_set:
                                transcript_nb = str(len(gene_tr_dict[gene_id]))
                                gtf_o.write(gene_chr_nb+
                                            '\t'+gene_annot_source+
                                            '\t'+gene_entity+
                                            '\t'+gene_beg+ 
                                            '\t'+gene_end+
                                            '\t'+gene_score+
                                            '\t'+gene_strand+
                                            '\t'+gene_phase+
                                            "\tgene_id "+'"'+gene_id+'"; '+
                                            'transcript_nb "'+transcript_nb+'"; '+
                                            'transcript_allbiotypes "'+transcript_allbiotypes+'; "'+
                                            'ensembl_name "'+ensembl_name+'"; '+
                                            'gene_name "'+gene_name+'"; '+
                                            'gene_biotype "'+gene_biotype+'";'+'\n')
                                
                                for tr in gene_tr_dict[gene_id]:
                                    if tr.startswith("CGA"):
                                        tr_annot_source = "CGAlcode" #2
                                    else:
                                        tr_annot_source = "Ensembl"
                                    tr_entity = "transcript" #3
                                    tr_beg = "" #4
                                    tr_end = "" #5
                                    tr_score = "." #6
                                    tr_strand = "" #7
                                    tr_phase = "." #8
                                    #9:
                                    tr_id = ""
                                    exon_nb = 0
                                    tr_name = ""
                                    tr_biotype = "protein_coding"
                                    
                                    out_exons = ""
                                    
                                    with open(rep+tr+".fasta", 'r') as tr_file:
                                        lecture_tr = tr_file.readlines()
                                        for line_tr in lecture_tr:
                                            line_tr = line_tr.replace('\n', '')
                                            if line_tr.startswith('>'):
                                                line_tr = line_tr.replace('\n', '').split()
                                                tr_beg = line_tr[3]
                                                tr_end = line_tr[4]
                                                if line_tr[5] == "1":
                                                    tr_strand = "+"
                                                elif line_tr[5] == "-1":
                                                    tr_strand = "-"
                                                tr_id = line_tr[1]
                                                tr_name = line_tr[7]
                                                
                                            else:
                                                exon_entity = "exon"
                                                exon_beg = ""
                                                exon_end = ""
                                                exon_score = "."
                                                line_tr = line_tr.replace('\n', '').split('], [')
                                                for event in line_tr:
                                                    event = event.replace('[', '').replace(']', '').replace("'", '').replace(' ', '').split(',')
                                                    if len(event) == 3:
                                                        if event[2] == "EXON":
                                                            exon_nb += 1
                                                            exon_beg = event[0]
                                                            exon_end = event[1]
                                                                
                                                            out_exons += str(gene_chr_nb+'\t'+tr_annot_source+'\t'+exon_entity+'\t'+exon_beg+'\t'+exon_end+'\t'+exon_score+'\t'+tr_strand+'\t'+tr_phase+"\tgene_id "+'"'+gene_id+'"; '+"transcript_id "+'"'+tr_id+'"; '+'exon_number "'+str(exon_nb)+'"; '+'class_code "U"; '+'transcript_name "'+tr_name+'"; '+'transcript_biotype "'+tr_biotype+'";'+'\n')
                                                
                                        gtf_o.write(gene_chr_nb+
                                                    '\t'+tr_annot_source+
                                                    '\t'+tr_entity+
                                                    '\t'+tr_beg+
                                                    '\t'+tr_end+
                                                    '\t'+tr_score+
                                                    '\t'+tr_strand+
                                                    '\t'+tr_phase+
                                                    "\tgene_id "+'"'+gene_id+'"; '+
                                                    "transcript_id "+'"'+tr_id+'"; '+
                                                    'exon_nb "'+str(exon_nb)+'"; '+
                                                    'class_code "U"; '+
                                                    'transcript_name "'+tr_name+'"; '+
                                                    'transcript_biotype "'+tr_biotype+'";'+'\n')
                                        
                                        gtf_o.write(out_exons)
               
def gene_to_gtf_from_ensembl_only(repository, gene_tr, gene_triplet_interest):
    '''
    Function description:
        This function use the set of Ensembl transcripts and their genes used
        by CG-Alcode programm and return a new file in gtf format and a gtf by
        gene.
    
    Input:
        - repository: the repository where are located genes and transcripts 
                      Ensembl files.
        - gene_tr: a file with the correspondance between genes and their 
                   transcripts from Ensembl.
        - gene_triplet_interest: a file with only interested genes.
        
    Output:
        - mm_ensembl.gtf: a file in gtf format. For the moment, this fonction is 
                          adapted to the species Mus musculus (mouse).
        - gene_id_ensembl.gtf: a gtf file for a gene_id
        
    '''
    genes_interested_set = list()
        
    with open(gene_triplet_interest, 'r') as gset_i:
        for gene_sp in gset_i.read().splitlines():
            genes_interested_set.append(gene_sp)
            
    gene_tr_dict = gene_tr_ensembl_relation(gene_tr)
    
    with open("mm_ensembl.gtf", 'w') as gtf_o_tot:
        for file in os.listdir(repository):
            gene_chr_nb = "" #1
            gene_annot_source = "Ensembl" #2
            gene_entity = "gene" #3
            gene_beg = "" #4
            gene_end = "" #5
            gene_score = "." #6
            gene_strand = "" #7
            gene_phase = "." #8
            #9:
            gene_id = ""
            transcript_nb = ""
            transcript_allbiotypes = "NA"
            ensembl_name = "NA"
            gene_name = "NA"
            gene_biotype = "NA"
                   
            if file.startswith("ENSG"):
                gene = (file)
                gene_id_out = gene.split('.fasta')[0]
                if gene_id_out in genes_interested_set:
                    with open(gene_id_out+"_ensembl.gtf", 'w') as gtf_o:
                        with open(rep+gene, 'r') as gene_file:
                            lecture_gene = gene_file.readlines()
                            for line_gene in lecture_gene:
                                if line_gene.startswith('>'):
                                    line_gene = line_gene.replace('\n', '').split(" ")
                                    gene_chr_nb = line_gene[3]
                                    gene_beg = line_gene[5]
                                    gene_end = line_gene[6]
                                    if line_gene[7] == "1":
                                        gene_strand = "+"
                                    elif line_gene[7] == "-1":
                                        gene_strand = "-"
                                    gene_id = line_gene[1]
                                    transcript_nb = str(len(gene_tr_dict[gene_id]))
                                    if int(transcript_nb) > 0:
                                        gtf_o.write(gene_chr_nb+
                                                    '\t'+gene_annot_source+
                                                    '\t'+gene_entity+
                                                    '\t'+gene_beg+ 
                                                    '\t'+gene_end+
                                                    '\t'+gene_score+
                                                    '\t'+gene_strand+
                                                    '\t'+gene_phase+
                                                    "\tgene_id "+'"'+gene_id+'"; '+
                                                    'transcript_nb "'+transcript_nb+'"; '+
                                                    'transcript_allbiotypes "'+transcript_allbiotypes+'; "'+
                                                    'ensembl_name "'+ensembl_name+'"; '+
                                                    'gene_name "'+gene_name+'"; '+
                                                    'gene_biotype "'+gene_biotype+'";'+'\n')
                                        
                                        gtf_o_tot.write(gene_chr_nb+
                                                        '\t'+gene_annot_source+
                                                        '\t'+gene_entity+
                                                        '\t'+gene_beg+ 
                                                        '\t'+gene_end+
                                                        '\t'+gene_score+
                                                        '\t'+gene_strand+
                                                        '\t'+gene_phase+
                                                        "\tgene_id "+'"'+gene_id+'"; '+
                                                        'transcript_nb "'+transcript_nb+'"; '+
                                                        'transcript_allbiotypes "'+transcript_allbiotypes+'; "'+
                                                        'ensembl_name "'+ensembl_name+'"; '+
                                                        'gene_name "'+gene_name+'"; '+
                                                        'gene_biotype "'+gene_biotype+'";'+'\n')
                                    
                                    for tr in gene_tr_dict[gene_id]:
                                        if not tr.startswith("CGA"):
                                            tr_annot_source = "Ensembl" #2
                                            tr_entity = "transcript" #3
                                            tr_beg = "" #4
                                            tr_end = "" #5
                                            tr_score = "." #6
                                            tr_strand = "" #7
                                            tr_phase = "." #8
                                            #9:
                                            tr_id = ""
                                            exon_nb = 0
                                            tr_name = ""
                                            tr_biotype = "protein_coding"
                                            
                                            out_exons = ""
                                            
                                            with open(rep+tr+".fasta", 'r') as tr_file:
                                                lecture_tr = tr_file.readlines()
                                                for line_tr in lecture_tr:
                                                    line_tr = line_tr.replace('\n', '')
                                                    if line_tr.startswith('>'):
                                                        line_tr = line_tr.replace('\n', '').split()
                                                        tr_beg = line_tr[3]
                                                        tr_end = line_tr[4]
                                                        if line_tr[5] == "1":
                                                            tr_strand = "+"
                                                        elif line_tr[5] == "-1":
                                                            tr_strand = "-"
                                                        tr_id = line_tr[1]
                                                        tr_name = line_tr[7]
                                                        
                                                    else:
                                                        exon_entity = "exon"
                                                        exon_beg = ""
                                                        exon_end = ""
                                                        exon_score = "."
                                                        line_tr = line_tr.replace('\n', '').split('], [')
                                                        for event in line_tr:
                                                            event = event.replace('[', '').replace(']', '').replace("'", '').replace(' ', '').split(',')
                                                            if len(event) == 3:
                                                                if event[2] == "EXON":
                                                                    exon_nb += 1
                                                                    exon_beg = event[0]
                                                                    exon_end = event[1]
                                                                        
                                                                    out_exons = out_exons + str(gene_chr_nb+'\t'+tr_annot_source+'\t'+exon_entity+'\t'+exon_beg+'\t'+exon_end+'\t'+exon_score+'\t'+tr_strand+'\t'+tr_phase+"\tgene_id "+'"'+gene_id+'"; '+"transcript_id "+'"'+tr_id+'"; '+'exon_number "'+str(exon_nb)+'"; '+'class_code "U"; '+'transcript_name "'+tr_name+'"; '+'transcript_biotype "'+tr_biotype+'";'+'\n')
                                                                                                    
                                                gtf_o.write(gene_chr_nb+
                                                            '\t'+tr_annot_source+
                                                            '\t'+tr_entity+
                                                            '\t'+tr_beg+
                                                            '\t'+tr_end+
                                                            '\t'+tr_score+
                                                            '\t'+tr_strand+
                                                            '\t'+tr_phase+
                                                            "\tgene_id "+'"'+gene_id+'"; '+
                                                            "transcript_id "+'"'+tr_id+'"; '+
                                                            'exon_nb "'+str(exon_nb)+'"; '+
                                                            'class_code "U"; '+
                                                            'transcript_name "'+tr_name+'"; '+
                                                            'transcript_biotype "'+tr_biotype+'";'+'\n')
                                                
                                                gtf_o_tot.write(gene_chr_nb+
                                                            '\t'+tr_annot_source+
                                                            '\t'+tr_entity+
                                                            '\t'+tr_beg+
                                                            '\t'+tr_end+
                                                            '\t'+tr_score+
                                                            '\t'+tr_strand+
                                                            '\t'+tr_phase+
                                                            "\tgene_id "+'"'+gene_id+'"; '+
                                                            "transcript_id "+'"'+tr_id+'"; '+
                                                            'exon_nb "'+str(exon_nb)+'"; '+
                                                            'class_code "U"; '+
                                                            'transcript_name "'+tr_name+'"; '+
                                                            'transcript_biotype "'+tr_biotype+'";'+'\n')
                                                
                                                gtf_o.write(out_exons)
                                                
                                                gtf_o_tot.write(out_exons)
                
                
def gene_to_gtf_with_only_prediction(repository, gene_tr, gene_triplet_interest):
    '''
    Function description:
        This function use the set of predicted transcripts and their genes used
        by CG-Alcode programm and return a new file in gtf format.
    
    Input:
        - repository: the repository where are located genes and predicted 
                      transcripts files.
        - gene_tr: a file with the correspondance between genes and their 
                   predicted transcripts.
        - gene_triplet_interest: a file with only interested genes.
        
    Output:
        - mm_cgalcode.gtf: a file in gtf format. For the moment, this fonction is 
                           adapted to the species Mus musculus (mouse).
        - gene_id_prediction.gtf: a gtf file for a gene_id
        
    '''
    genes_interested_set = list()
        
    with open(gene_triplet_interest, 'r') as gset_i:
        for gene_sp in gset_i.read().splitlines():
            genes_interested_set.append(gene_sp)
    
    gene_tr_dict = gene_tr_pred_relation(gene_tr)
    
    with open("mm_cgalcode.gtf", 'w') as gtf_o_tot:
        for file in os.listdir(repository):
            gene_chr_nb = "" #1
            gene_annot_source = "Ensembl" #2
            gene_entity = "gene" #3
            gene_beg = "" #4
            gene_end = "" #5
            gene_score = "." #6
            gene_strand = "" #7
            gene_phase = "." #8
            #9:
            gene_id = ""
            transcript_nb = ""
            transcript_allbiotypes = "NA"
            ensembl_name = "NA"
            gene_name = "NA"
            gene_biotype = "NA"
                   
            if file.startswith("ENSG"):
                gene = (file)
                gene_id_out = gene.split('.fasta')[0]
                if gene_id_out in genes_interested_set:
                    with open(gene_id_out+"_prediction.gtf", 'w') as gtf_o:
                        with open(rep+gene, 'r') as gene_file:
                            lecture_gene = gene_file.readlines()
                            for line_gene in lecture_gene:
                                if line_gene.startswith('>'):
                                    line_gene = line_gene.replace('\n', '').split(" ")
                                    gene_chr_nb = line_gene[3]
                                    gene_beg = line_gene[5]
                                    gene_end = line_gene[6]
                                    if line_gene[7] == "1":
                                        gene_strand = "+"
                                    elif line_gene[7] == "-1":
                                        gene_strand = "-"
                                    gene_id = line_gene[1]
    #                                print(gene_id)
                                    transcript_nb = str(len(gene_tr_dict[gene_id]))
                                    if gene_id in genes_interested_set:
                                        if int(transcript_nb) > 0:
                                            
                                            gtf_o.write(gene_chr_nb+
                                                        '\t'+gene_annot_source+
                                                        '\t'+gene_entity+
                                                        '\t'+gene_beg+ 
                                                        '\t'+gene_end+
                                                        '\t'+gene_score+
                                                        '\t'+gene_strand+
                                                        '\t'+gene_phase+
                                                        "\tgene_id "+'"'+gene_id+'"; '+
                                                        'transcript_nb "'+transcript_nb+'"; '+
                                                        'transcript_allbiotypes "'+transcript_allbiotypes+'; "'+
                                                        'ensembl_name "'+ensembl_name+'"; '+
                                                        'gene_name "'+gene_name+'"; '+
                                                        'gene_biotype "'+gene_biotype+'";'+'\n')
                                        
                                        gtf_o_tot.write(gene_chr_nb+
                                                        '\t'+gene_annot_source+
                                                        '\t'+gene_entity+
                                                        '\t'+gene_beg+ 
                                                        '\t'+gene_end+
                                                        '\t'+gene_score+
                                                        '\t'+gene_strand+
                                                        '\t'+gene_phase+
                                                        "\tgene_id "+'"'+gene_id+'"; '+
                                                        'transcript_nb "'+transcript_nb+'"; '+
                                                        'transcript_allbiotypes "'+transcript_allbiotypes+'; "'+
                                                        'ensembl_name "'+ensembl_name+'"; '+
                                                        'gene_name "'+gene_name+'"; '+
                                                        'gene_biotype "'+gene_biotype+'";'+'\n')
                                        
                                        for tr in gene_tr_dict[gene_id]:
                                            if tr.startswith("CGA"):
                                                tr_annot_source = "CGAlcode" #2
                                                tr_entity = "transcript" #3
                                                tr_beg = "" #4
                                                tr_end = "" #5
                                                tr_score = "." #6
                                                tr_strand = "" #7
                                                tr_phase = "." #8
                                                #9:
                                                tr_id = ""
                                                exon_nb = 0
                                                tr_name = ""
                                                tr_biotype = "protein_coding"
                                                
                                                out_exons = ""
                                                
                                                with open(rep+tr+".fasta", 'r') as tr_file:
                                                    lecture_tr = tr_file.readlines()
                                                    for line_tr in lecture_tr:
                                                        line_tr = line_tr.replace('\n', '')
                                                        if line_tr.startswith('>'):
                                                            line_tr = line_tr.replace('\n', '').split()
                                                            tr_beg = line_tr[3]
                                                            tr_end = line_tr[4]
                                                            if line_tr[5] == "1":
                                                                tr_strand = "+"
                                                            elif line_tr[5] == "-1":
                                                                tr_strand = "-"
                                                            tr_id = line_tr[1]
                                                            tr_name = line_tr[7]
                                                            
                                                        else:
                                                            exon_entity = "exon"
                                                            exon_beg = ""
                                                            exon_end = ""
                                                            exon_score = "."
                                                            line_tr = line_tr.replace('\n', '').split('], [')
                                                            for event in line_tr:
                                                                event = event.replace('[', '').replace(']', '').replace("'", '').replace(' ', '').split(',')
                                                                if len(event) == 3:
                                                                    if event[2] == "EXON":
                                                                        exon_nb += 1
                                                                        exon_beg = event[0]
                                                                        exon_end = event[1]
                                                                            
                                                                        out_exons = out_exons + str(gene_chr_nb+'\t'+tr_annot_source+'\t'+exon_entity+'\t'+exon_beg+'\t'+exon_end+'\t'+exon_score+'\t'+tr_strand+'\t'+tr_phase+"\tgene_id "+'"'+gene_id+'"; '+"transcript_id "+'"'+tr_id+'"; '+'exon_number "'+str(exon_nb)+'"; '+'class_code "U"; '+'transcript_name "'+tr_name+'"; '+'transcript_biotype "'+tr_biotype+'";'+'\n')
                                                            
                                                    gtf_o.write(gene_chr_nb+
                                                                '\t'+tr_annot_source+
                                                                '\t'+tr_entity+
                                                                '\t'+tr_beg+
                                                                '\t'+tr_end+
                                                                '\t'+tr_score+
                                                                '\t'+tr_strand+
                                                                '\t'+tr_phase+
                                                                "\tgene_id "+'"'+gene_id+'"; '+
                                                                "transcript_id "+'"'+tr_id+'"; '+
                                                                'exon_nb "'+str(exon_nb)+'"; '+
                                                                'class_code "U"; '+
                                                                'transcript_name "'+tr_name+'"; '+
                                                                'transcript_biotype "'+tr_biotype+'";'+'\n')
                                                    
                                                    gtf_o_tot.write(gene_chr_nb+
                                                                '\t'+tr_annot_source+
                                                                '\t'+tr_entity+
                                                                '\t'+tr_beg+
                                                                '\t'+tr_end+
                                                                '\t'+tr_score+
                                                                '\t'+tr_strand+
                                                                '\t'+tr_phase+
                                                                "\tgene_id "+'"'+gene_id+'"; '+
                                                                "transcript_id "+'"'+tr_id+'"; '+
                                                                'exon_nb "'+str(exon_nb)+'"; '+
                                                                'class_code "U"; '+
                                                                'transcript_name "'+tr_name+'"; '+
                                                                'transcript_biotype "'+tr_biotype+'";'+'\n')
                                                    
                                                    gtf_o.write(out_exons)
                                                    
                                                    gtf_o_tot.write(out_exons)
                                            
                                
def open_file_gene_tr_relation(gene_tr_file):
    '''
    Function description:
        This function read a file with the list of reference genes and return a
        list of this gene.
        
    Input:
        - gene_file_list: a file with reference genes.
        
    Return:
        - gene_tr_list: a list of reference genes.
    '''
    gene_tr_list = dict()
    
    with open(gene_tr_file, 'r') as ref_gene_f:
        lecture = ref_gene_f.readlines()
        for gene_tr in lecture:
            gene_tr = gene_tr.replace('\n', '').split('\t')
            try:
                gene_tr_list[gene_tr[0]].append(gene_tr[1])
            except KeyError:
                gene_tr_list[gene_tr[0]] = [gene_tr[1]]
#            gene_tr_list[gene_tr]
#            print(gene_tr)
#            gene_ref_list.append(gene)
        
    return gene_tr_list
                                
def gene_to_gtf_with_only_prediction_on_a_set_of_interest(repository, gene_tr, gene_tr_interest):
    '''
    Function description:
        This function use the set of predicted transcripts and their genes used
        by CG-Alcode programm and return a new file in gtf format. This 
        function uses a genes/transcripts correspondance of a specific set of
        interest.
    
    Input:
        - repository: the repository where are located genes and predicted 
                      transcripts files.
        - gene_tr: a file with the correspondance between genes and their 
                   predicted transcripts.
        - gene_tr_interest: a file with correspondance between genes and their 
                            predicted transcripts on a specific set of interest.
        
    Output:
        - cga_caf.gtf: a file in gtf format to a set of interest genes/
                     transcripts relationship. For the moment, this fonction is 
                     adapted to the species Canis lupus familiaris (dog).
    '''
    gene_tr_dict = gene_tr_pred_relation(gene_tr)
    
    gene_tr_interet = open_file_gene_tr_relation(gene_tr_interest)
    
    with open("clf_cgalcode_set_interest.gtf", 'w') as gtf_o:
        for file in os.listdir(repository):
            gene_chr_nb = "" #1
            gene_annot_source = "Ensembl" #2
            gene_entity = "gene" #3
            gene_beg = "" #4
            gene_end = "" #5
            gene_score = "." #6
            gene_strand = "" #7
            gene_phase = "." #8
            #9:
            gene_id = ""
            transcript_nb = ""
            transcript_allbiotypes = "NA"
            ensembl_name = "NA"
            gene_name = "NA"
            gene_biotype = "NA"
                   
            if file.startswith("ENSCAFG"):
                gene = (file)
                with open(rep+gene, 'r') as gene_file:
                    lecture_gene = gene_file.readlines()
                    for line_gene in lecture_gene:
                        if line_gene.startswith('>'):
                            line_gene = line_gene.replace('\n', '').split(" ")
                            gene_chr_nb = line_gene[3]
                            gene_beg = line_gene[5]
                            gene_end = line_gene[6]
                            if line_gene[7] == "1":
                                gene_strand = "+"
                            elif line_gene[7] == "-1":
                                gene_strand = "-"
                            gene_id = line_gene[1]
                            try:
                                transcript_nb = str(len(gene_tr_interet[gene_id]))
                            except KeyError:
                                transcript_nb = "0"
                            if int(transcript_nb) > 0:
                                gtf_o.write(gene_chr_nb+
                                            '\t'+gene_annot_source+
                                            '\t'+gene_entity+
                                            '\t'+gene_beg+ 
                                            '\t'+gene_end+
                                            '\t'+gene_score+
                                            '\t'+gene_strand+
                                            '\t'+gene_phase+
                                            "\tgene_id "+'"'+gene_id+'"; '+
                                            'transcript_nb "'+transcript_nb+'"; '+
                                            'transcript_allbiotypes "'+transcript_allbiotypes+'; "'+
                                            'ensembl_name "'+ensembl_name+'"; '+
                                            'gene_name "'+gene_name+'"; '+
                                            'gene_biotype "'+gene_biotype+'";'+'\n')
                            
                                for tr in gene_tr_dict[gene_id]:
                                    if tr in gene_tr_interet[gene_id]:
                                        if tr.startswith("CGA"):
                                            tr_annot_source = "CGAlcode" #2
                                            tr_entity = "transcript" #3
                                            tr_beg = "" #4
                                            tr_end = "" #5
                                            tr_score = "." #6
                                            tr_strand = "" #7
                                            tr_phase = "." #8
                                            #9:
                                            tr_id = ""
                                            exon_nb = 0
                                            tr_name = ""
                                            tr_biotype = "protein_coding"
                                            
                                            out_exons = ""
                                            
                                            with open(rep+tr+".fasta", 'r') as tr_file:
                                                lecture_tr = tr_file.readlines()
                                                for line_tr in lecture_tr:
                                                    line_tr = line_tr.replace('\n', '')
                                                    if line_tr.startswith('>'):
                                                        line_tr = line_tr.replace('\n', '').split()
                                                        tr_beg = line_tr[3]
                                                        tr_end = line_tr[4]
                                                        if line_tr[5] == "1":
                                                            tr_strand = "+"
                                                        elif line_tr[5] == "-1":
                                                            tr_strand = "-"
                                                        tr_id = line_tr[1]
                                                        tr_name = line_tr[7]
                                                        
                                                    else:
                                                        exon_entity = "exon"
                                                        exon_beg = ""
                                                        exon_end = ""
                                                        exon_score = "."
                                                        line_tr = line_tr.replace('\n', '').split('], [')
                                                        for event in line_tr:
                                                            event = event.replace('[', '').replace(']', '').replace("'", '').replace(' ', '').split(',')
                                                            if len(event) == 3:
                                                                if event[2] == "EXON":
                                                                    exon_nb += 1
                                                                    exon_beg = event[0]
                                                                    exon_end = event[1]
                                                                        
                                                                    out_exons += str(gene_chr_nb+'\t'+tr_annot_source+'\t'+exon_entity+'\t'+exon_beg+'\t'+exon_end+'\t'+exon_score+'\t'+tr_strand+'\t'+tr_phase+"\tgene_id "+'"'+gene_id+'"; '+"transcript_id "+'"'+tr_id+'"; '+'exon_number "'+str(exon_nb)+'"; '+'class_code "U"; '+'transcript_name "'+tr_name+'"; '+'transcript_biotype "'+tr_biotype+'";'+'\n')
                                                        
                                                gtf_o.write(gene_chr_nb+
                                                            '\t'+tr_annot_source+
                                                            '\t'+tr_entity+
                                                            '\t'+tr_beg+
                                                            '\t'+tr_end+
                                                            '\t'+tr_score+
                                                            '\t'+tr_strand+
                                                            '\t'+tr_phase+
                                                            "\tgene_id "+'"'+gene_id+'"; '+
                                                            "transcript_id "+'"'+tr_id+'"; '+
                                                            'exon_nb "'+str(exon_nb)+'"; '+
                                                            'class_code "U"; '+
                                                            'transcript_name "'+tr_name+'"; '+
                                                            'transcript_biotype "'+tr_biotype+'";'+'\n')
                                                
                                                gtf_o.write(out_exons)
                
def create_individual_gtf_from_a_ref_gtf(ref_gtf_file):
    '''
    Function description:
        This function use a reference gtf file and split by gene.
        
    Input: 
        - ref_gtf_file: a gtf reference file to split by gene.
        
    Output:
        - gene_id_reference.gtf: gtf file individual to each gene.
    '''
    with open(ref_gtf_file, 'r') as gtf_ref_i:
        lecture = gtf_ref_i.readlines()
        
        multi_ensembl_name = list()
        
        gene_chr_nb = "" #1
        gene_annot_source = "" #2
        gene_entity = "" #3
        gene_beg = "" #4
        gene_end = "" #5
        gene_score = "" #6
        gene_strand = "" #7
        gene_phase = "" #8
        #9
        gene_id = ""
        transcript_nb = 0
        transcript_allbiotypes = "NA"
        ensembl_name = "NA"
        gene_name = "NA"
        gene_biotype = "NA"
        
        out = ""
        
        for line_gtf in lecture:
            if "protein_coding" in line_gtf and not 'ensembl_name "NA"' in line_gtf:
                line = line_gtf.split('\t')
                
                if line[2] == "gene":
                    
                    if len(multi_ensembl_name) != 0 and transcript_nb != 0:
                        for gene_id in multi_ensembl_name:
                            with open(gene_id+"_reference.gtf", 'w') as gene_ref_gtf_o:
                                gene_ref_gtf_o.write(gene_chr_nb+
                                            '\t'+gene_annot_source+
                                            '\t'+gene_entity+
                                            '\t'+gene_beg+ 
                                            '\t'+gene_end+
                                            '\t'+gene_score+
                                            '\t'+gene_strand+
                                            '\t'+gene_phase+
                                            "\tgene_id "+'"'+gene_id+'"; '+
                                            'transcript_nb "'+str(transcript_nb)+'"; '+
                                            'transcript_allbiotypes "'+transcript_allbiotypes+'; "'+
                                            'ensembl_name "'+gene_id+'"; '+
                                            'gene_name "'+gene_name+'"; '+
                                            'gene_biotype "'+gene_biotype+'";'+'\n')
                                gene_ref_gtf_o.write(out)
                                
                                out = ""
                    
                    out = ""
                    transcript_nb = 0
                    
                    gene_chr_nb = line[0] #1
                    gene_annot_source = line[1] #2
                    gene_entity = line[2] #3
                    gene_beg = line[3] #4
                    gene_end = line[4] #5
                    gene_score = line[5] #6
                    gene_strand = line[6] #7
                    gene_phase = line[7] #8
                    #9
                    gene_id = ""
                
                    transcript_allbiotypes = "protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding"
                    ensembl_name = "NA"
                    gene_name = "NA"
                    gene_biotype = "protein_coding"
                    
                    caracteristics = line[8].replace('\n', '').replace(' ', '').replace('"', ';').replace(';;', ';').split(';') #9
                    if ',' in caracteristics[1]:
                        multi_ensembl_name = caracteristics[1].split(',')
                    else:
                        multi_ensembl_name = [caracteristics[1]]
                    
                else:
                    if line[2] == "transcript":
                        transcript_nb += 1
                    out += line_gtf
                
def correspondance_gtf_ref_test(repository):
    '''
    Function description:
        This function use a repository where gtf file are conserved and check
        if a gtf tested has a correspondance with a gtf reference. If not, gtf
        are removed.
        
    Input:
        - rep: the repository where gtf are conserved.
    '''
    
    file_list_ref = list()
    file_list_pred = list()
    file_list_ens = list()
    
    file_ref_conserved = 0
    file_pred_conserved = 0
    file_ens_conserved = 0
    
    file_list_ref_conserved = list()
    file_list_pred_conserved = list()
    file_list_ens_conserved = list()
    
    for file in os.listdir(repository):
        if not "_reference" in file and not "_ensembl" in file and not "gtf2bed.pl" in file: 
            if not file in file_list_pred and not ".py" in file and not "mm_cgalcode" in file:
                file_list_pred.append(file)
        else:
            if not "_ensembl" in file and file not in file_list_ref and not "gtf2bed.pl" in file:
                file_list_ref.append(file)
            else:
                if not file in file_list_ens and not "mm_ensembl" in file and not "gtf2bed.pl" in file:
                    file_list_ens.append(file)
    
    for file_ref in file_list_ref:
        if os.path.exists(file_ref[0:15]+'_prediction.gtf'):
            file_ref_conserved += 1
            file_list_ref_conserved.append(file_ref)
#            print(file_ref)
        else:
            os.remove(file_ref)
            pass
        
    for file_pred in file_list_pred:
        if os.path.exists(file_pred[0:15]+'_reference.gtf'):
            file_pred_conserved += 1
            file_list_pred_conserved.append(file_pred)
        else:
            os.remove(file_pred)
            pass
        
    for file_ens in file_list_ens:
        if os.path.exists(file_ens[0:15]+'_prediction.gtf'):
#            print(file_ens)
            file_ens_conserved += 1
            file_list_ens_conserved.append(file_ens)
        else:
            os.remove(file_ens)
            pass
        
    with open("correspondance_gene_pred_ref.txt", 'w') as correspondance:
        for file_ref in file_list_ref:
            out = file_ref
            for file_pred in file_list_pred:
    #            print(file_pred[:-4])
                if file_ref.startswith(file_pred[0:15]):
                    out += '\t'+file_pred
                    correspondance.write(file_ref+'\t'+file_pred+'\n')
                    
#            for file_ens in file_list_ens:
#                if file_ref.startswith(file_ens[0:18]):
#                    out += '\t'+file_ens
#            correspondance.write(out+'\n')
                    
    nb_tr = 0
    print(len(file_list_ens_conserved))
    for file_ref_conserved_file in file_list_ref_conserved:
#        print(os.path.exists(file_ref_conserved), file_ref_conserved)
        with open(file_ref_conserved_file, 'r') as f:
            for line in f.read().splitlines():
                if '\ttranscript\t' in line:
                    nb_tr += 1
#                    pass
#                    print(line)
                else:
                    if '\tgene\t' in line:
#                        print(line)
                        pass
##                if '\ttranscript\t' in line:
##                    print(line)
    
    print("Number of genes with predicted transcripts conserved:", str(file_pred_conserved)+"/"+str(len(file_list_pred)))
    print("Number of ref genes conserved:", str(file_ref_conserved)+"/"+str(len(file_list_ref)))
    print("Number of ens genes conserved:", str(file_ens_conserved))
    print("Number of ref transcripts conserved:", str(nb_tr))
    
    return

def gtf2bed(file, utr=False):
    '''
    Function description:
        This function convert a gtf file into bed file with 12 columns (bed12)
        with utr if exists (.tmp) and without utr (.bed).
    '''    
    path_to_script = "./gtf2bed.pl"
    
    if utr == True:
        out = file+".bed"
        cmd = "{0} {1}  >> {2}".format(path_to_script, file, out)
        subprocess.call(cmd, shell=True)
    
    if utr == False:
#        if "_reference" in file:
        out = file+".tmp"
#        else:
#            out = file+".bed"
        cmd = "{0} {1}  >> {2}".format(path_to_script, file, out)
        subprocess.call(cmd, shell=True)
#        if ".tmp" in out:
            
#            file_out = out+".bed"
#            path = "awk '{OFS='\t';split($11,a,','); split($12,b,','); A=''; B=''; if($7==$8) next; j=0; for(i=1;i<length(a);i++) if(($2+b[i]+a[i])>$7 && ($2+b[i])<$8) {j++; start=$2+b[i]-$7; size=a[i]; if(($2+b[i])<=$7) {start=0;size=size-($7-($2+b[i]));} if(($2+a[i]+b[i])>=$8) {size=size-($2+a[i]+b[i]-$8);} A=A''size',';B=B''start',';} print $1,$7,$8,$4,$5,$6,$7,$8,$9,j,A,B;}'"
#            cmd = "{0} >> {1}".format(path, file_out)
#            subprocess.call(cmd, shell = True)
    #         cmd = "awk '{OFS='\t';split($11,a,','); split($12,b,','); A=''; B=''; if($7==$8) next; j=0; for(i=1;i<length(a);i++) if(($2+b[i]+a[i])>$7 && ($2+b[i])<$8) {j++; start=$2+b[i]-$7; size=a[i]; if(($2+b[i])<=$7) {start=0;size=size-($7-($2+b[i]));} if(($2+a[i]+b[i])>=$8) {size=size-($2+a[i]+b[i]-$8);} A=A''size',';B=B''start',';} print $1,$7,$8,$4,$5,$6,$7,$8,$9,j,A,B;}' out > file.bed"
        with open(out, 'r') as file_tmp:
            with open(file+".bed", 'w') as file_bed:
                for line in file_tmp.read().splitlines():
                    line = line.split('\t')
#                        print(line)
                    a = line[10].split(',')
#                        print(a)
                    b = line[11].split(',')
                    A = ""
                    B = ""
                    size = ""
                    start = ""
#                        
                    if line[6] == line[7]:
                        pass
                    j=0
                    for i in range(len(a)-1):
                        if (int(line[1])+int(b[i])+int(a[i])) > int(line[6]) and (int(line[1])+int(b[i])) < int(line[7]):
                            j+=1
                            start = int(line[1])+int(b[i])-int(line[6])
                            size = int(a[i])
                            if (int(line[1])+int(b[i]) <= int(line[6])):
                                start = 0
                                size = size-(int(line[6])-(int(line[1])+int(b[i])))
                        if (int(line[1])+int(a[i])+int(b[i])) >= int(line[7]):
                            size = size-(int(line[1])+int(a[i])+int(b[i])-int(line[7]))
                        A = A+''+str(size)+','
                        B = B+''+str(start)+','
                        if j == int(line[9]):
                            new_line = line[0]+'\t'+line[6]+'\t'+line[7]+'\t'+line[3]+'\t'+line[4]+'\t'+line[5]+'\t'+line[6]+'\t'+line[7]+'\t'+line[8]+'\t'+str(j)+'\t'+str(A)+'\t'+str(B)
                            file_bed.write(new_line+'\n')
                                
                            
                    
#                    lecture = file_tmp.readlines()
#    #                print(lecture)
#                    for line in lecture:
##                        print(line)
#                        line = line.split()
##                        if line[3] == "ENST00000644493":
#                        new_coord = ""
#                        block_length = line[10].split(',')
#                        block_position = line[11].split(',')
#                        print(block_length)
#                        print(block_position)
#                        position_correction_first_block = int(line[6])-int(line[1])
#                        position_correction_last_block = int(line[2])-int(line[7])
#                        
##                        print("Nb of block:", len(block_length))
#                        length_test = len(block_length)
##                        print(type(length_test))
#                        
#                        #Correction 1st and last block length -> line[10]
#                        len_list_without_utr = ""
#                        for block_length_index in range(len(block_length)-1):
#                            print(block_length_index)
#                            print(len_list_without_utr)
#                            if block_length_index == 0:
#                                print(int(block_length[block_length_index]))
#                                print(int(position_correction_first_block))
#                                len_list_without_utr += str(int(block_length[block_length_index])-int(position_correction_first_block))+','
#                            else:
#                                if block_length_index < len(block_length) - 2 :
#    #                                print(int(block_length[i]))
#                                    len_list_without_utr += str(int(block_length[block_length_index]))+','
#                                else:
#    #                                print(int(block_length[i]))
#                                    len_list_without_utr += str(int(block_length[block_length_index])-int(position_correction_last_block))+','
#                        print(len_list_without_utr)
#                        
#                        #Correction block positions
#                        pos_list_without_utr = ""
#    #                    for 
#                        
#                        
#                        
#                        print("1st block correction:", position_correction_first_block)
#                        print(int(line[6]), int(line[1]))
#                        print("Last block correction:", position_correction_last_block)
#                        print(int(line[2]), int(line[7]))
#                        new_line = line[0]+'\t'+line[6]+'\t'+line[7]+'\t'+line[3]+'\t'+line[4]+'\t'+line[5]+'\t'+line[6]+'\t'+line[7]+'\t'+line[8]+'\t'+line[9]+'\t'+len_list_without_utr+'\t'+line[11]+'\n'
#                        file_bed.write(new_line)
                    
#                        break
    
def bed2fasta(file):
    '''
    Function description:
        This function convert a bed file into fasta file.
    '''
    chr_rep = "/home/niguilla/Documents/these_nguillaudeux/Data/CanFam3.1_assembly/"
    chromosome = ""
    with open(file, 'r') as file:
        lecture = file.readlines()
        for line in lecture:
            line = line.split('\t')
            chromosome = line[0]
#            break
#    print("chromosome:", chromosome)
    path_to_script = "bedtools getfasta"
    file_name=str(file).split()[1].split('.')[0][6:]
#    print(file_name)
    chr_file = chr_rep+"canfam3.1_chr"+chromosome+".fasta"
#    print(chr_file)
#    print(os.path.exists(chr_file))
    out_file = file_name+".fasta"
#    print(out_file)
#    print(path_to_script, chr_file, file_name, out_file)
#    print(str(file).split()[1][6:-1])
    bed_file = file_name+".gtf.bed"
#    print(bed_file)
    cmd = "{0} -fi {1} -bed {2} -split -s -name -fo {3}".format(path_to_script, chr_file, bed_file, out_file)
    subprocess.call(cmd, shell=True)
    to_del = chr_file+".fai"
    cmd = "rm {0}".format(to_del)
    subprocess.call(cmd, shell=True)

def analyseCorrespondanceBetweenBedWithBedTools(input_file):
    '''
    Function description:
        This function run a BedTools analysis betwenne two gtf genes files. We 
        use a file which contains a correspondance between a gene reference
        file and a gene test file. The file format is a tabular separator.
        This function convert GTF file into BED12 and FASTA file. BED12 
        contains exon information for each transcript. 
    '''
#    input_file = "/home/niguilla/Documents/these_nguillaudeux/Scripts/gtfGeneratorMouse/correspondance_gene_pred_ref.txt"
    path_to_script = "bedtools intersect"
    
    with open(input_file, 'r') as f:
        for line in f.read().splitlines():
            arg_a, arg_b = line.split("\t")
            #To convert GTF to BED12:
            gtf2bed(arg_a, False)
            gtf2bed(arg_b, True)
            #To convert BED to FASTA
            ## convert ref BED12
            bed_file_a = arg_a+".bed"
#            bed2fasta(bed_file_a)
            ## convert predict BED12
            bed_file_b = arg_b+".bed"
#            bed2fasta(bed_file_b)
            
            path_to_script = "bedtools intersect"
#            print(path_to_script)
#            cmd = "echo {0} {1} >> test.txt".format(bed_file_b, bed_file_a)
#            subprocess.call(cmd, shell=True)
            cmd = "{0} -a {1} -b {2} -f 1 -r -wa -wb -split >> result_bedtools_intersect_test.out".format(path_to_script, bed_file_b, bed_file_a)
            
#            bed_file_a = arg_a+".bed"
#            bed_file_b = arg_b+".bed"
#            #To realyze the line command
#            cmd = "{0} -a {1} -b {2} -f 1 -r -wa -wb -split >> result_bedtools_intersect.out".format(path_to_script, bed_file_b, bed_file_a)
#            break
            #To execute the line command
            subprocess.call(cmd, shell=True)
#            break
#            exit()
    
def analyseResultBedToolsOutput(result_output):
    '''
    Function description:
        This function check if transcripts predicted are found in reference 
        genes gtf file. A predicted transcript is valid if it is found with 
        BedTools intersect and if all of these exons also found too.
        
    Input:
        - result_output: the result file obtains after the run of BedTools 
                         intersect.
    
    Output:
        - result_final_bedtools_test.gtf: the result of analyse. The file 
                                          contain transcripts which are found 
                                          in a reference gtf file and when all 
                                          of their exons are also found.
    '''
    
    with open(result_output, 'r') as f:
#    with open("test.gtf", 'r') as f:
        with open("result_final_bedtools_test.gtf", 'w') as o:
            lecture = f.readlines()
            
            tr_ref_source = ""
            transcript_id = ""#"DDHD1_6"
            tr_exon_nb = 0
            nb_exon = 0
            out = ""
            
            for line in lecture:
                line = line.replace('\n', '')
                line_read = line.split('\t')
                if line_read[2] == "transcript":
                    out = ""
#                    print(line_read[10])
                    #source of original transcript reference
                    tr_ref_source = line_read[10] 
                    #number of known exons for this predicted transcript
                    line_tr_pred_caracteristics = line_read[8].replace('"', '').replace('; ', ';').replace(' ', ';').split(';')
                    tr_exon_nb = line_tr_pred_caracteristics[5]
                    #transcript id of predicted transcript
                    tr_pred_id = line_tr_pred_caracteristics[3]
                    #number of known exons for the reference transcript
                    line_tr_ref_caracteristics = line_read[17].replace('"', '').replace('; ', ';').replace(' ', ';').split(';')
                    tr_exon_ref_nb = line_tr_ref_caracteristics[5]
                    #reference transcript id
                    transcript_id = line_tr_ref_caracteristics[3]
                    #Control: check if predicted tr and reference tr are
                    #the same number of exons
                    if tr_exon_nb == tr_exon_ref_nb:
                            out = line+'\n'
                    nb_exon = 0
                    
                if line_read[2] == "exon":
                    #source of original exon reference
                    exon_ref_source = line_read[10]
                    # Control: test if the reference source is similar between 
                    #transcript and exon.
                    if tr_ref_source == exon_ref_source:
                        #the number associate to a predict exon
                        line_exon_pred_caracteristics = line_read[8].replace('"', '').replace('; ', ';').replace(' ', ';').split(';')
                        exon_pred_number = line_exon_pred_caracteristics[5]
                        #exon id of predicted transcript
                        exon_pred_tr_id = line_exon_pred_caracteristics[3]
                        #the number associate to a reference exon
                        line_exon_ref_caracteristics = line_read[17].replace('"', '').replace('; ', ';').replace(' ', ';').split(';')
                        exon_ref_number = line_exon_ref_caracteristics[5]
                        transcript_exon_id = line_exon_ref_caracteristics[3]
                        # Control: check if the numer associate to a predict 
                        #exon is similar to a reference exon.
#                        if exon_pred_number == exon_ref_number:
                        # Control: check if the transcript id is the same
                        #between transcript and exon
                        if transcript_id == transcript_exon_id:
                            # Control: check id ensembl name is similar
                            #between transcript and exon
                            if exon_pred_tr_id == tr_pred_id:
                                # Control: check if output is not empty.
                                if out != "":
                                    nb_exon += 1
                                    out += line+'\n'
                
                if int(tr_exon_nb) == int(nb_exon) and int(nb_exon) > 0:
                    o.write(out)
                    tr_exon_nb = 0
                    nb_exon = 0
                    

def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="")
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-rep", "--repository", 
                        default = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/mouse/", 
                        dest = "repository", 
                        metavar = "DATA REPOSITORY", 
                        help = "Enter the path repository where are located genes and transcripts files used by CG-Alcode programm.")
    parser.add_argument("-c", "--correspondance", 
                        default = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/0_total_genes_exact_without_cn_prediction-hs-mm-clf-2019-03-13/transcripts_with_prediction/transcripts_mm_withprediction.csv", 
                        dest = "correspondance", 
                        metavar = "GENES TRANSCRIPTS CORRESPONDANCE", 
                        help = "Give the file with the description relationship between genes and their transcripts.")
    parser.add_argument("-gti", "--genetrinterest", 
                        default= "/home/niguilla/Documents/these_nguillaudeux/Data/Resultats_finaux/gene_tr_to_test_read.txt", 
                        dest = "gTrInterest", 
                        metavar = "GENES AND THEIR TRANSCRIPTS INTERSTED", 
                        help = "Give an set of genes and their transcripts interested to create gtf on this informations.")
    parser.add_argument("-set", "--setgene",
                        default="/home/niguilla/Documents/these_nguillaudeux/set_genes/c135_mm.txt",
                        dest = "setgene",
                        metavar = "SET OF INTERESTED GENES")
    parser.add_argument("-ref", "--reference",
                        default="/home/niguilla/Documents/these_nguillaudeux/Data/gtf_mm/ensembl/Mus_musculus.GRCm38.96.chr.gtf",
                        dest = "reference",
                        metavar = "GTF REFERENCE FILE")
    
    
    return parser.parse_args()                     

if __name__ == "__main__":
    OPTIONS = parse()
    
    rep = OPTIONS.repository#+'/'
    gene_tr = OPTIONS.correspondance
#    gene_tr_interest = OPTIONS.gTrInterest
    gene_triplet_interest = OPTIONS.setgene
    ref_gtf_file = OPTIONS.reference
    
    #To create a gtf file with all information (Ensembl, prediction)
    gene_to_gtf(rep, gene_tr, gene_triplet_interest)
#    gene_to_gtf(rep, gene_tr)
    
    # To create gtf files (one by gene and a gtf with the total information)
    gene_to_gtf_from_ensembl_only(rep, gene_tr, gene_triplet_interest)
    
#    # To delete empty file
#    nb_empty_file = 0
#    for file in os.listdir():
#        if os.path.getsize(file) == 0:
#            nb_empty_file += 1
#            os.remove(file)
#    print("Number of empty gene Ensembl:", nb_empty_file)    

#    gene_to_gtf_with_only_prediction_on_a_set_of_interest(rep, gene_tr, gene_tr_interest)
    
    # To create gtf files (one by gene and a gtf with the total information)
    gene_to_gtf_with_only_prediction(rep, gene_tr, gene_triplet_interest)
    # To delete empty file
    nb_empty_file = 0
    for file in os.listdir():
        if os.path.getsize(file) == 0:
            nb_empty_file += 1
#            os.remove(file)
    print("Number of empty gene without prediction:", nb_empty_file)
#            
#    print("Number of empty gene (Ensembl or prediction):", nb_empty_file)
#     
    
    # To create a gtf by gene from a reference gtf
#    ref_gtf_file = "/home/niguilla/Documents/canfam3.2/canfam3_cons_annot_TritouPuppetMasterChief_23-03-2016_lncClasse_geneBiotype_withEnsCds_withSilicoCds_withGoodGeneSource.gtf"
#    ref_gtf_file = "/home/niguilla/Documents/these_nguillaudeux/Data/gtf_hs/ensembl/Homo_sapiens.GRCh38.96.chr.gtf"
    create_individual_gtf_from_a_ref_gtf(ref_gtf_file)
    
    # To delete empty file
    nb_empty_file = 0
    for file in os.listdir():
        if os.path.getsize(file) == 0:
            nb_empty_file += 1
#            os.remove(file)
    print("Number of empty gene reference:", nb_empty_file)
        
    # WARNING -> 15 for human (ENSG) and 18 for other species (ENSxxxG)
    correspondance_gtf_ref_test(os.getcwd()) 
    
    analyseCorrespondanceBetweenBedWithBedTools("correspondance_gene_pred_ref.txt")
    
#    #TEST NUMBER OF REFERENCE TRANSCRIPTS
#    with open("correspondance_gene_pred_ref.txt", 'r') as file_correspondance:
#        count = 0
#        lecture = file_correspondance.readlines()
#        for line in lecture:
#            line = line.split()
#            print(line[0]+".bed")
##            print(os.path.exists(line[0]+'.bed'))
#            count += len(open(line[0]+'.bed').readlines(  ))
#            print(count)
#            
#            with open(line[0]+'.bed', 'r') as file_bed:
#                lecture_bed = file_bed.readlines()
#                for line_bed in lecture_bed:
#                    line_bed = line_bed.replace('\r', '')
                    
#
#    #IF WE USE BEDTOOLS INTERSECT WITH GTF FORMAT FILES
### #   analyseResultBedToolsOutput("/home/niguilla/Documents/these_nguillaudeux/Scripts/gtfGenerator/test_real.gtf")
##    analyseResultBedToolsOutput("/home/niguilla/Documents/these_nguillaudeux/Scripts/gtfGenerator/result_bedtools_intersect.out")
### #   analyseResultBedToolsOutput("/home/niguilla/Documents/canfam3.2/test_bedtools_total_genes.out")
