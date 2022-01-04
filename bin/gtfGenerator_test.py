#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 11:57:38 2020

@author: niguilla
"""

import os
import argparse

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

def gene_to_gtf_all_genes(repository, gene_tr):
    '''
    Function description:
        This function uses the set of transcripts and their genes used by 
        CG-Alcode programm and return a new file in gtf format.
    
    Input:
        - repository: the repository where are located genes and transcripts
                      files.
        - gene_tr: a file with the correspondance between genes and their 
                   transcripts.
        
    Output:
        species_ensembl_cgalcode.gtf: a file in gtf format.
        
    '''
    species = repository.split('/')[-2]
    
    gene_tr_dict = gene_tr_relation(gene_tr)
    
    with open(species+"_ensembl_cgalcode_tot.gtf", 'w') as gtf_o:
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
                   
            if "ENS" in file and "G00" in file:#.startswith("ENSG"):
                gene = (file)
                with open(repository+gene, 'r') as gene_file:
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
                            if gene_id in gene_tr_dict:
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
                                    
                                    with open(repository+tr+".fasta", 'r') as tr_file:
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
                                        
def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="")
    # Create command-line parser for all options and arguments to give
    
    parser.add_argument("-rep", "--repository", 
#                        default = "../../examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/mouse/", 
                        #default = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/human/", 
                        dest = "repository", 
                        metavar = "DATA REPOSITORY", 
                        help = "Enter the path repository where are located genes and transcripts files used by CG-Alcode programm.")
    
    parser.add_argument("-c", "--correspondance", 
#                        default = "../../examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/transcripts_with_prediction/transcripts_mm_withprediction.csv",
                        #default = "/home/niguilla/Documents/these_nguillaudeux/Data/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/0_total_genes_exact_without_cn_prediction-hs-mm-clf-2019-03-13/transcripts_with_prediction/transcripts_hs_withprediction.csv", 
                        dest = "correspondance", 
                        metavar = "GENES TRANSCRIPTS CORRESPONDANCE", 
                        help = "Give the file with the description relationship between genes and their transcripts.")
    
#    parser.add_argument("-gti", "--genetrinterest", 
#                        default= "/home/niguilla/Documents/these_nguillaudeux/Data/Resultats_finaux/gene_tr_to_test_read.txt", 
#                        dest = "gTrInterest", 
#                        metavar = "GENES AND THEIR TRANSCRIPTS INTERSTED", 
#                        help = "Give an set of genes and their transcripts interested to create gtf on this informations.")
    
    parser.add_argument("-set", "--setgene",
#                        default = "../../examples/set_genes/c135_mouse.txt",
                        #default="/home/niguilla/Documents/these_nguillaudeux/set_genes/c135_human.txt",
                        dest = "setgene",
                        metavar = "SET OF INTERESTED GENES",
                        help = "give a list of interested genes.")
    
    parser.add_argument("-ref", "--reference",
#                        default = "../../examples/reference_gtf/human/ensembl/Ensembl98/Homo_sapiens.GRCh38.98.gtf",
#                        default = "../../examples/reference_gtf/mouse/ensembl/Ensembl96/Mus_musculus.GRCm38.96.gtf",
                        dest = "reference",
                        metavar = "GTF REFERENCE FILE",
                        help = "give a reference gtf file to validate prediction realize with CG-Alcode programm.")
    
    
    return parser.parse_args()                     

if __name__ == "__main__":
    
    print("########################################")
    print("########## PROGRAMM BEGINNING ##########")
    print("")
    
    OPTIONS = parse()
    
    ###########################################################################
    ### To check input data
    
    rep = OPTIONS.repository#+'/'
    if os.path.exists(rep):
        print(">>> Directory containing genes and transcripts files:", os.path.exists(rep))
    else:
        print(">>> Directory containing genes and transcripts files:", os.path.exists(rep))
        print("############ PROGRAMM STOP #############")
        print("########################################")
        exit()
        
    gene_tr = OPTIONS.correspondance
    if os.path.exists(gene_tr):
        print(">>> File containing genes and transcripts correspondance:", os.path.exists(gene_tr))
    else:
        print(">>> File containing genes and transcripts correspondance:", os.path.exists(gene_tr))
        print("############ PROGRAMM STOP #############")
        print("########################################")
        exit()
    
        
    ###########################################################################
    print("")
    ###########################################################################
    ### To create a gtf file with all information (Ensembl, prediction) 
    ### from CG-Alcode output
    print(">>> Creation: a gtf file with all information (Ensembl, prediction)")
    gene_to_gtf_all_genes(rep, gene_tr)
    ###########################################################################