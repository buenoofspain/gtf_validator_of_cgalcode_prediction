#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 11:37:05 2020

@author: Nicolas Guillaudeux

This script create a table from a set of gene with all informations.
"""
import urllib.request
import requests
import time
from bs4 import BeautifulSoup

import os
import csv

#from urllib.parse import urlparse
#from bs4 import BeautifulSoup

def listOfInterestedGene(file):
    """
    Function description:
        This function extract the list of ensembl id genes contains in a file 
        to give a list.
        
    Input:
        - file: a list with ensembl id genes.
        
    Return:
        - listGenes: a list of ensembl id genes.
    """
    listGenes = list()
    with open(file, 'r') as f_i:
        lecture = f_i.readlines()
        listGenes = [ensemblId.replace('\n', '') for ensemblId in lecture]

    return listGenes

def geneNameFromEnsemblId(ensembl_id, url):
    """
    Function description:
        This function uses ensembl id to obtain gene id and gene name.
        
    Input:
        - ensembl_id: a file with all ensembl id of interest
        - url: the url of the web site to obtain information
    
    Return:
        - dict_geneName_ensemblId: a dictionary to store information
    """
#    data = (urllib.request.urlopen(url).read())
#    data = (urllib.request.urlopen(url).read().decode('utf-8'))
#    print(str(data).split("\n   "))
    
    dict_geneName_ensemblId = dict()
    
    for ensemblId in ensembl_id:
        

        ensemblIdUrl = url+ensemblId
        request = requests.get(ensemblIdUrl)
        page = request.content
        soup = BeautifulSoup(page, "html.parser")
        title = soup.find("title")
        
        geneInformation = title.string.strip()
        geneId = geneInformation.split()[0]
        
        
        desc = geneInformation.split("[")
        desc_name = desc[0].split(" ")
        desc_name_detail = ""
        
        for i in range(1, len(desc_name)):
            if desc_name_detail == "":
#                if "," in desc_name[i]:
#                    naming = desc_name[i].replace(",", "")
#                    desc_name_detail += naming
#                else:
                desc_name_detail += desc_name[i]
            else:
#                if "," in desc_name[i]:
#                    naming = desc_name[i].replace(",", "")
#                    desc_name_detail += " "+naming
#                else:
                desc_name_detail += " "+desc_name[i]
        print("ENSEMBL ID:", ensemblId, "GENE ID:", geneId, "GENE NAME:", desc_name_detail)
        
        
        dict_geneName_ensemblId[ensemblId] = [geneId, desc_name_detail]
            
    return dict_geneName_ensemblId
            
def species1Specie2GeneCorrespondance(setGeneSp1Sp2, ensemblIdSp1List, ensemblIdSp2List):
    """
    Function description:
        This function creates two dictionaries to store the pairwise 
        orthologic relationship between genes between two species.
    
    Input:
        - setGeneSp1Sp2: a file with the relationship between genes
        - ensemblIdSp1List: the species 1 list genes
        - ensemblIdSp2List: the species 2 list genes
        
    Return:
        - dict_correspondance_sp1_sp2: a dictionary with the correspondance 
                                       between genes
        - dict_correspondance_sp2_sp1: a dictionary with the correspondance 
                                       between genes
    """
    dict_correspondance_sp1_sp2 = dict()
    dict_correspondance_sp2_sp1 = dict()
    with open(setGeneSp1Sp2, 'r') as correspondance_i:
        lecture = correspondance_i.readlines()
        lecture
        for line in lecture:
            line = line.replace('\n', '').split(",")
            ensemblIdSp1 = line[0]
            ensemblIdSp2 = line[1]
            if ensemblIdSp1.startswith("ENSG"):
                if ensemblIdSp1 in ensemblIdSp1List:
                    dict_correspondance_sp1_sp2[ensemblIdSp1] = ensemblIdSp2
#                    print(ensemblIdSp1, ensemblIdSp2)
                    dict_correspondance_sp2_sp1[ensemblIdSp2] = ensemblIdSp1
            else:
#                print(ensemblIdSp1, ensemblIdSp2)
                if ensemblIdSp2 in ensemblIdSp2List:
                    dict_correspondance_sp1_sp2[ensemblIdSp1] = ensemblIdSp2
#                    print(ensemblIdSp1, ensemblIdSp2)
                    dict_correspondance_sp2_sp1[ensemblIdSp2] = ensemblIdSp1
                
    return dict_correspondance_sp1_sp2, dict_correspondance_sp2_sp1

def geneTranscriptCorrespondace(file, geneList):
    """
    Function description:
        This function creates a dictionary with a gene and all its transcripts
        associated and return a list of all transcripts studied also.
        
    Input:
        - file: a file with a gene - transcript correspondance.
        - geneList: the list of genes studied.
        
    Return:
        - dictGeneTranscript: a dictionary with all gene - transcript 
                              correspondance.
        - listTranscript: the list of all transcripts studied.
    """
    dictGeneTranscript = dict()
    listTranscript = list()
    
    with open(file, 'r') as f_i:
        lecture = f_i.readlines()
        for line in lecture:
            if not line.startswith("g"):
                gene = line.replace('\n', '').split(",")[0]
                transcript = line.replace('\n', '').split(",")[1]
                if gene in geneList:
                    listTranscript.append(transcript)
                    if not gene in dictGeneTranscript:
                        dictGeneTranscript[gene] = [transcript]
                    else:
                        dictGeneTranscript[gene].append(transcript)
                        
    return dictGeneTranscript, listTranscript

def transcriptOrthologRelationship(file, listTranscriptHS, listTranscriptMM, listTranscriptCLF):
    """
    Function description:
        This function creates dictionaries for each species to store all 
        ortholog transcript in others species.
        
    Input:
        - file: a file with the correspondance between two transcripts.
        - listTranscriptHS: list of all transcripts studied in human.
        - listTranscriptMM: list of all transcripts studied in mouse.
        - listTranscriptCLF: list of all transcripts studied in dog.
        
    Return:
        - dict_hs: a dictionary with all ortholog transcript for each human 
                   transcript
        - dict_mm: a dictionary with all ortholog transcript for each mouse 
                   transcript
        - dict_clf: a dictionary with all ortholog transcript for each dog 
                   transcript
    """
    dict_hs = dict()
    dict_mm = dict()
    dict_clf = dict()
    
    with open(file, 'r') as f_i:
        lecture = f_i.readlines()
        for line in lecture:
            line = line.replace('\n', '').split()
            tr = line[0]
            otr = line[1]
            
            if tr.startswith("ENST") or tr.startswith("CGAT"):
                if tr in listTranscriptHS:
                    if otr.startswith("ENSMUST") or otr.startswith("CGAMUST") or otr.startswith("ENSCAFT") or otr.startswith("CGACAFT"):
                        if not tr in dict_hs.keys():
                            dict_hs[tr] = [otr]
                        else:
                            dict_hs[tr].append(otr)
                    else:
                        if not tr in dict_hs.keys():
                            dict_hs[tr] = ["Na"]
                        else:
                            dict_hs[tr].append("Na")
                        
            elif tr.startswith("ENSMUST") or tr.startswith("CGAMUST"):
                if tr in listTranscriptMM:
                    if otr.startswith("ENST") or otr.startswith("CGAT") or otr.startswith("ENSCAFT") or otr.startswith("CGACAFT"):
                        if not tr in dict_mm.keys():
                            dict_mm[tr] = [otr]
                        else:
                            dict_mm[tr].append(otr)
                    else:
                        if not tr in dict_mm.keys():
                            dict_mm[tr] = ["Na"]
                        else:
                            dict_mm[tr].append("Na")
                        
            elif tr.startswith("ENSCAFT") or tr.startswith("CGACAFT"):
                if tr in listTranscriptCLF:
                    if otr.startswith("ENST") or otr.startswith("CGAT") or otr.startswith("ENSMUST") or otr.startswith("CGAMUST"):
                        if not tr in dict_clf.keys():
                            dict_clf[tr] = [otr]
                        else:
                            dict_clf[tr].append(otr)
                    else:
                        if not tr in dict_clf.keys():
                            dict_clf[tr] = ["Na"]
                        else:
                            dict_clf[tr].append("Na")
                        
    return dict_hs, dict_mm, dict_clf

def globalGeneInformation(species, species2, species3, geneList, geneListIdName, geneListIdNameSp2, geneListIdNameSp3, geneOrthologySp1_2, geneOrthologySp1_3, geneTranscriptCorrespondance, transcriptOrtholy):

    """
    """
    with open("resume_gene_information_"+species+"_"+species2+"_"+species3+".csv", 'w') as resume_file:
        resume_file.write(species+"_ensemble_id_gene"+"\t"+"gene_id"+"\t"+"gene_name"+"\t"+species2+"_ensembl_id"+"\t"+species2+"_gene_id"+"\t"+species2+"_gene_name"+"\t"+species3+"_ensembl_id"+"\t"+species3+"_gene_id"+"\t"+species3+"_gene_name"+"\t"+"transcript"+"\t"+"status"+"\t"+species2+"_ortholog_transcript"+"\t"+"status"+"\t"+species3+"_ortholog_transcript"+"\t"+"status"+"\n")
        for gene in geneList:
            gene_id = geneListIdName[gene][0]
            gene_name = geneListIdName[gene][1]
            ortholog_gene_sp2 = geneOrthologySp1_2[gene]
            ortholog_gene_sp3 = geneOrthologySp1_3[gene]
            for transcript in geneTranscriptCorrespondance[gene]:
                status = "known"
                if transcript.startswith("CGA"):
                    status = "predicted"
                
                transcriptOrtholySpecies2 = transcriptOrtholy[transcript][0]
                status_sp2 = "known"
                if transcriptOrtholySpecies2.startswith("CGA"):
                    status_sp2 = "predicted"
                
                transcriptOrtholySpecies3 = transcriptOrtholy[transcript][1]
                status_sp3 = "known"
                if transcriptOrtholySpecies3.startswith("CGA"):
                    status_sp3 = "predicted"
                
                gene_id_sp2 = geneListIdNameSp2[ortholog_gene_sp2][0]
                gene_name_sp2 = geneListIdNameSp2[ortholog_gene_sp2][1]
                
                gene_id_sp3 = geneListIdNameSp3[ortholog_gene_sp3][0]
                gene_name_sp3 = geneListIdNameSp3[ortholog_gene_sp3][1]
                
                resume_file.write(gene+"\t"+gene_id+"\t"+gene_name+"\t"+ortholog_gene_sp2+"\t"+gene_id_sp2+"\t"+gene_name_sp2+"\t"+ortholog_gene_sp3+"\t"+gene_id_sp3+"\t"+gene_name_sp3+"\t"+transcript+"\t"+status+"\t"+transcriptOrtholySpecies2+"\t"+status_sp2+"\t"+transcriptOrtholySpecies3+"\t"+status_sp3+"\n")
    
    return
                

if __name__ == "__main__":
    
    setHSGenes = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/set_genes/set_253_human.txt"#test_hs.txt"
    listGenesFromHSFile = listOfInterestedGene(setHSGenes)
    setMMGenes = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/set_genes/set_253_mouse.txt"#test_mm.txt"
    listGenesFromMMFile = listOfInterestedGene(setMMGenes)
    setCLFGenes = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/set_genes/set_253_dog.txt"#test_clf.txt"
    listGenesFromCLFFile = listOfInterestedGene(setCLFGenes)
    
    geneOrtholyHSMM = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/ortho_hs_mm_filtered2.csv"
    geneOrtholyHSCLF = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/ortho_hs_clf_filtered2.csv"
    geneOrtholyMMCLF = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/genes_one2one_hs_mm_clf_with_at_least_2_CCDS/ortho_mm_clf_filtered2.csv"
    
    geneTranscriptCorrespondanceHS = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/transcripts_hs_withprediction_final.csv"
    geneTranscriptCorrespondanceMM = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/transcripts_mm_withprediction_final.csv"
    geneTranscriptCorrespondanceCLF = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/transcripts_clf_withprediction_final.csv"
    
    pairTortho = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/paireTorthocsv"
    
    termeGO = ""
    
#    url = "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000109943"
    url = "https://www.ncbi.nlm.nih.gov/gene/?term="#ENSG00000109943"
    
    print("STEP 1: STORE GENE ID AND NAME FROM ENSEMBL ID")
    dict_geneIdNameHS = geneNameFromEnsemblId(listGenesFromHSFile, url)
    dict_geneIdNameMM = geneNameFromEnsemblId(listGenesFromMMFile, url)
    dict_geneIdNameCLF = geneNameFromEnsemblId(listGenesFromCLFFile, url)
    
    print("STEP 2: STORE ORTHOLOGY GENE CORRESPONDACE BETWEEN SPECIES")
    dict_HS_MM_gene_correspondance, dict_MM_HS_gene_correspondance = species1Specie2GeneCorrespondance(geneOrtholyHSMM, listGenesFromHSFile, listGenesFromMMFile)
    dict_HS_CLF_gene_correspondance, dict_CLF_HS_gene_correspondance = species1Specie2GeneCorrespondance(geneOrtholyHSCLF, listGenesFromHSFile, listGenesFromCLFFile)
    dict_MM_CLF_gene_correspondance, dict_CLF_MM_gene_correspondance = species1Specie2GeneCorrespondance(geneOrtholyMMCLF, listGenesFromMMFile, listGenesFromCLFFile)
    
    print("STEP 3: STORE TRANSCRIPT LIST ASSOCIATED TO EACH GENE")
    dict_gene_transcript_correspondance_HS, list_transcript_HS = geneTranscriptCorrespondace(geneTranscriptCorrespondanceHS, listGenesFromHSFile)
    dict_gene_transcript_correspondance_MM, list_transcript_MM = geneTranscriptCorrespondace(geneTranscriptCorrespondanceMM, listGenesFromMMFile)
    dict_gene_transcript_correspondance_CLF, list_transcript_CLF = geneTranscriptCorrespondace(geneTranscriptCorrespondanceCLF, listGenesFromCLFFile)
    
    print("STEP 4: STORE ALL ORTHOLOG RELATIONSHIP BETWEEN TRANSCRIPTS")
    dict_transcript_orthology_HS, dict_transcript_orthology_MM, dict_transcript_orthology_CLF = transcriptOrthologRelationship(pairTortho, list_transcript_HS, list_transcript_MM, list_transcript_CLF)
    
    print("STEP 5: CREATION OF THE COMPLETE FILE")
    globalGeneInformation("human", "mouse", "dog", listGenesFromHSFile, dict_geneIdNameHS, dict_geneIdNameMM, dict_geneIdNameCLF, dict_HS_MM_gene_correspondance, dict_HS_CLF_gene_correspondance, dict_gene_transcript_correspondance_HS, dict_transcript_orthology_HS)
#    globalGeneInformation("mouse", "human", "dog", listGenesFromMMFile, dict_geneIdNameMM, dict_MM_HS_gene_correspondance, dict_MM_CLF_gene_correspondance, dict_gene_transcript_correspondance_MM, dict_transcript_orthology_MM)
#    globalGeneInformation("dog", "human", "mouse", listGenesFromCLFFile, dict_geneIdNameCLF, dict_CLF_HS_gene_correspondance, dict_CLF_MM_gene_correspondance, dict_gene_transcript_correspondance_CLF, dict_transcript_orthology_CLF)
    
#    dict_transcriptGeneHS = 
#    dict_transcriptGeneMM = 
#    dict_transcriptGeneCLF = 
