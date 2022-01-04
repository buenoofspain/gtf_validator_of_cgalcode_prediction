#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 11:33:46 2020

@author: niguilla
"""

import argparse

def species_gene_list(file):
    """
    Function description:
        This function extract the gene list in a list object.
        
    Input:
        file: a file with the list of interestes genes.
        
    Return: 
        list_gene: the list of interested genes.
    """
    list_gene = list()
    
    with open(file, 'r') as input_f:
        lecture = input_f.readlines()
        for line in lecture:
            line = line.replace('\n', '')
            list_gene.append(line)
    
    return list_gene

def extract_Found_CorrectORF(cga_file, list_species1, list_species2, list_species3):
    """
    Function description:
        This function extract the list of predicted transcripts already found 
        and correct-orf in the CG-Alcode output"
        
    Input:
        - cga_file: the CG-Alcode output
        - list_speciesX: the gene list for each species.
    """
    gene_tr_found = dict()
    
    with open(cga_file, 'r') as input_f:
        lecture = input_f.readlines()
        for line in lecture:
            if line.startswith("CO "):
                line_split = line.replace('\n', '').split(" ")
                try:
                    source_gene = line_split[1]
                    source_tr = line_split[3]
                
                    target_gene = line_split[2]
                    target_tr = line_split[4]
                
                    known_status = line_split[5]
                    orf_status = line_split[6]
                    
                    if known_status == "Found" and orf_status == "Correct-ORF":
                        if (source_gene in list_species1 or source_gene in list_species2 or source_gene in list_species3) and (target_gene in list_species1 or target_gene in list_species2 or target_gene in list_species3) and not source_tr.startswith("CGA") and not target_tr.startswith("CGA"):
                            try:
                                if not target_tr in gene_tr_found[target_gene]:
                                    gene_tr_found[target_gene].append(target_tr)
                            except KeyError:
                                gene_tr_found[target_gene] = [target_tr]
                
                except IndexError:
                    pass

    return gene_tr_found

def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="A SCRIPT TO OBTAIN THE NUMBER OF PREDICTED TRANSCRIPTS ALREADY FOUND IN ANOTHER SPECIES.")
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-s1", "--set1",
                        dest = "set_genes_hs", 
#                        default = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/set_genes/set_253_human.txt",
                        metavar = "GENE LIST OF HUMAN SPECIES FILE", 
                        help = "Give the file with the set of genes with an interest.",
                        required = True)
    
    parser.add_argument("-s2", "--set2",
                        dest = "set_genes_mm", 
#                        default = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/set_genes/set_253_mouse.txt",
                        metavar = "GENE LIST OF MOUSE SPECIES FILE", 
                        help = "Give the file with the set of genes with an interest.",
                        required = True)
    
    parser.add_argument("-s3", "--set3",
                        dest = "set_genes_clf", 
#                        default = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/set_genes/set_253_dog.txt",
                        metavar = "GENE LIST OF DOG SPECIES FILE", 
                        help = "Give the file with the set of genes with an interest.",
                        required = True)
    
    
    return parser.parse_args()  

if __name__ == "__main__":
    
#    human_file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/set_genes/set_253_human.txt"
#    mouse_file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/set_genes/set_253_mouse.txt"
#    dog_file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/set_genes/set_253_dog.txt"
    
    OPTIONS = parse()
    
    human_file = OPTIONS.set_genes_hs
    mouse_file = OPTIONS.set_genes_mm
    dog_file = OPTIONS.set_genes_clf

    human_gene_list = species_gene_list(human_file)
    mouse_gene_list = species_gene_list(mouse_file)
    dog_gene_list = species_gene_list(dog_file)
    
    cgalcode_file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/output_cgalcode_total_genes-hs-mm-clf-2019-03-13/output_run2.txt"
    
    tr_predicted_alrealdy_known = extract_Found_CorrectORF(cgalcode_file, dog_gene_list, human_gene_list, mouse_gene_list)
    
    tr_already_known_human = list()
    tr_already_known_mouse = list()
    tr_already_known_dog = list()
    
    for gene in tr_predicted_alrealdy_known:
        for tr_predicted in tr_predicted_alrealdy_known[gene]:
            if gene.startswith("ENSG") and not tr_predicted in tr_already_known_human:
                tr_already_known_human.append(tr_predicted)
                #print("HUMAN", gene, tr_predicted)
            elif gene.startswith("ENSMUSG") and not tr_predicted in tr_already_known_mouse:
                tr_already_known_mouse.append(tr_predicted)
#                print("MOUSE", gene, tr_predicted)
            elif gene.startswith("ENSCAFG") and not tr_predicted in tr_already_known_dog:
                tr_already_known_dog.append(tr_predicted)
#                print("DOG", gene, tr_predicted)
    
    print("Number of predicted transcripts already known in human:", len(tr_already_known_human))
    print("Number of predicted transcripts already known in mouse:", len(tr_already_known_mouse))
    print("Number of predicted transcripts already known in dog:", len(tr_already_known_dog))