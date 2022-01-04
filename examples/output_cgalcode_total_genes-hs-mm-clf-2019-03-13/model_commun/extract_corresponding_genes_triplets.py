#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 11:12:14 2018

@author: Nicolas Guillaudeux

With Python 3.6
"""

import argparse
import re

def splitModelLetters(model):
    '''
    # Function definition:
    This function return a letters list of a model
    
    # INPUT:
    @model: the model to split
    
    # OUTPUT:
    @list_independent_letters: a list with letters of a model only 
    '''
    regex = r"([\]\[<>.]?){0,}"
    model_split = re.split(regex, model, flags=re.IGNORECASE)
    list_letters = [i for i in model_split if i] 
    # if i is not "", ex: ["A", "", "", "", "", "B", ""] from [A<.>[B]] 
    # we conserve only: ["A", "B"]
    
    list_independent_letters = list()
    pattern = r'[0-9+]'
    for letters in list_letters:
        # To separate letters collapse, ex: AB frorm [AB<.>C]
        if len(letters) > 1 and not re.search(pattern, letters):
            for letter in letters:
                list_independent_letters.append(letter)
        else:
            list_independent_letters.append(letters)
            
    return list_independent_letters, len(list_independent_letters)

def defineCorrectGenes(gene_1, gene_2, gene_3):
    '''
    # Function definition:
    This function return species associated genes.
    
    # INPUT: 
    three genes to classify
    
    # OUTPUT:
    @gene_hs, @gene_mus, @gene_caf
    '''
    if gene_1.startswith("ENSG"):
        gene_hs = gene_1
    elif gene_1.startswith("ENSCAFG"):
        gene_caf = gene_1
    elif gene_1.startswith("ENSMUSG"):
        gene_mus = gene_1
    
    if gene_2.startswith("ENSG"):
        gene_hs = gene_2
    elif gene_2.startswith("ENSCAFG"):
        gene_caf = gene_2
    elif gene_2.startswith("ENSMUSG"):
        gene_mus = gene_2
        
    if gene_3.startswith("ENSG"):
        gene_hs = gene_3
    elif gene_3.startswith("ENSCAFG"):
        gene_caf = gene_3
    elif gene_3.startswith("ENSMUSG"):
        gene_mus = gene_3
    
    return gene_hs, gene_mus, gene_caf

def extractGOWithinOutputFile(file):
    '''
    # Function definition:
    This function creates a dictionary in this format: 
        {pair_number:{gene:[models_list]}}
        where: 
            * pair_number is the pair number compare in file
            * gene is the gene name
            * models_list is a list within Gene Oracle model
            
    # INPUT:
    @file: a file output of CG-Alcode programm
    => Three dictionaries output of extractGOWithinOutputFile(file) where file is an output of CG-Alcode programm
    '''
    print(">>>> START FUNCTION TO EXTRACT GENE MODELS IN DICTIONARY")
    pair_number = ""
    pairs = 0
    gene_model_dict = dict()
    transcripts_genes = dict()
    with open(file, 'r') as filin:
    #        print("===> FILE OK")
        lines = filin.readlines()
        for line in lines:
            # To define genes pairs
            if line.startswith("OP"):
                pairs += 1
                regex = r"([ ]){0,}"
                line_OP = line.rstrip()
                line_OP = re.split(regex, line_OP, flags=re.IGNORECASE)
                
            # To search transcripts number
            elif line.startswith("TL"):
                line = line.rstrip().split(" ")
                # To define the number of genes transcripts
                transcript_number = line[4]
                # To define the number of genes transcripts known
                transcript_number = 0
                for test in line:
                    if ("ENST" in test or "ENSCAFT" in test or "ENSMUST" in test) and "CGA" not in test:
                        transcript_number += 1
                transcript_number = str(transcript_number)
                
                dict_intermediate = dict()
                # To define the list of genes name
                list_genes = list() 
                # To define genes pair number
                pair_number = str(line_OP[2])
                # To define gene name of specie 1
                gene_sp_1 = line_OP[4]
                # To define gene name of specie 2                
                gene_sp_2 = line_OP[6]
                # To add genes name in list_genes
                list_genes.append(gene_sp_1)
                list_genes.append(gene_sp_2)
                # To create an intermediate dictionary
                dict_intermediate = dict([(gene, list()) for gene in list_genes])
                # To prepare the dictionary
                for gene in list_genes:
                    try:
                        gene_model_dict[pair_number] = dict_intermediate
                    except KeyError:
                        pass
                if line[2] not in transcripts_genes:
                    transcripts_genes[line[2]] = transcript_number
            
            # To define model by gene
            elif line.startswith("GO"):
                line = line.rstrip().split(" ")
               # print(line[0], line[1], line[3])
                model = line[3]
                try:
                    gene_model_dict[pair_number][line[1]].append(model)
                except KeyError:
                    pass
                
    print "number of pairs genes:", pairs
    return gene_model_dict, transcripts_genes

def dictCompareModels(dict1, dict2, dict3):
    """
    # Function definition:
    This function creates two objects : 
        - a list of pair model different by one letter in this format: 
            [[pair_value,gene_sp_1,transcripts_number_sp_1,gene_sp_2,transcripts_number_sp_2,gene_sp_3,transcripts_number_sp_3,model_1,model_2,model_3],...]

        - a dictonary of genes with the same Gene Oracle model in this format:
            {(pair_value, gene_sp_1, transcript_number_sp1, gene_sp2..., transcript_number_sp3):[common_model],...}
            
                * pair_value is the pair number
                * gene_sp_1/2/3 are genes name with common model
                * transcripts_number_sp_1/2/3 are the number of transcripts associated to the gene
                * model_1/2/3 are illustrated model
                * common_model is a list within Gene Oracle model
            
    # INPUT:
    @dict1, dict_2, dict_3: dictionary with pair_number, gene and model associated.
    
    # OUTPUT:
    @genes_models_dict: a dictionay within commons models.
    @genes_models_list: a list with models different by one letter between genes pairs.
    """
    genes_models_dict = dict() # To initialize intermediaite dictionary
    genes_models_list = list() # To initialize intermediaite dictionary
    #genes_list = list()
    count_1 = 0 # To initialize a counter to know how many comparisons are kept
    gene_C = "", "", "" # To initialize intermediate gene
    gene_hs, gene_mus, gene_caf = "", "", "" # To initialize genes to the output
    
    dict_1, tr_1 = dict1
    dict_2, tr_2 = dict2
    dict_3, tr_3 = dict3
        
    # To treat first dictionary
    for pair_value_1, gene_dict_1 in dict_1.items(): # To treat pair number and genes dictionary in two variable
        #model_A = "" # To initialize the model to compare with others
        gene_1, gene_2 = gene_dict_1.keys() # To store the two genes
        # To verify that each gene has only one model
        if len(gene_dict_1[gene_1]) == 1 and len(gene_dict_1[gene_2]) == 1: 
            # To standardize the model names (remove the ".")
            for model_1 in gene_dict_1[gene_1]:
                if model_1.startswith("."):
                    model_1 = model_1[1:]
            for model_2 in gene_dict_1[gene_2]:
                if model_2.startswith("."):
                    model_2 = model_2[1:]
            # To compare gene models between the two genes and to store the common model if it exists
            if model_1 == model_2 :
                
                # To treat second dictionary
                for pair_value_2, gene_dict_2 in dict_2.items():
                    gene_3, gene_4 = gene_dict_2.keys()
                    # To verify if pair numbers are same and to verify that each gene has only one model
                    if pair_value_1 == pair_value_2 and len(gene_dict_2[gene_3]) == 1 and len(gene_dict_2[gene_4]) == 1 : 
                        # To standardize the model names (remove the ".")
                        for model_3 in gene_dict_2[gene_3]:
                            if model_3.startswith("."):
                                model_3 = model_3[1:]
                        for model_4 in gene_dict_2[gene_4]:
                            if model_4.startswith("."):
                                model_4 = model_4[1:]
                        # To validate if models are different
                        if model_3 == model_4 and model_3 != model_1 :
                            letters_list_model_1, len_model_1 = splitModelLetters(model_1)
                            letters_list_model_3, len_model_3 = splitModelLetters(model_3)
                            if abs(len_model_1-len_model_3) == 1:
                                if gene_3 != gene_1 and gene_3 != gene_2:
                                    gene_C = gene_3
                                elif gene_4 != gene_1 and gene_4 != gene_2:
                                    gene_C = gene_4
                                if len_model_1 < len_model_3:
                                    min_letter_difference, max_letter_difference = len_model_1, len_model_3
                                else:
                                    min_letter_difference, max_letter_difference = len_model_3, len_model_1
            
                                # To treat third dictionary
                                for pair_value_3, gene_dict_3 in dict_3.items():
                                    gene_5, gene_6 = gene_dict_3.keys()
                                    # To verify if pair numbers are same and to verify that each gene has only one model
                                    if pair_value_1 == pair_value_3 and len(gene_dict_3[gene_5]) == 1 and len(gene_dict_3[gene_6]) == 1:
                                        # To standardize the model names (remove the ".")
                                        for model_5 in gene_dict_3[gene_5]:
                                            if model_5.startswith("."):
                                                model_5 = model_5[1:]
                                        for model_6 in gene_dict_3[gene_6]:
                                            if model_6.startswith("."):
                                                model_6 = model_6[1:]
                                        # To validate if the models are different
                                        if model_5 == model_6 and ((model_5 != model_1 and model_5 == model_3) or (model_5 != model_3 and model_5 == model_1)):
                                            letters_list_model_5, len_model_5 = splitModelLetters(model_5)
                                            count_1 += 1
                                            
                                            # To set the correct order for the output file
                                            gene_hs, gene_mus, gene_caf = defineCorrectGenes(gene_1, gene_2, gene_C)

                                           #print pair_value_1, "\n",gene_hs, gene_mus, gene_caf,"\n", model_1, "\n", model_2, "\n", model_3, "\n", model_4, "\n", model_5, "\n", model_6, "\n"
                                           # To write in a list the exit line
                                            if gene_hs in tr_1 and gene_mus in tr_1 and gene_hs in tr_2:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_1[gene_hs]+","+gene_mus+","+tr_1[gene_mus]+","+gene_caf+","+tr_2[gene_caf]+","+str(len_model_1)+","+str(len_model_3)+","+str(len_model_5)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                            elif gene_hs in tr_1 and gene_mus in tr_1 and gene_hs in tr_3:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_1[gene_hs]+","+gene_mus+","+tr_1[gene_mus]+","+gene_caf+","+tr_2[gene_caf]+","+str(len_model_1)+","+str(len_model_5)+","+str(len_model_3)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                            elif gene_hs in tr_2 and gene_mus in tr_2 and gene_hs in tr_1:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_2[gene_hs]+","+gene_mus+","+tr_2[gene_mus]+","+gene_caf+","+tr_1[gene_caf]+","+str(len_model_3)+","+str(len_model_1)+","+str(len_model_5)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                            elif gene_hs in tr_2 and gene_mus in tr_2 and gene_hs in tr_3:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_2[gene_hs]+","+gene_mus+","+tr_2[gene_mus]+","+gene_caf+","+tr_1[gene_caf]+","+str(len_model_3)+","+str(len_model_5)+","+str(len_model_1)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                            elif gene_hs in tr_3 and gene_mus in tr_3 and gene_hs in tr_1:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_3[gene_hs]+","+gene_mus+","+tr_3[gene_mus]+","+gene_caf+","+tr_1[gene_caf]+","+str(len_model_5)+","+str(len_model_1)+","+str(len_model_3)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                            elif gene_hs in tr_3 and gene_mus in tr_3 and gene_hs in tr_2:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_3[gene_hs]+","+gene_mus+","+tr_3[gene_mus]+","+gene_caf+","+tr_1[gene_caf]+","+str(len_model_5)+","+str(len_model_3)+","+str(len_model_1)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                                
                                            # To store the list
                                            genes_models_list.append(genes_list)                                            
                                            
                        # To validate if the models are identical to the common model             
                        elif model_3 == model_4 and model_3 == model_1 :
                            letters_list_model_3, len_model_3 = splitModelLetters(model_3)
                            if gene_3 != gene_1 and gene_3 != gene_2:
                                gene_C = gene_3
                            elif gene_4 != gene_1 and gene_4 != gene_2:
                                gene_C = gene_4
                                
                            # To treat third dictionary
                            for pair_value_3, gene_dict_3 in dict_3.items():
                                gene_5, gene_6 = gene_dict_3.keys()
                                # To verify if pair numbers are same and to verify that each gene has only one model
                                if pair_value_1 == pair_value_3 and len(gene_dict_3[gene_5]) == 1 and len(gene_dict_3[gene_6]) == 1:
                                    # To standardize the model names (remove the ".")
                                    for model_5 in gene_dict_3[gene_5]:
                                        if model_5.startswith("."):
                                            model_5 = model_5[1:]
                                    for model_6 in gene_dict_3[gene_6]:
                                        if model_6.startswith("."):
                                            model_6 = model_6[1:]
                                    # To validate if the models are different
                                    if model_5 == model_6 and model_5 != model_1:
                                        letters_list_model_1, len_model_1 = splitModelLetters(model_1)
                                        letters_list_model_5, len_model_5 = splitModelLetters(model_5)
                                        if abs(len_model_1-len_model_5) == 1:
                                            count_1 += 1
                                            if len_model_1 < len_model_5:
                                                min_letter_difference, max_letter_difference = len_model_1, len_model_5
                                            else:
                                                min_letter_difference, max_letter_difference = len_model_5, len_model_1
                                            
                                            # To set the correct order for the output file
                                            gene_hs, gene_mus, gene_caf = defineCorrectGenes(gene_1, gene_2, gene_C)
                                            
                                           # print pair_value_1, "\n",gene_hs, gene_mus, gene_caf,"\n", model_1, "\n", model_2, "\n", model_3, "\n", model_4, "\n", model_5, "\n", model_6, "\n"
                                           # To write in a list the exit line
                                            if gene_hs in tr_1 and gene_mus in tr_1 and gene_hs in tr_2:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_1[gene_hs]+","+gene_mus+","+tr_1[gene_mus]+","+gene_caf+","+tr_2[gene_caf]+","+str(len_model_1)+","+str(len_model_3)+","+str(len_model_5)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                            elif gene_hs in tr_1 and gene_mus in tr_1 and gene_hs in tr_3:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_1[gene_hs]+","+gene_mus+","+tr_1[gene_mus]+","+gene_caf+","+tr_2[gene_caf]+","+str(len_model_1)+","+str(len_model_5)+","+str(len_model_3)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                            elif gene_hs in tr_2 and gene_mus in tr_2 and gene_hs in tr_1:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_2[gene_hs]+","+gene_mus+","+tr_2[gene_mus]+","+gene_caf+","+tr_1[gene_caf]+","+str(len_model_3)+","+str(len_model_1)+","+str(len_model_5)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                            elif gene_hs in tr_2 and gene_mus in tr_2 and gene_hs in tr_3:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_2[gene_hs]+","+gene_mus+","+tr_2[gene_mus]+","+gene_caf+","+tr_1[gene_caf]+","+str(len_model_3)+","+str(len_model_5)+","+str(len_model_1)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                            elif gene_hs in tr_3 and gene_mus in tr_3 and gene_hs in tr_1:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_3[gene_hs]+","+gene_mus+","+tr_3[gene_mus]+","+gene_caf+","+tr_1[gene_caf]+","+str(len_model_5)+","+str(len_model_1)+","+str(len_model_3)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                            elif gene_hs in tr_3 and gene_mus in tr_3 and gene_hs in tr_2:
                                                genes_list = str(pair_value_1+","+gene_hs+","+tr_3[gene_hs]+","+gene_mus+","+tr_3[gene_mus]+","+gene_caf+","+tr_1[gene_caf]+","+str(len_model_5)+","+str(len_model_3)+","+str(len_model_1)+","+str(min_letter_difference)+" / "+str(max_letter_difference))
                                               
                                            # To store the list
                                            genes_models_list.append(genes_list)
                                            
                                    # To validate if the models are identical to the common model
                                    elif model_5 == model_6 and model_6 == model_1:
                                        count_1 += 1
                                                  
                                        # To set the correct order for the output file
                                        gene_hs, gene_mus, gene_caf = defineCorrectGenes(gene_1, gene_2, gene_C)
                                        
                                       # print pair_value_1, "Common model", "\n",gene_hs, gene_mus, gene_caf,"\n", model_1, "\n", model_2, "\n", model_3, "\n", model_4, "\n", model_5, "\n", model_6, "\n"
                                       # To write the output in a dictionary
                                        if gene_hs in tr_1 and gene_mus in tr_1 :
                                            genes_commons_list = str(pair_value_1+","+gene_hs+","+tr_1[gene_hs]+","+gene_mus+","+tr_1[gene_mus]+","+gene_caf+","+tr_2[gene_caf])
                                        elif gene_hs in tr_2 and gene_mus in tr_2 :
                                            genes_commons_list = str(pair_value_1+","+gene_hs+","+tr_2[gene_hs]+","+gene_mus+","+tr_2[gene_mus]+","+gene_caf+","+tr_1[gene_caf])
                                        elif gene_hs in tr_3 and gene_mus in tr_3 :
                                            genes_commons_list = str(pair_value_1+","+gene_hs+","+tr_3[gene_hs]+","+gene_mus+","+tr_3[gene_mus]+","+gene_caf+","+tr_1[gene_caf])
                                            
                                        try:
                                            genes_models_dict[genes_commons_list] = (gene_dict_1[gene_1])
                                        except KeyError:
                                            pass
                                            
    print "Number of comparisons conserved:", count_1
    print "Number of genes triplets conserved:", len(genes_models_dict)
    print "Number of genes triplets conserved:", len(genes_models_list)
    return genes_models_dict, genes_models_list

if __name__ == "__main__":
    
    # Parser
    parser = argparse.ArgumentParser(prog = "extract_corresponding_triplets_with_same_letter_number.py")
    
    # Output required parameters
    parser.add_argument("-a", "--output_1", dest = "output_1", metavar = "CG_ALCODE_OUTPUT_1", help = "Output of CG-alcode programm", required = True)
    parser.add_argument("-b", "--output_2", dest = "output_2", metavar = "CG_ALCODE_OUTPUT_2", help = "Output of CG-alcode programm", required = True)
    parser.add_argument("-c", "--output_3", dest = "output_3", metavar = "CG_ALCODE_OUTPUT_3", help = "Output of CG-alcode programm", required = True)

    # Args definition
    args = parser.parse_args()
    
    output_1 = args.output_1#"output_hs_clf.txt"  #args.output_1
    output_2 = args.output_2#"output_hs_mm.txt"  #args.output_2
    output_3 = args.output_3#"output_mm_clf.txt"  #args.output_3
    
 #   output_1 = "output_hs_clf.txt"  #args.output_1
 #   output_2 = "output_mm_clf.txt"  #args.output_2
 #   output_3 = "output_hs_mm.txt"  #args.output_3
    
 #   output_1 = "output_hs_mm.txt"  #args.output_1
 #   output_2 = "output_hs_clf.txt"  #args.output_2
 #   output_3 = "output_mm_clf.txt"  #args.output_3

 #   output_1 = "output_hs_mm.txt"  #args.output_1
 #   output_2 = "output_mm_clf.txt"  #args.output_2
 #   output_3 = "output_hs_clf.txt"  #args.output_3
    
 #   output_1 = "output_mm_clf.txt"  #args.output_1
 #   output_2 = "output_hs_clf.txt"  #args.output_2
 #   output_3 = "output_hs_mm.txt"  #args.output_3
    
 #   output_1 = "output_mm_clf.txt"  #args.output_1
 #   output_2 = "output_hs_mm.txt"  #args.output_2
 #   output_3 = "output_hs_clf.txt"  #args.output_3
    
    dict_1 = extractGOWithinOutputFile(output_1)
    print("===> FIRST DICTIONARY: " + str(len(dict_1) > 1))
    print("********************************************")
    print ("")
    
    dict_2 = extractGOWithinOutputFile(output_2)
    print("===> SECOND DICTIONARY: " + str(len(dict_2) > 1))
    print("********************************************")
    print ("")
    
    dict_3 = extractGOWithinOutputFile(output_3)
    print("===> THIRD DICTIONARY: " + str(len(dict_3) > 1))
    print("********************************************")
    print("")
    
    print(">>>> START COMPARISON")
    genes_models_dict, interest_genes_models = dictCompareModels(dict_1, dict_2, dict_3)
    print("********************************************")
    print("")
    
    print(">>>> START OUTPUT FOR COMMON MODEL")
    with open('output_common_comparison.csv', 'w') as filout:
        title="GENE NUMBER,HUMAN GENE,NUMBER OF HUMAN TRANSCRIPTS,MOUSE GENE,NUMBER OF MOUSE TRANSCRIPTS,DOG GENE,NUMBER OF DOG TRANSCRIPTS,COMMON MODEL\n"
        filout.write(title)
        for genes, model_list in genes_models_dict.items():
            for model in model_list:
                if model.startswith("."):
                    model = model[1:]
                line = genes+','+model
                filout.write(line+'\n')
        print("===> OUTPUT : OK")
        print("********************************************")
        print ("")
    
    with open('output_one_letter_difference_comparison.csv', 'w') as filout:
        title="GENE NUMBER,HUMAN GENE,NUMBER OF HUMAN TRANSCRIPTS,MOUSE GENE,NUMBER OF MOUSE TRANSCRIPTS,DOG GENE,NUMBER OF DOG TRANSCRIPTS,NUMBER OF BLOCKS IN HUMAN-MOUSE,NUMBER OF BLOCKS IN HUMAN-DOG,NUMBER OF BLOCKS IN MOUSE-DOG,NUMBER OF MINIMUM / MAXIMUM BLOCKS\n"
        filout.write(title)
        for model_list in interest_genes_models:
            line = model_list
            filout.write(line+'\n')
        print("===> OUTPUT : OK")
        print("********************************************")
