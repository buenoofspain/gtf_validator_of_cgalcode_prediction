#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:02:23 2019

@author: Nicolas Guillaudeux
"""

import argparse

def orthologous_genes(cga_file, genes_set, cga_output_file):
    """
    Function description:
        
    """
    genes_interested = list()
    with open(genes_set, 'r') as f:
        lecture = f.readlines()
        for line in lecture:
            line = line.replace('\n', '')
            genes_interested.append(line)
    
    orthologous_genes = list()
    with open(cga_file, 'r') as f:
        lecture = f.readlines()
        for line in lecture:
            line = line.replace('\n', '')
            if line.startswith(">>> Genes"):
                line = line.replace(">>> Genes: ", "").split(" ")
                if line[0] in genes_interested:
                    orthologous_genes.append([line[0], line[1], line[2]])
    
    gene_model = dict()
    error_model = list()
    with open(cga_output_file, 'r') as cga_f:
        lecture = cga_f.readlines()
        for line in lecture:
            if line.startswith("GO"):
                line = line.replace('\n', '').split(" ")
                for genes_ortholog in orthologous_genes:
                    if line[1] in genes_ortholog:
                        if not line[1] in gene_model:
                            gene_model[line[1]] = line[3]
                            
                            
                            
                        elif gene_model[line[1]] == line[3]:
                            pass
                        else:
                                
                            error_model.append(line[1])
#                            print(line[1], gene_model[line[1]])
#                            print(line[1], line[3])
#                            print()
#                            print(genes_ortholog)
#                            print(line)
                            
    gene_model_common = dict()
    value = 0
    for gene in gene_model:
        if gene not in error_model:
            value += 1
            
    print(value)
    
    return
        
def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="Extract common gene model.")
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-g", "--genes",
                        dest = "set_genes_interested", 
                        metavar = "A set of genes.", 
                        help = "")
    
    parser.add_argument("-cga2",
                        dest = "output2_cgalcode_file",
                        metavar = "A CG-alcode output.",
                        help = "")
    
    parser.add_argument("-cga5",
                        dest = "output5_cgalcode_file",
                        metavar = "A CG-alcode output.",
                        help = "")
    
    return parser.parse_args()
    


if __name__ == "__main__":
    
    OPTIONS = parse()
    
    genes = OPTIONS.set_genes_interested
    output2 = OPTIONS.output2_cgalcode_file
    output5 = OPTIONS.output5_cgalcode_file
    
    orthologous_genes(output5, genes, output2)
