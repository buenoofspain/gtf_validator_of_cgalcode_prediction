#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:30:53 2019

@author: Aurélie Nicolas
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from math import *
import itertools

def create_graph(file, number_of_junctions):
    """
    Function description: 
        This fonction create a diagramm to analyse specific exon couples.
    
    Input: 
        file: a file with th enumber of reads aligned with a tab separator
    
    Output: a pdf file with the diagramm.
    """
    file = pd.read_csv(file, sep='\t')
    tissue = list(file.columns)
    if "Junctions" in tissue:
        tissue.remove("Junctions")
    
    file['total'] = file[tissue].sum(axis=1)
    name_junction = file.Junctions.values
    liste_junction=[]
    for junction in name_junction :
        liste_junction.append(junction[-9:])
    file['Junctions_light'] = liste_junction
    
    


    # Verification que les données par figure font 100 jonctions
    val_max=int(number_of_junctions) # Nombre max de jonction par figure
    int_jonc = ceil(len(file.values)/val_max) # Nombre de figure qui sera réalisé
    
    for int_fig in range(int_jonc):
        datafile = pd.DataFrame(file.iloc[int_fig * val_max: min(len(file.values), val_max * (int_fig + 1))])
        # Re indexage du tableau
        datafile.index = range(len(datafile.index))
        # Création de la figure
        plt.figure(figsize=(20,5)) #20,5
        # Préparation des données pour les plots
        ind = np.arange(len(datafile.values))
        
        liste_plot=[]
        bottom_data = [0] * len(list(datafile.values))
        
        liste_couleur = ['orangered', 'deepskyblue','lime', 'deeppink', 'gold', 'Brown', 'silver', 'purple']
                #'skyblue','red','green', 'pink', 'gray', 'Indigo', 'Navy', 'Lime', 'Olive', 'Teal', 'aqua']
        
        nb_color = len(liste_couleur) #Nombre d'éléments dans la liste des couleurs
        liste_hatch = ['-', '.', '/', 'O', '//', 'o']
        
        #[liste_hatch.extend(list(itertools.repeat(elt, nb_color))) for elt in ['-', '.', '/', 'O', '//', 'o']]
        
        
        
        
        
       # + list(itertools.repeat('--', nb_color)) + list(itertools.repeat(':', nb_color))
        
        
        #linescycler = itertools.cycle(liste_lines)
        hatchcycler = itertools.cycle(liste_hatch)
        couleurcycler = itertools.cycle(liste_couleur)
        
#        print(len(liste_couleur))
#        print(len(liste_hatch))

        
        # Création des plots
        for data in tissue:
            #print(next(linescycler), next(markercycler), next(couleurcycler))
#            plot = plt.bar(ind, datafile[data], bottom= bottom_data)
            plot = plt.bar(ind, datafile[data], bottom= bottom_data, color = next(couleurcycler), hatch = next(hatchcycler), linewidth = 3)
            liste_plot.append(plot[0])
            bottom_data = [datafile[data][int_data] + bottom_data[int_data] for int_data in range(len(list(datafile[data])))]
            
        # Annotation du graphique
        for int_value in range(len(list(datafile.total))):
            plt.annotate(datafile.total[int_value], (int_value - min(len(datafile)/200, 0.5) , datafile.total[int_value]+max(datafile.total)/100*2), fontsize=5)
        
        plt.legend(liste_plot, tissue, loc=9, bbox_to_anchor=(0.5, 1.00), fontsize=8, ncol=8, title="Tissue") #ncol=8
        int_legend = ceil(len(tissue)/8)+2
        plt.ylim(0,max(datafile.total) + max(datafile.total)/100 * 10 * int_legend)
        plt.title("Reads aligned on the specific exon couples applied to the predicted transcripts")#Jonctions en fonction du nombre de reads alignés sur les couples d'exons spécifiques au jeu de données des 111 clf en fonction des tissus.")
        plt.xticks(ind, datafile['Junctions_light'], rotation=90, fontsize=6)
        plt.ylabel("Number of reads aligned")
        plt.xlabel("Specific exon couples")
        
        plt.subplots_adjust(left = 0.05, bottom = 0.18, right = 0.92, top = 0.90) #0.08, bottom = 0.25, right = 0.92, top = 0.90)
        plt.savefig('figure_'+str(int_fig)+'.png')
        plt.close("all")
    
    
def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="A SCRIPT TO OBTAIN SPECIFIC EXON JUNCTION NOT KNOWN IN ENSEMBL.")
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-f", "--file",
                        dest = "nb_reads_mapped", 
                        metavar = "NUMBER OF READS MAPPED IN EXON-EXON JUNCTIONS FILE", 
                        help = "Give the file with the number of reads aligned in exon-exon junctions.",
                        required = False)
    
    parser.add_argument("-n", "--nb",
                        dest = "nb_reads_by_figure", 
                        metavar = "NUMBER OF JUNCTIONS WITHIN EACH FIGURE", 
                        help = "Give a number of juncions for each figure.",
                        required = False)
    
    
    return parser.parse_args()  
    
    
    
if __name__ == "__main__":
    OPTIONS = parse()
    
    nb_read_file = OPTIONS.nb_reads_mapped
    nb_junction_total = OPTIONS.nb_reads_by_figure
    
#    file_i = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/junction_validation_dog_2019-09-27/analysed/nb_reads_mapped_tissue.csv"
    create_graph(nb_read_file, nb_junction_total)
    
    
    