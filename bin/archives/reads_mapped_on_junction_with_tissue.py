#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:30:53 2019

@author: Aurélie Nicolas
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from math import *


    
    
if __name__ == "__main__":
    file_i = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/validation_part2/dog_junction_exons/junction_validation_dog_2019-09-27/analysed/nb_reads_mapped_tissue_9.csv"
    
    file = pd.read_csv(file_i, sep='\t')
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
    #if len(file.values) > 50:
    val_max=100
    int_jonc = ceil(len(file.values)/val_max)
    
    for i in range(int_jonc):
        file_temp = pd.DataFrame(file.iloc[i*val_max: min(len(file.values), val_max*(i+1))])
        
        
    # Création de figure pour max 100 jonction
        datafile = file_temp
        
        plt.figure(figsize=(13,5))
        # Préparation des données pour les plots
        ind = np.arange(len(datafile.values))
      
        
        liste_plot=[]
        bottom_data = [0] * len(list(datafile.values))
        # Création des plots
        for data in tissue:
            plot = plt.bar(ind, datafile[data], bottom= bottom_data)
            liste_plot.append(plot[0])
            bottom_data = [datafile[data][int_data] + bottom_data[int_data] for int_data in range(len(list(datafile[data])))]
       
        # Annotation du graphique
        for int_value in range(len(list(datafile.total))):
            plt.annotate(file.total[int_value], (int_value - 0.5 , datafile.total[int_value]+max(datafile.total)/100*2), fontsize=6)
        
        plt.legend(liste_plot, tissue, loc=9, bbox_to_anchor=(0.5, 1.00), fontsize=8, ncol=5, title="Tissue")
        int_legend = ceil(len(tissue)/8)+2
        plt.ylim(0,max(datafile.total) + max(datafile.total)/100 * 10 * int_legend)
        plt.title("Reads aligned on the specific exon couples applied to the predicted transcripts")#Jonctions en fonction du nombre de reads alignés sur les couples d'exons spécifiques au jeu de données des 111 clf en fonction des tissus.")
        plt.xticks(ind, datafile['Junctions_light'], rotation=90, fontsize=6)
        plt.ylabel("Number of reads aligned")
        plt.xlabel("Specific exon couples")
        
        plt.subplots_adjust(left = 0.06, bottom = 0.18, right = 0.98, top = 0.90)
        print('figure_'+str(int_jonc)+'.pdf')
        plt.savefig('figure_'+str(int_jonc)+'.png')
        plt.close()
    
    
    
    