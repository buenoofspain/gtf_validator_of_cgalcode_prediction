#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:16:09 2019

@author: niguilla
"""

from matplotlib_venn import venn3, venn2, venn2_circles, venn3_circles
import matplotlib.pyplot as plt
#from sympy import FiniteSet # to use sets
import numpy as np

plt.figure(figsize=(4.2,4)) # largeur puis hauteur

## TO TEST
#figure, axes = plt.subplots(2, 2)
#venn2(subsets={'10': 1, '01': 11, '11': 1}, set_labels = ('A', 'B'), ax=axes[0][0])
#venn2_circles((1, 2, 3), ax=axes[0][1])
#venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('A', 'B', 'C'), ax=axes[1][0])
#venn3_circles({'001': 10, '100': 20, '010': 21, '110': 13, '011': 14}, ax=axes[1][1])
#plt.show()


    
if __name__ == "__main__":
    path = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/set_genes/"
    clique_file1 = path+"set_255_signal.txt"
    clique_file2 = path+"c135_human.txt"
    
    c_list = []
    c_list_str = ""
    db_list = []
    db_list_str = ""
    s_list = []
    s_list_str = ""
    
    with open(clique_file1, 'r') as c_f:
        c_f_lecture = c_f.readlines()
        for gene_c in c_f_lecture:
            gene_c = gene_c.replace('\n', '')
            c_list.append(gene_c)
            
    with open(clique_file2, 'r') as db_f:
        db_f_lecture = db_f.readlines()
        for gene_db in db_f_lecture:
            gene_db = gene_db.replace('\n', '')
            db_list.append(gene_db)
    
    for i in c_list:
        if c_list_str == "":
            c_list_str = str(i)+","
        else:
            c_list_str += str(i)+","
    
    for j in db_list:
        if db_list_str == "":
            db_list_str = str(j)+","
        else:
            db_list_str += str(j)+","
            
    set1 = set([i for i in c_list])
    set2 = set([j for j in db_list])
    
    v = venn2([set1, set2], set_labels = ("", ""))
#    v.get_label_by_id('A').set_text('Set 222')
#    v.get_label_by_id('B').set_text('Set 387')
    #v.get_label_by_id('A').set_text('$x^2$') # Those are set labels
    v.get_label_by_id('A').set_fontsize(18)
    v.get_label_by_id('B').set_fontsize(18)
    
    
    
    v.get_patch_by_id('01').set_alpha(0.4)
    v.get_patch_by_id('01').set_color('green')
    v.get_patch_by_id('11').set_alpha(0.4)
    v.get_patch_by_id('11').set_color('orange')
    
    # To define size
    label = v.get_label_by_id('11')          # Those are subset labels (i.e. numbers)
    label.set_fontsize(18) 
    #label.set_family('serif')
    label.set_x(label.get_position()[0])
    
    label = v.get_label_by_id('01')
    label.set_fontsize(18)
    label.set_x(label.get_position()[0])
    
    label = v.get_label_by_id('10')
    label.set_fontsize(18)
    label.set_x(label.get_position()[0])
    
    label = v.get_label_by_id('A')
    label.set_fontsize(18)
    label.set_x(label.get_position()[0]-0.2)
    
    label = v.get_label_by_id('B')
    label.set_fontsize(18)
    label.set_x(label.get_position()[0])
    
    
    
#    plt.annotate('Interested set', xy=v.get_label_by_id('11').get_position() - np.array([0, 0.05]), xytext=(-10,-30),
#ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
#arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
#    
#    plt.annotate('Error set', xy=v.get_label_by_id('001').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
#ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
#arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
#    
#    plt.annotate('Error set', xy=v.get_label_by_id('011').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
#ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
#arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    
#    plt.title('Venn Diagram').set_fontsize(22)
    plt.savefig("venn_255vs135.svg") #To save the graph drawing
    plt.show()
    
