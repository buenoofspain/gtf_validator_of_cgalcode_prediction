#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:53:31 2019

@author: niguilla
"""

import argparse

def gtf2bed(file):
    '''
    Function description:
        This function convert a gtf file into bed file with 12 columns (bed12)
        with utr if exists (.tmp) and without utr (.bed).
    
    Input:
        - file: a gtf file to convert in bed format.
        - utr: if you want to del utr information choose True.
    '''    

    with open(file, 'r') as file_tmp:
        with open(file[:-4]+"_without_utr.bed", 'w') as file_bed:
            for line in file_tmp.read().splitlines():
                line = line.split('\t')
#                        print(line)
                a = line[10].split(',')
#                        print(a)
                b = line[11].split(',')
                A = ""
                B = ""
                size = 0
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
                            

def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="")
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-i", "--input",
                        dest = "input", 
                        metavar = "INPUT FILE WITH UTR", 
                        help = "Give a file with utr information to obtain a file without utr.")
    
    
    return parser.parse_args()                     

if __name__ == "__main__":
    
    print("########################################")
    print("########## PROGRAMM BEGINNING ##########")
    print("")
    
    OPTIONS = parse()
    
    file = OPTIONS.input
#    file = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/examples/reference_gtf/mouse/xbseq/mm10_copie.bed"
    
    gtf2bed(file)