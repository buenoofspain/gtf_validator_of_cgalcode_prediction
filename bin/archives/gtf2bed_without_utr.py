#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 16:00:01 2019

@author: niguilla
"""

import argparse
import subprocess

def gtf2bed(file, gtf2bedscriptconvertor, utr=False):
    '''
    Function description:
        This function convert a gtf file into bed file with 12 columns (bed12)
        with utr if exists (.tmp) and without utr (.bed).
    
    Input:
        - file: a gtf file to convert in bed format.
        - utr: if you want to del utr information choose True.
    '''    
    path_to_script = "./"+gtf2bedscriptconvertor
    
    if utr == True:
        out = file+".bed"
        cmd = "{0} {1}  >> {2}".format(path_to_script, file, out)
        subprocess.call(cmd, shell=True)
    
    if utr == False:
        out = file+".tmp"
        cmd = "{0} {1}  >> {2}".format(path_to_script, file, out)
        subprocess.call(cmd, shell=True)
        with open(out, 'r') as file_tmp:
            with open(file+".bed", 'w') as file_bed:
                for line in file_tmp.read().splitlines():
                    line = line.split('\t')
                    a = line[10].split(',')
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
                                
def parse():
    """
    Method to create a parser for command-line options
    """
    parser = argparse.ArgumentParser(description="")
    # Create command-line parser for all options and arguments to give
    parser.add_argument("-gtf2bed", "--gtf2bedconvertion",
                        default = "../../dependencies/gtf2bed.pl", 
                        dest = "gtfbedscript", 
                        metavar = "GTF2BED SCRIPT CONVERTOR", 
                        help = "Enter the path repository where are located genes and transcripts files used by CG-Alcode programm.")
    
    parser.add_argument("-ref", "--reference",
                        default = "../../examples/reference_gtf/human/ensembl/Ensembl98/Homo_sapiens.GRCh38.98.gtf",
                        # A voir : default = "/home/niguilla/Documents/these_nguillaudeux/Data/gtf_hs/new/hg38.gtf",
                        #default="../../examples/reference_gtf/human/ensembl/Homo_sapiens.GRCh38.96.chr.gtf",
                        #default = "../reference_gtf/human/ensembl/Homo_sapiens.GRCh38.96.chr.gtf",
                        #default="/home/niguilla/Documents/these_nguillaudeux/Data/gtf_hs/ensembl/Homo_sapiens.GRCh38.96.chr.gtf",
                        dest = "reference",
                        metavar = "GTF REFERENCE FILE",
                        help = "give a reference gtf file to validate prediction realize with CG-Alcode programm.")
              
    return parser.parse_args()                     

if __name__ == "__main__":
    OPTIONS = parse()
    gtf2bedscriptconvertor = OPTIONS.gtfbedscript
    ref_gtf_file = OPTIONS.reference
    
    gtf2bed(ref_gtf_file, gtf2bedscriptconvertor, False)