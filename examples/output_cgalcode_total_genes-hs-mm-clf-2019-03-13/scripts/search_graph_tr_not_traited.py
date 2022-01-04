#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 11:04:10 2019

@author: niguilla
"""

with open("output_run5.txt", 'r') as o_5:
    lecture = o_5.readlines()
    for line in lecture:
        line = line.replace('\n', '')
        if line.startswith(">>> Genes:") or line.startswith(">>> !!!Graph not conserved!!!"):
            print (line)