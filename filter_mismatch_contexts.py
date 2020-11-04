# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 18:14:00 2020

@author: s4426986

Look for complimentarity/palindromic artifacts in mismatch tally 
(output from mismatch.py)
"""


import os
# from collections import OrderedDict

script = os.path.basename(__file__).split('.',1)[0]
script_dir = os.path.dirname(__file__)
## Assign path to FASTA file and SAM file

textResults = "filteredTriNuc.txt"

samfile = "FAO02899_pass_7b15f441_0.sam"

with open(textResults, 'r') as tfile:
    for i, line in enumerate(tfile):
        if line.startswith('A<-CGG'):
            print (line)
        elif line.startswith('T<-CCG'):
            print (line)
        elif line.startswith('A<-AGC'):
            print (line)
        elif line.startswith('G<-AAC'):
            print (line)
        
        
with open(samfile, 'r') as sfile:
    for i, line in enumerate(sfile):
    
        line = line.rsplit()
        if line[1] == "2064":
            sFlag = line[1]
            print (sFlag)
            
# test = "877be2c4-fd7e-4f4d-82f9-235f03910f3a"

# print (len(test))