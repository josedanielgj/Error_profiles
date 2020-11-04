# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 22:01:50 2020

@author: s4426986

Version1.2 notes:
Script to demultiplex FLDXXXX barcodes from fastq filea after guppy HAC basecalling

Mission aborted!

V1.3 (current version):

Script gets readLength for all reads in fastq, adds it to dict, then writes dict as .csv
 
Order reads chronologically,
extract quality string of each read, translate ascii to ordinal and
concatenate a string that refelects change of quality throughout seq run

# NB: number of keys in timeDict is NOT equal to the number of reads in fastq because some reads sequenced at the same time,
so keys are not unique. It still reflects change of quality over time.

# NB : if working with minKNOW's fastQs, account for extra 'flow_cell_id' line
in first line of each read (not present in standalone Guppy output). Start_time becomes the
5th (instead of 6th) element of line (so [4] instead of [5] after rsplit()) 

"""


import os 
import csv

############ TURN QSCORE TO LIST OF CHARACTERS IN DICT VALUE
def split(word): 
    return [char for char in word] 
###############################################
    
############################## FUNCTION TAKEN FROM EUGENE'S SCRIPT TO REV.COMP. SEQUENCES
def reverse_complement(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)
    
    seq = seq[::-1].translate(translation_table)
    
    return seq
###################################################
     

script = os.path.basename(__file__).split('.',1)[0]
script_dir = os.path.dirname(__file__)

fastq = 'BC97.fastq'

primerInfo = {'FLD0097':'TGCTACATCA', 'FLD0098':'AGTGTGTCTA', 
                'FLD0099':'TCATATCGCG', 'FLD0100':'TACGTATAGC',
                'FLD0101':'CAGCTATAGC', 'FLD0102':'TCGATGCGCT', 
                'FLD0103':'GCACGCGTAT','FLD0104':'GCAGTATGCG',
                'FLD0105':'TGATAGAGAG','FLD0106':'GCTACTAGCG',
                'FLD0107':'TGCGAGACGT','FLD0108':'CGATGACAGA',
                'FLD0109':'GACTCATGCT','FLD0110':'GTCTGATACG',
                'FLD0111':'ACTAGCTGTC','FLD0112':'GCGTAGACGA',
                'P5_28b_6X_CS1':'ACACTGACGACATGGTTCTACANVNVNVGTATGGGATTTTGCTAAACAAC',
                'P5_28b_8X_CS1':'ACACTGACGACATGGTTCTACANBNVNVNVGTATGGGATTTTGCTAAACAAC',
                'P7_28gIII_CS2':'TACGGTAGCAGAGACTTGGTCTCCTCGAAAGCAAGCTGAT',
                'CS1':'ACACTGACGACATGGTTCTACA','CS2':'TACGGTAGCAGAGACTTGGTCT',
                'P5_IllumFor_HPLC':'AATGATACGGCGACCACCGA','P7_IllumRev_HPLC':'CAAGCAGAAGACGGCATACGA'}

barcode_list = list(primerInfo.values()) # turn dict into list


## make a dictionary of the reverse complementsa of the primers
revCompDict = {}
for k, v in primerInfo.items():
    revCompDict[k + '_rc'] = reverse_complement(v)
    
count1 = 0 #name        
count2 = 1 #sequence
count3 = 2 #separator
count4 = 3 #quality

##### SECOND BLOCK OF FOUR LINES ######
#### count1 = 4 #name            #####
#### count2 = 5 #sequence       #####
#### count3 = 6 #separator      #####
#### count4 = 7 #quality        #####

with open(fastq, 'r') as file:
    timeDict = {}
    lengthDict = {}
    headDict = {}
    tailDict = {}
    
    for i, line in enumerate(file):
        
        # check for name
        if i == count1: #name and other info
            
            line = line.rsplit()
            
            name = line[0]
            time = line[5]
            
            # print(line, i, count1)

            count1+=4
            
        elif i == count2: #sequence
            
            head = line[:30] #take first 30 bases in sequence by slicing
            tail = line[-29:] #take last 30 bases in sequence with negatuve index
            headDict[name] = head #add seq head to headDict
            tailDict[name] = tail #add tail to tailDict
            
            
            readLength = len(line)
            
            if readLength in lengthDict:
                lengthDict[readLength]+=1 #tally +1 if length already in dic
            else:
                lengthDict[readLength]=1 #else, count first occurrence
                
            # print (line.split()[0], i, count1)
            # pass
                
            count2+=4

        elif i == count3: #separator '+'
            
            # print (line.split()[0], i, count1)
            # pass
            count3+=4
            
        elif i == count4: #quality string
            
            # print (line.split()[0], i, count1) #count1/4 = number of lines in fastq
            # pass
            timeDict[time] = split(line.split()[0])
            "add qScore as value to timeDict"
            "and use split function from above to return each character in quality string as a list"
            "to then use ord() to translate each character from ascii to ordinal"
            count4+=4          
            
results = open("lengthResults.csv", "w")
writer = csv.writer(results)
for key, value in lengthDict.items():
    writer.writerow([key,value])
results.close()

results = open("headDict.csv", "w")
writer = csv.writer(results)
for key, value in headDict.items():
    writer.writerow([key,value])
results.close() 

results = open("tailDict.csv", "w")
writer = csv.writer(results)
for key, value in tailDict.items():
    writer.writerow([key,value])
results.close() 
            
####### Translate ASCII values from timeDict to ordinal, reports 2 list (l1 is ASCII, l2 is ordinal)  ###

d1 = timeDict
for v in d1:

    l1 = d1[v]
    l2 = []
    for c in l1:
        # l1.append(c) 
        l2.append(ord(c))
    # print(l1)
    # print(l2)       

##### TURN L1 AND L2 TO KEYS AND VALUES IN QUALDICT #######     
keys = list(timeDict.keys())
values = l2
qualDict = dict(zip(keys, values))
# print(qualDict) 
############################################################

results = open("qualResults.csv", "w")
writer = csv.writer(results)
for key, value in qualDict.items():
    writer.writerow([key,value])
results.close()