# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 18:51:37 2020

@author: s4426986

Go through SAM, look for duplicates and see if they map to same position in the reference.
Similar to testv1, but doesn't separate 0s and 16s into different files and reports only RNAME,
FLAG and POS, for ease of visualization'
"""


import os
# from collections import OrderedDict
import csv 

script = os.path.basename(__file__).split('.',1)[0]
script_dir = os.path.dirname(__file__)

samfile = "FAO02899_pass_7b15f441_24.sam"

with open(samfile, 'r') as sfile:
    unique=set()
    dupList = []
    for i, line in enumerate(sfile):
        
    
        line = line.rsplit() #takes line and returns a list of strings found in line
        #ignore the @SQ and @PG tag
        if line[0] == "@PG":
            pass
        elif line[0] == "@SQ":
            pass
        elif line[0] not in unique:
            # print(i)
            unique.add(line[0])
            
            rName = line[0]
            sFlag = line[1]
        
        elif rName in unique:   #if name already seen, append rName to DupList
            print("in line", i, "read ", rName, "already seen")
            dupList.append(rName)
              
            
with open(samfile, 'r') as sfile:
       
    for i, line in enumerate(sfile):
        line = line.rsplit()
        if line[0]=="@PG":
            pass
        elif line[0]== "@SQ":
            pass
        elif line[1]=="4":
            print(line[0], "unmapped")
            
        elif line[0] in dupList and line[1]=="0": #Stpre info for duplicates with FLAG 0
            # print(rName, "is forward sense alignment")
            
            
            
            rName = line[0]
            sFlag = line[1]
            # print(rName, "has 2ry alignment, with FLAG", sFlag)
            
            mPosition = line[3]
            
            flag0Dict = {"QNAME":rName, 
                                 "FLAG":sFlag,
                                 "POS":mPosition}
            
            ###########now we write info from dict into tab delimeted (SAM format) file
            fieldnames = ('QNAME', 'FLAG', 'POS')
            
            flag0result = [
                {'QNAME': rName, 'FLAG': sFlag, 'POS':mPosition}
            ]
            
            with open ('flag0easy.csv', 'a+') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
                # writer.writeheader() #don't print headers since SAMs don't have them, just values
                writer.writerows(flag0result)    
            with open('flag0easy.csv', 'r') as f:
                for line in f:
                    print(line.rstrip())
            ############################   
                        
with open(samfile, 'r') as sfile:
       
    for i, line in enumerate(sfile):
        line = line.rsplit()
        if line[0]=="@PG":
            pass
        elif line[0]== "@SQ":
            pass
        elif line[1]=="4":
            print(line[0], "unmapped")
            
        elif line[0] in dupList and line[1]=="16": #store infor for duplicates with flag 16
            
            rName = line[0]
            sFlag = line[1]
            print(rName, "has 2ry alignment, with FLAG", sFlag)
          
            mPosition = line[3]
            
            
            flag16Dict = {"QNAME":rName, 
                                 "FLAG":sFlag,
                                 "POS":mPosition}
                                
            
            ###########now we write info from dict into tab delimeted (SAM format) file
            fieldnames = ('QNAME','FLAG','POS')
            
            flag16result = [
                {'QNAME': rName, 'FLAG': sFlag, 'POS':mPosition}
            ]
            
            with open ('flag16easy.csv', 'a+') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
                # writer.writeheader() #don't print headers since SAMs don't have them, just values
                writer.writerows(flag16result)    
            with open('flag16easy.csv', 'r') as f:
                for line in f:
                    print(line.rstrip())
            #################