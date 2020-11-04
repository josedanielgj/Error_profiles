# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 16:49:16 2020

@author: s4426986

Go through SAM and look for duplicates, create list of dups and use the list to parse SAM again.
If ReadName in DUPLIST, store SAM info in dictionaries, one for FLAG:0 another one fior FLAG:16. 
Then report 0s and 16s into different txt output (a filtered SAM, tab delim.) that can be used 
with mismatch.py to extract error info from new SAM 

TO DO LIST:
Turn "with open" bit into function and use that for to open SAM twice (first to write 
flag0verbose.csv, then to write flag16verbose.csv)
"""


import os
from collections import OrderedDict
import csv 

script = os.path.basename(__file__).split('.',1)[0]
script_dir = os.path.dirname(__file__)

samfile = "sample10000_pm4_pm1.sam"

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
            
            # if sFlag == "16":
            #     print (sFlag)
        
            
        
            # print ("with sequence   ", line[9])
            # print ("and FLAG ", sFlag, "already seen")
    
   

with open(samfile, 'r') as sfile:
    unmapped = 0
    # Dict0={}
    # Dict16={}
    
    for i, line in enumerate(sfile):
        line = line.rsplit()
        if line[0]=="@PG":
            pass
        elif line[0]== "@SQ":
            pass
        elif line[1]=="4":
            
            print(line[0], "unmapped")
            unmapped += 1
            
        elif line[0] in dupList and line[1]=="0": #Store info for duplicates with FLAG 0
            # print(rName, "is forward alignment")
            
            
            
            rName = line[0]
            sFlag = line[1]
            print(rName, "has 2ry alignment, with FLAG", sFlag)
           
            cName = line[2]
            mPosition = line[3]
            mQuality = line[4]
            cigar = line[5]
            mName = line[6]
            mPosition = line[7]
            tLen = line[8]
            rSeq = line[9]
            rQuality = line[10]
            
            addInfoDict0 = OrderedDict() #call empty dict of additional info for 0-flagged reads
            addInfo0 = [aline.split(":") for aline in line[11:]] #lines 11-20 have format tag:type:value
            for aline in addInfo0:
                #addTag, addType, addValue = aline.split(":")
                if aline[1] == "i": #type = int
                    addInfoDict0[aline[0]] = int(aline[2])
                elif aline[1] == "f": #type = float
                    addInfoDict0[aline[0]] = float(aline[2])
                else: #all others are sting
                    addInfoDict0[aline[0]] = aline[2]
                    
                    NM = line[11]
                    ms = line[12]
                    AS = line[13]
                    nn = line[14]
                    tp = line[15]
                    cm = line[16]
                    s1 = line[17]
                    s2 = line[18]
                    dv = line[19]
                    MD = line[20]
            
            flag0Dict = {"QNAME":rName, 
                                 "FLAG":sFlag,
                                 "RNAME":cName,
                                 "POS":mPosition,
                                 "MAPQ":mQuality,
                                 "CIGAR":cigar,
                                 "RNEXT":mName,
                                 "PNEXT":mPosition,
                                 "TLEN":tLen,
                                 "SEQ":rSeq,
                                 "QUAL":rQuality,
                                 "NM":NM,
                                 "ms":ms,
                                 "AS":AS,
                                 "nn":nn,
                                 "tp":tp,
                                 "cm":cm,
                                 "s1":s1,
                                 "s2":s2,
                                 "dv":dv,
                                 "MD":MD}
            
            ###########now we write info from dict into tab delimeted (SAM format) file
            fieldnames = ('QNAME', 'FLAG','RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT',
                               'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NM', 'ms', 'AS',
                               'nn', 'tp', 'cm', 's1', 's2', 'dv', 'MD')
            
            flag0result = [
                {'QNAME': rName, 'FLAG': sFlag, 'RNAME':cName, 'POS':mPosition, 'MAPQ': mQuality, 'CIGAR':cigar,
                 'RNEXT':mName, 'PNEXT':mPosition, 'TLEN':tLen, 'SEQ':rSeq, 'QUAL':rQuality, 'NM':NM, 'ms':ms,
                 'AS':AS, 'nn':nn, 'tp':tp, 'cm':cm, 's1':s1, 's2':s2, 'dv':dv, 'MD':MD}
            ]
            
            with open ('flag0verbose.csv', 'a+') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
                # writer.writeheader() #don't print headers since SAMs don't have them, just values
                writer.writerows(flag0result)    
            with open('flag0verbose.csv', 'r') as f:
                for line in f:
                    print(line.rstrip())
            ############################   
                        
with open(samfile, 'r') as sfile:
    
    # Dict0={}
    # Dict16={}
    
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
           
            cName = line[2]
            mPosition = line[3]
            mQuality = line[4]
            cigar = line[5]
            mName = line[6]
            mPosition = line[7]
            tLen = line[8]
            rSeq = line[9]
            rQuality = line[10]
            
            addInfoDict16 = OrderedDict() #call empty dict of additional info for 0-flagged reads
            addInfo16 = [aline.split(":") for aline in line[11:]] #lines 11-20 have format tag:type:value
            for aline in addInfo16:
                #addTag, addType, addValue = aline.split(":")
                if aline[1] == "i": #type = int
                    addInfoDict16[aline[0]] = int(aline[2])
                elif aline[1] == "f": #type = float
                    addInfoDict16[aline[0]] = float(aline[2])
                else: #all others are sting
                    addInfoDict16[aline[0]] = aline[2]
                    
                    NM = line[11]
                    ms = line[12]
                    AS = line[13]
                    nn = line[14]
                    tp = line[15]
                    cm = line[16]
                    s1 = line[17]
                    s2 = line[18]
                    dv = line[19]
                    MD = line[20]
            
            flag16Dict = {"QNAME":rName, 
                                 "FLAG":sFlag,
                                 "RNAME":cName,
                                 "POS":mPosition,
                                 "MAPQ":mQuality,
                                 "CIGAR":cigar,
                                 "RNEXT":mName,
                                 "PNEXT":mPosition,
                                 "TLEN":tLen,
                                 "SEQ":rSeq,
                                 "QUAL":rQuality,
                                 "NM":NM,
                                 "ms":ms,
                                 "AS":AS,
                                 "nn":nn,
                                 "tp":tp,
                                 "cm":cm,
                                 "s1":s1,
                                 "s2":s2,
                                 "dv":dv,
                                 "MD":MD}
            
            ###########now we write info from dict into tab delimeted (SAM format) file
            fieldnames = ('QNAME', 'FLAG','RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT',
                               'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NM', 'ms', 'AS',
                               'nn', 'tp', 'cm', 's1', 's2', 'dv', 'MD')
            
            flag16result = [
                {'QNAME': rName, 'FLAG': sFlag, 'RNAME':cName, 'POS':mPosition, 'MAPQ': mQuality, 'CIGAR':cigar,
                 'RNEXT':mName, 'PNEXT':mPosition, 'TLEN':tLen, 'SEQ':rSeq, 'QUAL':rQuality, 'NM':NM, 'ms':ms,
                 'AS':AS, 'nn':nn, 'tp':tp, 'cm':cm, 's1':s1, 's2':s2, 'dv':dv, 'MD':MD}
            ]
            
            with open ('flag16verbose.csv', 'a+') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
                # writer.writeheader() #don't print headers since SAMs don't have them, just values
                writer.writerows(flag16result)    
            with open('flag16verbose.csv', 'r') as f:
                for line in f:
                    print(line.rstrip())
            #############################
        
        
          