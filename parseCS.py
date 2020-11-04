# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 14:53:27 2020

@author: s4426986

Script to count the number of insertions and deletions.
Indels are counted as occurences first (regardless of the lentgh)
The magnitude (length) of deletions is assessed later.
Approach should be similar to the read length tally on lengthYqual.py
To avoid ambiguity, the script will ideally report everything with referenvce 
to the query (basecalled reads) rather than the reference (e.g. lambda genome)
In other words, an 'insert' reported by indel.py means that the basecalled read 
has a sequence that the reference doesn't. A 'deletion' menas that the read lacks
a sequence that the reference doesn't'

Use the --cs command from minimap2 to encode a flag similar to CIGAR that
reports indels, and use the script to parse the aligned SAM, look for the CS flag
and compute indel info from it

https://github.com/lh3/minimap2
The cs optional tag

The cs SAM/PAF tag encodes bases at mismatches and INDELs. It matches regular 
expression /(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)+/. Like CIGAR, cs consists of series of operations.
Each leading character specifies the operation; the following sequence is the one involved in the operation.

The cs tag is enabled by command line option --cs. The following alignment, for example:

CGATCGATAAATAGAGTAG---GAATAGCA
||||||   ||||||||||   |||| |||
CGATCG---AATAGAGTAGGTCGAATtGCA

is represented as :6-ata:10+gtc:4*at:3, where :[0-9]+ represents an identical block, -ata 
represents a deltion, +gtc an insertion and *at indicates reference base a is substituted with 
a query base t. It is similar to the MD SAM tag but is standalone and easier to parse.

If --cs=long is used, the cs string also contains identical sequences in the alignment. 
The above example will become =CGATCG-ata=AATAGAGTAG+gtc=GAAT*at=GCA. 
The long form of cs encodes both reference and query sequences in one string. 
The cs tag also encodes intron positions and splicing signals 
(see the minimap2 manpage for details).

*** Long story short *****
Use --cs=long on minimap2 (from guppy4 package) to keep bases that matched (capital ACGT), 
that way I don't have to index my reference or translate (unlike mismatch.py)
I basically parse the CS tag and count '+' and '-' (if I'm only interested in numbr of occurrences)
I can count the lower case letters that follow if I want to
compute length of ins. or del.
"""


from collections import OrderedDict
# import re
import csv

# Python code to convert string to list character-wise 
def Convert(string): 
    list1=[] 
    list1[:0]=string 
    return list1 
# # Driver code 
# str1="ABCD"
# print(Convert(str1)) 

    

samfile = "first10.sam"


with open(samfile, 'r') as sfile:
    csDic = {}
    tallyIns = {}
    tallyDel = {}
    tallyMis = {}
    
    for i, line in enumerate(sfile):
        line = line.rsplit()
        #ignore the @SQ and @PG tag
        if line[0] == "@SQ":
            pass
        elif line[0] == "@PG":
            pass
        else:
            #print (i)
            
            rName = line[0]
            sFlag = line[1]
            cName = line[2]
            mPosition = line[3]
            mQuality = line[4]
            cigar = line[5]
            mName = line[6]
            mPosition = line[7]
            tLen = line[8]
            rSeq = line[9]
            rQuality = line[10]
            
            ## additional info from line[11] onwards
            addInfoDict = OrderedDict() #call empty dict of additional information
            addInfo = [aline.split(":") for aline in line[11:]]
            for aline in addInfo:
                #addTag, addType, addValue = aline.split(":")
                if aline[1] == "i": #type = int
                    addInfoDict[aline[0]] = int(aline[2])
                elif aline[1] == "f": #type = float
                    addInfoDict[aline[0]] = float(aline[2])
                else: #all others are sting
                    addInfoDict[aline[0]] = str(aline[2])
                    
            if "cs" in addInfoDict and len(rSeq) != 1:
                 csTag = addInfoDict["cs"]
                 length = len(csTag)
                 match = 0
                 special = ['+', '-', '*', '=']
                 indel = ""
                 
                 for i in csTag:
                     if i.isupper():
                         match = match + 1  
                        
                     elif i.isupper() == False and i not in special:
                         indel = indel + i
                             
                 print(match, "out of", length, "bases match the reference", "\n")       
                 print(indel, "\n")
                 print("indel length = ", len(indel), "\n")
                 
                    
                 csList = Convert(csTag) #conver CS to list
                 csDic[rName] = csList #add list to dic
                 tallyIns[rName] = csList.count("+") #tally of ins in read
                 tallyDel[rName] = csList.count("-") #tally of dels in read
                 tallyMis[rName] = csList.count("*") #tally of mismatches in read
    

results = open("outInsertions.csv", "w")
writer = csv.writer(results)
for key, value in tallyIns.items():
    writer.writerow([key,value])
results.close()

results = open("outDeletions.csv", "w")
writer = csv.writer(results)
for key, value in tallyDel.items():
    writer.writerow([key,value])
results.close() 

results = open("outMismatches.csv", "w")
writer = csv.writer(results)
for key, value in tallyMis.items():
    writer.writerow([key,value])
results.close()

########## NOW COMPUTE LENGTH OF INDELS #############
