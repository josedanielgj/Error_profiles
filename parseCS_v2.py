# -*- coding: utf-8 -*-
"""
Created on Sat Sep  5 17:46:42 2020

@author: s4426986

Version 2 notes:

Script to compute the number of indel
Computes length of insertions and deletions captured in lengthDic_insertions 
and lengthDic_deletions (output as csv files). Also length of matches.

By adding up LENGTH of ins, dels, matches + the number of mismatches (one * 
equals one op that consumes ref and query), you get the number of CS operations, 
to calculate % of total later (with excel).

Finds the position of the CS operation and reports the context.

Mismatches from CS tag come in *tg*at=ATGC... format, so computing length with 
same implementation would result in all mismatches were of length 1, so
think of another way to compute length of mismatches when consecutive "*" are 
found (e.g. *gt*ct*cg)
"""


from collections import OrderedDict
import re
import csv

# Python code to convert string to list character-wise 
def Convert(string): 
    list1=[] 
    list1[:0]=string 
    return list1 
# # Driver code 
# str1="ABCD"
# print(Convert(str1)) 

    

samfile = "hg38_align_cs.sam"


with open(samfile, 'r') as sfile:
    csDic = {} #number of keys in dic = number of reads with csTag
    tally_ins = {}
    tally_del = {}
    tally_mis = {}
    lengthDic_insertions = {}
    lengthDic_deletions = {}
    lengthDic_matches = {}
    context_dic_insertions = OrderedDict()
    context_dic_deletions = OrderedDict()
    context_dic_mismatches = OrderedDict()
    context_dic_matches = OrderedDict()
    
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
                 # cs_ops = set('+-*=')

                 for i, letter in enumerate(csTag):
                     if letter == "+":
                         ins_position = csTag.find("+") #find position of cs operation
                         ins_context = csTag[ins_position-30:ins_position+30] #report context of CS op
                         if ins_context in context_dic_insertions:
                             context_dic_insertions[ins_context]+=1 #count +1 to tally if ins_context already in dictionary
                         else:
                            context_dic_insertions[ins_context] = 1 #count first occurence
                         
                     elif letter == "-":
                         del_position = csTag.find("-") #find position of cs operation
                         del_context = csTag[del_position-30:del_position+30] #report context of CS op
                         if del_context in context_dic_deletions:
                             context_dic_deletions[del_context]+=1 #count +1 to tally if del_context already in dictionary
                         else:
                             context_dic_deletions[del_context] = 1 #count first occurence
                     elif letter == "*":
                         mis_position = csTag.find("*") #find position of cs operation
                         mis_context = csTag[mis_position-30:mis_position+30] #report context of CS op
                         if mis_context in context_dic_mismatches:
                             context_dic_mismatches[mis_context]+=1 #count +1 to tally if mis_context already in dictionary
                         else:
                            context_dic_mismatches[mis_context] = 1 #count first occurence
                            
                     # #### this loop doesn't make sense, just here for consistency (in,del,mis,match,in,del,mis,match...)    
                     # elif letter == "=":
                     #     match_position = csTag.find("=") #find position of cs operation
                     #     match_context = csTag[match_position-30:match_position+30] #report context of CS op
                     #     if match_context in context_dic_matches:
                     #         context_dic_matches[match_context]+=1 
                     #     else:
                     #        context_dic_matches[match_context] = 1
                     # ####################    
                     else:
                          pass
                         
                         
                 
                 ###### V2 addition ##################################################################
                 
                 cs_to_strings = re.split("(\=|\+|\-|\*)", csTag) #split into lst of strings based on delims
                 # print(infoStrings)
                 cs_to_strings.remove("") #remove leftmost empty string after aplitting
                 N = 2 #used to join consecutive strings in infoStrings (e.g. '+' with 'gcct' to '+gcct')
                 temp = '{}' * N  
                 sorted_cs = [temp.format(*ele) for ele in zip(*[iter(cs_to_strings)] * N)] #new list
                 for value in sorted_cs:
                     if value.startswith("+"):
                         
                         insertion_length = len(value)-1 # -1 because of "+" in string
                         
                         if insertion_length in lengthDic_insertions:
                             lengthDic_insertions[insertion_length]+=1 #tally +1 if ins_length already in dic
                         else:
                             lengthDic_insertions[insertion_length]=1 #else, count first occurrence
                         
                     elif value.startswith("-"):
                         
                          deletion_length = len(value)-1
                          
                          if deletion_length in lengthDic_deletions:
                              lengthDic_deletions[deletion_length]+=1 #tally +1 if del_length already in dic
                          else:
                              lengthDic_deletions[deletion_length]=1 #else, count first occurrence
                    
                     elif value.startswith("="):
                         
                          match_length = len(value)-1
                          
                          if match_length in lengthDic_matches:
                              lengthDic_matches[match_length]+=1 #tally +1 if match_length already in dic
                          else:
                              lengthDic_matches[match_length]=1 #else, count first occurrence
                              
                 #############################################################################             
                              
                 csList = Convert(csTag) #convert CS to list for counting
                 csDic[rName] = csList #add list to dic
                 tally_ins[rName] = csList.count("+") #tally of ins in read
                 tally_del[rName] = csList.count("-") #tally of dels in read
                 tally_mis[rName] = csList.count("*") #tally of mismatches in read       

### COMMENT OUT OUTPUT THAT YOU DON'T NEED ####
                 
results = open("out_Insertions_tally.csv", "w")
writer = csv.writer(results)
for key, value in tally_ins.items():
    writer.writerow([key,value])
results.close()

results = open("out_Insertions_length.csv", "w")
writer = csv.writer(results)
for key, value in lengthDic_insertions.items():
    writer.writerow([key,value])
results.close()

results = open("out_Insertions_context.csv", "w") 
writer = csv.writer(results)
for key, value in context_dic_insertions.items():
    writer.writerow([key,value])
results.close()

results = open("out_Deletions_tally.csv", "w")
writer = csv.writer(results)
for key, value in tally_del.items():
    writer.writerow([key,value])
results.close()

results = open("out_Deletions_length.csv", "w")
writer = csv.writer(results)
for key, value in lengthDic_deletions.items():
    writer.writerow([key,value])
results.close()

results = open("out_Deletions_context.csv", "w") 
writer = csv.writer(results)
for key, value in context_dic_deletions.items():
    writer.writerow([key,value])
results.close()

results = open("out_Mismatches_tally.csv", "w") 
writer = csv.writer(results)
for key, value in tally_mis.items():
    writer.writerow([key,value])
results.close()

results = open("out_Mismatches_context.csv", "w") 
writer = csv.writer(results)
for key, value in context_dic_mismatches.items():
    writer.writerow([key,value])
results.close()

results = open("out_Matches_length.csv", "w")
writer = csv.writer(results)
for key, value in lengthDic_matches.items():
    writer.writerow([key,value])
results.close()
