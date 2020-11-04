# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 13:21:22 2020

@author: s4426986

Use regex to find the flanking regions of UMIs (6x/8x) and the library 
insert (XXX)7. UMIs and inserts are added to their respective dic (only
unique UMIs and inserts, i.e. all keys in dic have value 1, not a tally).

Works with the reads from a single barcoded sample (text file, just reads)

A translate function is used to translate the inserts into their peptide 
(7 aa) sequence and added to the peptide_dictionary 

Due to known indels in MinION reads (usually 1-3 bases), a list is used 
to filter for random sequences of roughly the expected length (i.e. 6, 8 and 21). 
This filters extremely long/short sequences from the dictionaries, and is more 
lenient than exact matching. To make filtering more stringent, change the list
to select for exact length matches only.

At the end, the script reports the normalised numbers of UMIs, inserts and peptides.
Normalised = divided by number of reads in sample

"""



import re
# import csv


#function to reverse complement sequences, where for_or_r of oligo = 'r'
def ReverseComplement(oligo):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N'}
    return (str("".join([seq_dict[base] for base in reversed(oligo)])))

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def translate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
    return protein 


just_reads="good_reads_just_reads.txt"

umi_dictionary={}
umi_length_dictionary={}
insert_dictionary={}
insert_length_dictionary={}

## 1-, 2- and 3-base indels common in MinION reads
allow_umiLength=[8] #ideally just 6x and 8x
allow_insertLength=[21] #ideally just 21x 

for line in open(just_reads):
    s=line
    result1= re.search('ACACTGACGACATGGTTCTACA(.*)GTATG',s) #umi in sense strand
    # result2= re.search('CCATA(.*)TGTAG',s) #umi in antisense strand
    result3= re.search('CCGCA(.*)ACAAG',s) #insert in sense strand
    # result4= re.search('CTTGT(.*)TGCGG',s) #insert in antisense strand
    if result1:
        umi_pattern=result1.group(0)
        umi=umi_pattern[22:-5]
        umi_length=len(umi)
        print(umi_pattern)
        print(umi)
        print(umi_length)
        
        if umi not in umi_dictionary and len(umi) in allow_umiLength:
            umi_dictionary[umi]=1
        elif umi in umi_dictionary and len(umi) in allow_umiLength:
            umi_dictionary[umi]+=1
        else:
            pass
        
        if umi_length not in umi_length_dictionary:
            umi_length_dictionary[umi_length]=1
        elif umi_length in umi_length_dictionary:
            umi_length_dictionary[umi_length]+=1
            
    # elif result2:
    #     umi_pattern2=result2.group(0)
    #     print(umi_pattern2)
    #     umi2=ReverseComplement(umi_pattern2[5:-5]) #REVCOMP and check in dic
    #     print(umi2)
        
    #     if umi2 not in umi_dictionary and len(umi2) in umiLength:
    #         umi_dictionary[umi2]=1
    #     # elif umi2 in umi_dictionary and len(umi2) in umiLength:
    #     #     umi_dictionary[umi2]+=1
    #     else:
    #         pass
            
    elif result3:
        insert_pattern=result3.group(0)
        insert=insert_pattern[5:-5]
        insert_length=len(insert)
        print(insert_pattern)
        print(insert)
        print(insert_length)
        
        
        if insert not in insert_dictionary and len(insert) in allow_insertLength:
            insert_dictionary[insert]=1
        elif insert in insert_dictionary and len(insert) in allow_insertLength:
            insert_dictionary[insert]+=1
        else:
            pass
        
        if insert_length not in insert_length_dictionary:
            insert_length_dictionary[insert_length]=1
        elif insert_length in insert_length_dictionary:
            insert_length_dictionary[insert_length]+=1
    
    # elif result4:
    #     insert_pattern2=result4.group(0)
    #     print(insert_pattern2)
    #     insert2=ReverseComplement(insert_pattern2[5:-5]) #REVCOMP and check in dic
    #     print(insert2)
        
    #     if insert2 not in insert_dictionary and len(insert2) in insertLength:
    #         insert_dictionary[insert2]=1
    #     # elif insert2 in insert_dictionary and len(insert2) in insertLength:
    #     #     insert_dictionary[insert2]+=1
    #     else:
    #         pass
        


read_number= file_len(just_reads)
umi_per_kread= (1000*len(umi_dictionary))/read_number
insert_per_kread= (1000*len(insert_dictionary))/read_number

print('##############################################','\n')
print("Number of reads: ", read_number)
print("UMIs per 1000 reads: ", umi_per_kread)
print("Unique inserts (not peptides) per 1000 reads: ", insert_per_kread)


peptide_dictionary={translate(k): v for k, v in insert_dictionary.items()}
peptide_unique= len(peptide_dictionary)
peptide_per_kread= (1000*peptide_unique)/read_number
print("Unique peptides: ", peptide_unique)       
print("Number of unique peptides per 1000 reads: ", peptide_per_kread)


with open('good_reads_umi_dic.csv', 'w') as f:
    for key in umi_dictionary.keys():
        f.write("%s,%s\n"%(key,umi_dictionary[key]))
                   

# with open('105_insert_dic.csv', 'w') as f:
#     for key in insert_dictionary.keys():
#         f.write("%s,%s\n"%(key,insert_dictionary[key]))