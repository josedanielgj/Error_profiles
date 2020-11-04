# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 18:31:29 2020

@author: s4426986

Go through the text file with reads line by line, find position 
of cs1/2 and their rev.comp., then use their positions to splice reads from
cs1 to cs2_rc (sense) and from cs2 to cs1_rc (antisense). Then apply rev.comp.
function to antisense reads to write all the reads in the sense orientation
in the output fasta file. 

As we do this, get the 13 bases before/after CS2 and add them to the 
barcode_dictionary. For antisense reads, the 13 bases are first revComp'd
then checked against barcode_dictionary

"""


import fileinput

#function to reverse complement sequences, where for_or_r of oligo = 'r'
def ReverseComplement(oligo):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G', 'a':'t', 't':'a', 'c':'g', 'g':'c', 'N':'N'}
    return (str("".join([seq_dict[base] for base in reversed(oligo)])))


just_reads="300-400bp_just_reads.txt"

output='good_reads.fasta'

just_reads_length_dictionary={}
with open(just_reads, 'r') as infile:
    for line in infile:
        read_length=len(line)
        # print(read_length)
        if read_length not in just_reads_length_dictionary:
            just_reads_length_dictionary[read_length]=1
        elif read_length in just_reads_length_dictionary:
            just_reads_length_dictionary[read_length]+=1




barcode_dictionary={}
good_read_length_dictionary={}

good_read_count=0

with open(just_reads,'r') as infile, open(output, 'w') as outfile:
    for line in infile:
        cs1_pos=line.find('ACACTGACGACATGGTTCTACA')
        cs2_rc_pos=line.find('AGACCAAGTCTCTGCTACCGTA')
        cs2_pos=line.find('TACGGTAGCAGAGACTTGGTCT')
        cs1_rc_pos=line.find('TGTAGAACCATGTCGTCAGTGT')
        
        if cs1_pos != -1 and cs2_rc_pos != -1: #sense reads
            good_read_count+=1
            good_read=line[cs1_pos:cs2_rc_pos+22] #cs1-to-cs2_rc
            # print(good_read, '\n')
            good_read_length=len(good_read)
            barcode=line[cs2_rc_pos+23:cs2_rc_pos+36] #BC 13 bases after cs2_rc
            
            if 248 < good_read_length < 261:
                outfile.write('>')
                outfile.write(str(good_read_count))
                outfile.write('\n')
                outfile.write(good_read.replace(' ',''))
                outfile.write(barcode.replace(' ',''))
                outfile.write('\n')
            else:
                pass
            
            if good_read_length not in good_read_length_dictionary:
                good_read_length_dictionary[good_read_length]=1
            elif good_read_length in good_read_length_dictionary:
                good_read_length_dictionary[good_read_length]+=1
            
            if barcode not in barcode_dictionary:
                barcode_dictionary[barcode]=1
            elif barcode in barcode_dictionary:
                barcode_dictionary[barcode]+=1
                
        elif cs2_pos != -1 and cs1_rc_pos != -1: #antisense reads
            good_read_count+=1
            good_read2=line[cs2_pos:cs1_rc_pos+22] #cs2-to-cs1_rc
            good_read2_rc=ReverseComplement(good_read2)
            # print(good_read2, '\n')
            # print(good_read2_rc, '\n')
            good_read2_length=len(good_read2_rc)
            barcode2=line[cs2_pos-14:cs2_pos-1] #BC 13 bases before cs2
            barcode2_rc=ReverseComplement(barcode2)
            
            if 248 < good_read2_length < 261:
                outfile.write('>')
                outfile.write(str(good_read_count))
                outfile.write('\n')
                outfile.write(good_read2_rc.replace(' ',''))
                outfile.write(barcode2_rc.replace(' ',''))
                outfile.write('\n')
            else:
                pass
            
            if good_read2_length not in good_read_length_dictionary:
                good_read_length_dictionary[good_read2_length]=1
            elif good_read2_length in good_read_length_dictionary:
                good_read_length_dictionary[good_read2_length]+=1
        
            
            
            if barcode2_rc not in barcode_dictionary:
                barcode_dictionary[barcode2_rc]=1
            elif barcode2_rc in barcode_dictionary:
                barcode_dictionary[barcode2_rc]+=1
        
        else:
            pass


for line in fileinput.FileInput(output,inplace=1):
    if line.rstrip():
        newline=line.rstrip("\n")
        print(newline)
       
# with open('barcode_dic.csv', 'w') as f:
#     for key in barcode_dictionary.keys():
#         f.write("%s,%s\n"%(key,barcode_dictionary[key]))