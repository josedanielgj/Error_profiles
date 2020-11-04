# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 01:04:07 2020

@author: JDGJ

Version 1
   
Script to trim bases from the end of reads while keeping FASTQ format that
works with Galaxy tools for further data processing

Thanks John from BioStars for the 'toggle' approach. Useful to exploit 
the 4 lines per read format of FASTQ files. Keep in mind for future scripts.

https://www.biostars.org/p/188878/


version 1.2 (current version)

Process multiple fastQs at the same time (*.fastq) and output tailsBCXX for
each one

https://education.molssi.org/python_scripting_cms/03-multiple_files/index.html

Uses 'os' and 'glob' to:
    
Get multiple files from a directory, process them and write ouput to a new file (for each file
in directory) 
"""
# ############################# VERSION 1 ###########
# fastq = 'BC97.fastq' #change each time
# results = 'tailsBC97.fastq' #name should match BC above


# with open(fastq,'r') as infile, open(results, 'w') as outfile:
#     toggle = False
#     for line in infile:
#         if toggle:
#             print(line[-51:-1])
#             outfile.write(line[-51:-1])
#             outfile.write('\n')
#         else:
#             print(line[0:-1])
#             outfile.write(line[0:-1])
#             outfile.write('\n')
#         toggle = not toggle
# ################################################################

import os
import glob

## look for all FASTQs in a directory 
file_location = os.path.join('C:\\Users\\s4426986\\Desktop\\Jose', 'python', '*.fastq')
print(file_location, '\n')

filenames = glob.glob(file_location)
# print(filenames)

trimmedFiles = []
print('Trimmed files are:', '\n')
for file in filenames:
    file_name = os.path.basename(file)
    prefix = "tails"+file_name #adds 'tails' to file name
    outfile = "".join(prefix.split('.')[:-1]) + '.fastq' 
    ## redundant; removes file format (.fastq) and 
    ## adds ".fastq". I'll keep it in case Galaxy gives me .fastqsanger
    print(outfile)
    trimmedFiles.append(outfile)

    with open(file_name,'r') as infile, open(outfile, 'w+') as ofile:
        toggle = False
        for line in infile:
            if toggle:
                # print(line[-51:-1])
                ofile.write(line[-51:-1])
                ofile.write('\n')
            else:
                # print(line[0:-1])
                ofile.write(line[0:-1])
                ofile.write('\n')
            toggle = not toggle