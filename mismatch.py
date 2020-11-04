"""
Script to extract mismatch information from a SAM file. 

Make sure mapping program includes CIGAR string and MD optional tag in the SAM for script to work.

It also needs a FASTA reference file

Outputs two files, one with single-base mismatches, the other includes the flanking 
bases on either side of the mismatch (i.e. the context of the mismatch).

You can chnage the number of bases on the context (e.g. a trinucleotide or a 
pentanucleotide) by changing the splicing limits of triMismatchTally function.

For word = "banana", splicing indexing goes as follows;

| b | a | n | a | n | a |
0   1   2   3   4   5   6

where word[0:2] gives you "ba" and word[3:6] gives you "ana"

"""


import os
from collections import OrderedDict

script = os.path.basename(__file__).split('.',1)[0]
script_dir = os.path.dirname(__file__)
## Assign path to FASTA file and SAM file

fasta_reference_file = "sequence.fasta"
samfile = "mar27minimap.out.withMD.sam"


singleResults = "out_mismatch_single.txt"
trinucResults = "out_mismatch_pentanuc.txt"


samDict = {} 
dnaDict = {}    

def readFasta2Dict(Fasta):
    fastaDict = OrderedDict()
    
    with open(Fasta, 'r') as file:
        for line in file:
            line = line.rsplit()
            if line.startswith(">"):
                name = line
                fastaDict[name] = ""
            else:
                dna = line.upper()
                fastaDict[name] = dna
    return fastaDict


def translateCIGAR(CIGAR):
    cigarList = [
            "M", 
            "I", 
            "D", 
            "N", 
            "S", 
#            "H", 
            "P", 
            "=", 
            "X"
            ]
    
    Variables = [] 
    
#    print ("input CIGAR:", CIGAR)
            
    var = []
    for letter in CIGAR:

        
        if letter.isdigit() == True:
            var.append(letter)
        elif letter.isdigit() == False:
    #        var.append(letter) 
            if letter in cigarList: 
                subtotal = "".join(var)
                subtotal = int(subtotal)
#                Variables.append([subtotal, letter]) #add to list for debugging
                Variables.append("".join([letter for n in range(subtotal)]))
                
                var = [] #turn var back into blank list
            elif letter == "H": #deal with hardclip
                var = []
            elif letter == "P":
                pass
    
#    cigarSeq = Variables
    cigarSeq = "".join(Variables)    
    
#    return Variables
    return cigarSeq
	


def translateRefSeq(cigarSeq, refSeq):
    
    newSequence = []
    
    seqIndex = 0
    
    #walk along each letter of the cigar sequence
    for i, letter in enumerate(cigarSeq):
        if letter == "D":
            newSequence.append("*")
        else:
            try:
                newSequence.append(refSeq[seqIndex])
                seqIndex+=1
            except  IndexError:
                print ("cigar seq len:", len(cigarSeq,) , "seq Index:", seqIndex)
                print ("ref seq len:", len(refSeq), "seq Index:", seqIndex)
                break
    
    newSequence = "".join(newSequence)
	
    return newSequence
    


def translateMD(cigarSeq, mdtag):
    
    MDSeq = []
    
    #translate mdtag to seqence
    #remove all "^"
    mdtag = mdtag.replace("^","")
#    print (mdtag)
    
    var = []
    for letter in mdtag:
        #go through each variable in mdtag until you hit a letter
        if letter.isdigit() == True:        
            var.append(letter)
        elif letter.isdigit() == False:
            if var:
                subtotal = "".join(var)
#                print (subtotal)
#                print (letter)
                subtotal = int(subtotal)
                MDSeq.append("".join(["M" for n in range(subtotal)]))
                MDSeq.append(letter)
                
                var = []
            else: #if var is empty
                MDSeq.append(letter)
    #check if var is present
    if var:
        subtotal = "".join(var)
        subtotal = int(subtotal)
        MDSeq.append("".join(["M" for n in range(subtotal)]))
        
    
    MDSeq = "".join(MDSeq)
    
    
    newMDSeq = []
    seqIndex = 0
    
    #walk along each letter of the cigar sequence
    for i, letter in enumerate(cigarSeq):
        #first look for softclips
        if letter == "S":
            newMDSeq.append("*") #clips
        elif letter == "M":
            newMDSeq.append(MDSeq[seqIndex])
            seqIndex+=1
        elif letter == "D":
            newMDSeq.append(MDSeq[seqIndex].lower())
            seqIndex+=1
        elif letter == "I":
            newMDSeq.append("I")
            
            
       
    newMDSeq = "".join(newMDSeq)


    return newMDSeq



def dnaBaseCheck(base):
    
    dna = ["A", "T", "C", "G"]
    
    if base in dna:
        result = True
    else:
        result = False
    return result



def triMismatchTally (nRefSeq, nMDSequence):
    triTallyDict = OrderedDict()
    #parse through MD sequence to find index and base
    for i, base in enumerate(nMDSequence):
        if dnaBaseCheck(base) == True: #use dnaBaseCheck to find mismatches from MD sequence
            "if Mimsatch found in MD sequence (i.e. capital AGCT), move to refSeq AT THAT POSITION and"
            "report trinucleotide, where middle base is the mismatch"
            midBase = base 
            triNuc = nRefSeq[i-2:i+3] #report the 3 matching nRefSeq using [i-1:i+2] e.g. ATC
            mismatch = "".join([midBase, '<-', triNuc])
            
            if mismatch in triTallyDict:  #check if trinucleotide bases are in dictionary.
                triTallyDict[mismatch]+=1 #count +1 to tally if triNuc already in dictionary
            else:
                triTallyDict[mismatch] = 1 #count first occurence
                
    return triTallyDict   #Add to dictionary

def mismatchTally(nRefSeq, nMDSequence):
    tallyDict = OrderedDict()
    for ref_base, md_base in zip(nRefSeq, nMDSequence):
        
        if dnaBaseCheck(md_base) == True:
            
            bChange = "".join([ref_base, '<-', md_base])
            if bChange in tallyDict:
                tallyDict[bChange]+=1
            else:
                tallyDict[bChange] = 1
    return tallyDict



MisMatchDict = OrderedDict()
TriNucDict = OrderedDict()



with open(samfile, 'r') as sfile:
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
            
            #additional info from line[11] onwards
            addInfoDict = OrderedDict() #call empty dict of additional information
            addInfo = [aline.split(":") for aline in line[11:]]
            for aline in addInfo:
                #addTag, addType, addValue = aline.split(":")
                if aline[1] == "i": #type = int
                    addInfoDict[aline[0]] = int(aline[2])
                elif aline[1] == "f": #type = float
                    addInfoDict[aline[0]] = float(aline[2])
                else: #all others are sting
                    addInfoDict[aline[0]] = aline[2]
					
                    
            if "MD" in addInfoDict and len(rSeq) != 1:
                mdTag = addInfoDict["MD"]
                cigarSeq = translateCIGAR(cigar)
                newSeq = translateRefSeq(cigarSeq, rSeq)
                nMDSequence = translateMD(cigarSeq, mdTag)
                tallyDict = mismatchTally(newSeq, nMDSequence)
                triTallyDict = triMismatchTally(newSeq, nMDSequence)
                
				
                #add the mismatch tally from each 
                for cbase in tallyDict:
                    if cbase in MisMatchDict:
                        MisMatchDict[cbase]+=tallyDict[cbase]
                    else:
                        MisMatchDict[cbase] =tallyDict[cbase]

                      
                for cbase in triTallyDict:
                    if cbase in TriNucDict:
                        TriNucDict[cbase]+=triTallyDict[cbase]
                    else:
                        TriNucDict[cbase] =triTallyDict[cbase]

            elif len(rSeq) == 1: #when sequence is not stored, rSeq == "*"
                print (rName)    


with open(singleResults, 'w') as rfile:
    for cbase in MisMatchDict:
        rline = [cbase, str(MisMatchDict[cbase])]
        rline = "\t".join(rline)+"\n"
        rfile.write(rline)

with open(trinucResults, 'w') as rfile:
    
    for cbase in TriNucDict:
        rline = [cbase, str(TriNucDict[cbase])]
        rline = "\t".join(rline)+"\n"
        rfile.write(rline)
