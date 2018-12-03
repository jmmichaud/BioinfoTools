#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 14:43:49 2018

@author: jennifermichaud
"""

#Find DnaA box within a genome

import sys

#Functions
def GCSkew(text): 
    """Inputs DNA as a text string.  Scores each G base +1 and each C as -1. 
    Outputs score at each base as a list of integers.
    """
    skew = []
    skew_calc = 0
    skew.append(skew_calc)
    text = text.upper()
    for base in text:
        if base == 'C':
            skew_calc -= 1
        elif base == 'G':
            skew_calc += 1
        skew.append(skew_calc) 
    return skew   

def MinSkewPos(text):
    """finds the position of minimum G-C skew in a genome 
    in order to find the ORI.
    G-C skew turns from decreasing to increasing at the ORI
    the position outputs the position starting at 1 (not std python 
    zero indexing).
    """
    skew = GCSkew(text)
    min_skew = min(skew)
    position = []
    for i in range(len(skew)):
        if skew[i] == min_skew:
            position.append(i)
    return position

def PatternCount(text, pattern):
    """Inputs a text string and a pattern string to be found within the text.
    Outputs an integer representing the number of times the pattern is found.
    """
    count = 0
    test = ""
    for i in range(len(text) - len(pattern)+1):
        test = text[i:i + len(pattern)]
        if test == pattern:
            count += 1
        test = ""
    return count

def ReverseComplement(text): 
    """Input a DNA text string.  Outputs the reverse complement of the string.
    """
    compd = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    revtext = text[::-1]
    revcomp = []
    for let in revtext:
        revcomp.append(compd[let])
    return ''.join(revcomp)

def PatternToNumber(DNApattern): 
    """Inputs a DNA text string.  Outputs a integer representation of the DNA
    sequence using quaternary encoding.
    """
    DNAcode = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    k = len(DNApattern)
    total = 0 
    for let in DNApattern:
        total += DNAcode[let] * (4 ** (k-1))
        k -= 1
    return total 

def NumberToPattern(dnadex, k): 
    """Inputs a integer dnadex that is a numerical representation of a DNA 
    sequence, and an interger kmer length k.
    Outputs the DNA sequence encodes by the quarternary code).
    """
    pos = k
    divfour = int(dnadex / 4)
    remfour = int(dnadex % 4)
    DNAcode = {0: 'A',1: 'C', 2: 'G',3: 'T'}
    dexcode = []
    DNAout = []
    while pos > 0:
        dexcode.insert(0, remfour)
        remfour = int(divfour % 4)
        divfour = int(divfour / 4)
        pos -= 1
    try:
        for item in dexcode:
            DNAout.append(DNAcode[item])
        return "".join(DNAout)
    except:
        return ""
    

def HammingDistance(text1,text2):
    """Computes the hamming distance between
    two equal length dna strings.  Each mismatch 
    is a score of 1"""
    hamming = 0
    for i in range(len(text2)):
        if text1[i] != text2[i]:
            hamming += 1
    return hamming

def Neighbors(pattern, d):
    """Inputs a pattern as a text string and a number of allowed mismatches, d. 
    Outputs all permutations of the pattern than have no more than d mismatches 
    of original pattern as a list of strings.
    """
    if d == 0:
        return pattern
    elif len(pattern)== 1:
        return ['A', 'C', 'G', 'T']
    code = ['A', 'C', 'G', 'T']
    neighborhood = []
    suffixneighbors = Neighbors(pattern[1:], d)
    for suff in suffixneighbors:
        if HammingDistance(pattern[1:], suff) < d:
            for base in code:
                neighborhood.append(base + suff)
        else:
            neighborhood.append(pattern[0] + suff)
    return neighborhood

def ComputingFrequenciesWithMismatches(text, k, d):
    """Inputs a DNA text string, a kmer length (k), and allowed number of 
    mismatches (d).  Using quaternary encoding of DNA strands generates a frequency
    array of all kmers its representations with d mismatches for each position
    of the DNA. Outputs the frequency array as a list of intergers where the index
    reprsents the quaternary code and the value at the index represents the 
    number of occurences of that kmer.
    """
    freq_array = []
    for i in range((4**k)):
        freq_array.append(0)
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        neighborhood = Neighbors(pattern, d)
        for item in neighborhood:
            j = PatternToNumber(item)
            if j != "":
                freq_array[j] += 1
    return freq_array
    
def FreqWordsMismatchRC(text, k, d):
    """Inputs a DNA text string, a kmer length (k), and allowed number of 
    mismatches (d). Using frequency arrays to compute the occurrence of all kmers
    of length k, in both forward and reverse complement orientations.  
    Outputs a list of frequent kmers
    """
    freq_array_F = ComputingFrequenciesWithMismatches(text, k, d)
    freq_array_RC = ComputingFrequenciesWithMismatches(ReverseComplement(text), k, d)
    freq_array_union = [freq_array_F[i] + freq_array_RC[i] for i in range(len(freq_array_F))]
    max_count = max(freq_array_union)
    freq_kmers = []
#    kmerdict = {}
    for i in range(len(freq_array_union)):
        if freq_array_union[i] == max_count:
            freq_kmers.append(NumberToPattern(i, k))  
    return freq_kmers

#Program to find DnaA box within a genome.  
#Program was tested with a Salmonella enterica genome.
    
#Welcome and initial prompt.
print("\n\n\n===================================================================================\nHello")
print("This program takes a small genome and finds the highest occuring 9-mers near the origin of replication allowing for 1 mismatch.  It generates candidate sequences that could represent DnaA boxes.  It looks in both the forward and reverse directions.")
print("===================================================================================")
print("\n\nThe genome you provide should be in .txt format.")
path = input("Please enter the genome filename.\nInclude full file path or type 'exit' to leave program: ")

    
#load genome
genome = ""
while genome == "":
    if path == 'exit':
        sys.exit("\n\nGoodbye\n\n")
    try:
        f = open(path, 'r') 
        #f = open('/Users/jennifermichaud/OneDrive/Coding Reference/HiddenDNA Coursera/Salmonella_enterica.txt', 'r') 
        genome = f.read()
        f.close()
        genome = genome.replace("\n", "")
        genome = genome.upper()
    except:
        if genome == "":
            print("\nThe filename or path you entered seems to be incorrect.")
            path = input("Please enter the genome filename including full file path: " )




# Find position of minimum skew to find possible ORI region, 1 indexed
min_skew = MinSkewPos(genome)

print("\n\n\nThe position of minimum G-C skew occurs at position(s): " + str(min_skew))
print("\nThe minimum skew generally indicates the origin of replication and is a good choice for starting position when choosing a window to look for DnaA boxes.")

#Promt to use default minimum skew position as begining of search window or to 
#enter manually. Default uses the putative origin of replication.  Also checks 
#to ensure valid input and allows program exit.
glen = len(genome)
autowindow = ""
while autowindow != 'y' and autowindow != 'n':
    print("\nWould you like to use the first minimum skew position as your starting position for the search window?")
    autowindow = input("Type 'y', 'n', or 'exit: ")
    if autowindow == 'y':
        startpos = min_skew[0]
    elif autowindow == 'n':
        starpos = ''
        try:       
            startpos = int(input("Please enter an integer start position between 0 and " + str(glen) + ": ")) 
        except:
            pass
        while type(startpos) != int or startpos > glen - 100: #check for valid input
            if type(startpos) != int:
                print("\nYou did not enter an integer.")
                try:       
                    startpos = int(input("Please enter an integer start position between 0 and " + str(glen) + ": ")) 
                except:
                    pass
            if startpos > glen - 100: #check if starpos is within genome
                print("\n" +"The start position is outside the length of the genome.")
                try:       
                    startpos = int(input("Please enter an integer start position between 0 and " + str(glen) + ": ")) 
                except:
                    pass
            if startpos == 'exit':
                sys.exit("\n\nGoodbye\n\n")
    elif autowindow == 'exit':
        sys.exit("\n\nGoodbye\n\n")
    else:
        print("You did not type a valid command.")
        

#Prompt to use default window size of 500bps or enter manually. Currently allows 
#for window sizes between 100-3000bps.  This can be adjusted as necessary.  
#Checks for valid input and allows for program exit.        
windowsize = 500
print("\n\n\nThe window size is 500 bp by default.")
autosize = ""
while autosize != 'y' and autosize != 'n':
    print("\nWould you like to use the 500 bp default window size?")   
    autosize = input("Type 'y', 'n', or 'exit: ")
    if autosize == 'y':
        endpos = startpos + windowsize
    elif autosize == 'n':
        windowsize = ""
        try:
            windowsize = int(input("\nPlease enter an integer window size between 100 - 3000 bp: "))
        except:
            pass
        while type(windowsize) != int or windowsize >3000 or windowsize < 100 or startpos + windowsize > glen: #check for valid input
            if type(windowsize) != int:
                print("\nYou did not enter an integer.")
                try:
                    windowsize = int(input("\nPlease enter an integer window size between 100 - 3000 bp: "))
                except:
                    pass
            if windowsize >3000 or windowsize < 100: 
                print("\nYou did not enter a valid integer.")
                try:
                    windowsize = int(input("\nPlease enter an integer window size between 100 - 3000 bp: "))
                except:
                    pass
            if startpos + windowsize > glen:
                print("\nThis window size is greater than the genome length.")
                print("\nThe max window size for the given position is " + str(glen - startpos) +".")
                try:
                    windowsize = int(input("\nPlease enter an integer window size between " + str(glen - startpos) + " bps: "))
                except:
                    pass
            if windowsize == 'exit':
                sys.exit("\n\nGoodbye\n\n")
                
        endpos = startpos + windowsize
    elif autowindow == 'exit':
        sys.exit("\n\nGoodbye\n\n")
    else:
        print("You did not type a valid command.")

#Define the search region (defaults will approximate a putative origin of replicaiton)
ORIregion = genome[startpos -1:endpos -1]

#Find frequent 9-mers(k = 9) in the search region with 1 mismatch (d= 1)
#includes reverse complements.  To find kmers in other scenarios the parameters are easily changable.
Freq_9mers = FreqWordsMismatchRC(ORIregion, 9, 1)
print("\n\n\n=======================================================================")
print("The most frequent 9-mers in this " + str(windowsize) + " region are: " + " ".join(Freq_9mers))



#find exact counts and reverse complement counts for the freq kmers
kmer_counts = {}
count_list = []
for kmer in Freq_9mers:
    count_list = [PatternCount(ORIregion, kmer),PatternCount(ORIregion, ReverseComplement(kmer))]
    kmer_counts[kmer] = count_list
print("\n\n\nKMER OCCURENCES OF EXACT MATCHES IN WINDOW\nUsing the format - 'kmer': [#forward occurences, #reverse occurences]:\n")
print(kmer_counts)

#find occurences of kmers with mismatches in forward and reverse complements.
mismatchkmer_counts = {}
count_list = []
for kmer in Freq_9mers:
    fcount = 0
    rccount = 0
    fkmerneighbors = Neighbors(kmer, 1)
    rckmerneighbors = Neighbors(ReverseComplement(kmer), 1)
    for fnei in fkmerneighbors:
        fcount += PatternCount(ORIregion, fnei)
    for rcnei in rckmerneighbors:
        rccount += PatternCount(ORIregion, rcnei)
    mismatchkmer_counts[kmer] =[fcount, rccount]
print("\n\n\nKMER OCCURENCES WITH ALLOWED MISMATCHES IN WINDOW\nUsing the format - 'kmer': [#forward occurences, #reverse occurences]:\n")
print (mismatchkmer_counts)   

print("\n\n\nKmers with most occurences in forward and reverse directions and are palindromes of one another are likely candidates for DnaA boxes. Multiple DnaA boxes and origins of replication may be present in a genome.") 

    
#TESTING NOTES - Salmonella enterica genome
#When using default settings there are two positions 3818639 3818641 of minimum skew.
#The first position is used at the beginning of the search window.

#The most frequent 9-mers are: CCAGGATCC CCCGGATCC CCGGATCCG CGGATCCGG 
#GGATCCGGG GGATCCTGG. 

#KMER OCCURENCES OF EXACT MATCHES IN WINDOW
#Using the format - 'kmer': [#forward occurences, #reverse occurences]:
#kmer_counts = {'CCAGGATCC': [0, 1], 'CCCGGATCC': [1, 1], 'CCGGATCCG': [0, 0],
# 'CGGATCCGG': [0, 0], 'GGATCCGGG': [1, 1], 'GGATCCTGG': [1, 0]}

#KMER OCCURENCES WITH ALLOWED MISMATCHES IN WINDOW
#Using the format - 'kmer': [#forward occurences, #reverse occurences]:
#mismatchkmer_counts = {'CCAGGATCC': [2, 3], 'CCCGGATCC': [2, 3], 'CCGGATCCG': [2, 3],
# 'CGGATCCGG': [3, 2], 'GGATCCGGG': [3, 2], 'GGATCCTGG': [3, 2]}

#The most frequent 9-mers are: CCAGGATCC CCCGGATCC CCGGATCCG CGGATCCGG 
#GGATCCGGG GGATCCTGG in the region of 500 bases starting at the position 
#of lowest G-C skew accounting for 1 mismatch and reverse complements.  
#These all contain GATCC which which seems an important motif and the similarity
#between all the found sequence may mean they all represent DnaA boxes. 'CCCGGATCC' 
#with its reverse complement 'GGATCCGGG' is found without mismatches  
# making it a prime candidate.



