#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 11:25:42 2018

@author: jennifermichaud
"""



def KmerComposition(dnastring, k):
    """Inputs a DNA string (dnastring) and a length(k).
    Outputs a list of DNA strings that deconvolute the string into a list of 
    k-lenth strings where each string has a k-1 base overlap with at least one 
    other DNA string and they are ordered in the list according to their overlap 
    such that the DNA string [i][1:] matches DNA string [i+1][:k].
    """
    kmers = []
    for n in range(len(dnastring)-k +1):
        kmer = dnastring[n:n+k]
        kmers.append(kmer)
    return kmers
    

def GenomePathReconstruct(genomepathlist):
    """Inputs a list of k-length DNA strings (genomepathlist) that have a k-1 
    base overlap with at least one other DNA string and they are ordered in the 
    list according to their overlap such that the DNA string [i][1:] matches DNA 
    string [i+1][:k].
    Outputs a DNA string that represents the reconstruction of the DNA strings
    using overlaps to determine sequence.
    """
    sequence = []
    sequence.append(genomepathlist[0])
    for i in range(1,len(genomepathlist)):
        sequence.append(genomepathlist[i][-1])
    return ''.join(sequence)
 
#TEST DATA

#k_in= 5
#dnastring_in='ACCGAAGCT'
#print("\n".join(KmerComposition(dnastring_in, k_in)))  
#print(GenomePathReconstruct(KmerComposition(dnastring_in, k_in)))  

          
##correct_out_KmerComposition = """ACCGA
#CCGAA
#CGAAG
#GAAGC
#AAGCT"""
#dna_in= inputdna.split("\n")
##correct_out_GenomePathReconstruct = "ACCGAAGCT"



