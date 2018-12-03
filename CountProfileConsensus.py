#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 10:51:24 2018

@author: jennifermichaud
"""

import numpy as np
import math
import copy

def LaplaceCount(motif_list, k):
    """Inputs a list of k-length DNA strings (dna_list) and a kmer length (k).  
    Outputs a count matrix as a numpy array.  In the count matrix row 0 through
    3 represent the A, C, G, or T respectively and the columns represent the 
    position of a potential motif of length k. Initalizes a count of one to avoid 
    zero values that would interfere in the calculation of zero probabilities.
    """
    count_ACGT = np.ones((4, k))
    for i in range(len(motif_list)):  #generate counts for each base at position j for each motif[i]
        for j in range(k):
            if motif_list[i][j] == 'A':
                count_ACGT[0, j] += 1
            elif motif_list[i][j] == 'C':
                count_ACGT[1, j] += 1
            elif motif_list[i][j] == 'G':
                count_ACGT[2, j] += 1
            elif motif_list[i][j] == 'T':
                count_ACGT[3, j] += 1 
    return count_ACGT



def Profile(count_ACGT, motif_list, k):
    """Inputs a count matrix as a numpy array (see LaplaceCount()), a list of 
    k-length DNA strings, and a kmer length (k).  
    Outputs a probability matrix for each base being at each position, 0-k.
    The probability matrix is a numpy array where row 0 through 3 
    represent the probablity of finding an A, C, G, or T, respectively, and 
    each column corresponds to a position in the k-lenthg motif.
    """
    prob_ACGT = np.ones((4, k))
    for n in range(4): #generate probability for each base
        for l in range(k):
            prob_ACGT[n, l] = np.divide(count_ACGT[n,l],len(motif_list))
    return prob_ACGT


def Consensus(profile, k):
    """Inputs a probability matrix (profile) and a kmer length (k). 
    The probability matrix is a numpy array where row 0 through 3 
    represent the probablity of finding an A, C, G, or T, respectively, and 
    each column corresponds to a position in the k-length motif.
    Outpus a consensus sequence that is the highest probable motif given the 
    probability matrix.
    """
    consensus = []
    for i in range(k):
        maxpos = np.argmax(profile[0:4,i])
        if maxpos == 0:
            consensus.append('A')
        elif maxpos == 1:
            consensus.append('C')
        elif maxpos == 2:
            consensus.append('G')
        elif maxpos == 3:
            consensus.append('T')
    return ''.join(consensus)