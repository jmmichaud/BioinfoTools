#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 08:44:11 2018

@author: jennifermichaud
"""
import numpy as np

def CountDNAnucleotides(DNAsequence):
    """ Inputs a string of DNA.  Outputs the number of each nucleotide
    ACGT respectively as a space seperated string of integers. """
    seq = DNAsequence.upper()
    countACGT = np.zeros((4), dtype = int)
    
    for let in seq:
        if let == 'A':
            countACGT[0]  += 1
        elif let == 'C':
            countACGT[1]  += 1
        elif let == 'G':
            countACGT[2]  += 1
        elif let == 'T':
            countACGT[3]  += 1
    return " ".join(map(str, countACGT))

