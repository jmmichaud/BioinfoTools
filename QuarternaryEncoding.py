#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 16:20:41 2018

@author: jennifermichaud
"""

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
    Outputs the DNA sequence encoded by the quarternary code.
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