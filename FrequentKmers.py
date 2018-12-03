#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 22:08:03 2018

@author: jennifermichaud
"""

def PatternCount(text, Pattern):
    """Inputs a text string and a pattern or kmer to be found within the given text as a text string.
    Outputs a count of how often the given pattern is found in the text.
    """
    count = 0
    test = ""
    for i in range(len(text) - len(Pattern)+1): #pattern length accounted for to avoid out of range
        for n in range (len(Pattern)): #establishes piece of text to test pattern against of same length. slides 1 letter forward per loop
            test += text[i+n] 
        if test == Pattern:
            count += 1
        test = ""
    return count


def FrequentKmers(text, k):
    """Inputs a text string and an integer kmer length, k.  
    Outputs the most most frequent kmers of length k within the text using an 
    algorithm contained within the PatternCount() function to count the occurences 
    of kmers or patterns within the text. Will return ties (all kmers that are 
    repeated the discovered maximum occurences. Output is joined list as space 
    seperated text strings.
    """
    pattern = ""
    count = 0
    f_pattern = []
    max_count = 0
    for i in range(len(text)-k+1): #iterates through test accounting for k length to avoid out of range
        for n in range(k): #establishes pattern of length k
            pattern += text[i+n]
        count = PatternCount(text, pattern)#counts occurences of test pattern in text
        if count == max_count:  #adds pattern to list if there are multiple freq kmers 
            f_pattern.append(pattern)
        elif count > max_count: #if more frequent kmer deletes other previously found kmers and resets max count
            max_count = count
            f_pattern = []
            f_pattern.append(pattern)
        pattern = ""   #reset pattern 
    pattern_lst = list(set(f_pattern))  #eliminates duplicates 
    pattern_lst.sort()
    return " ".join(pattern_lst) 

#TEST DATA
#txt_in = 'TAACACGATCGGGCAAAGCGCACAACGGGCAAATAACACGATTATAGGGCTGTATAGGGCTGCGGGCAAATAACACGATCACGCGCTCGGGCAAACGGGCAAACACGCGCTGCGCACAACGGGCAAATAACACGATCGGGCAAACGGGCAAATAACACGATTATAGGGCTGCGGGCAAATATAGGGCTGCGGGCAAATAACACGATTAACACGATCGGGCAAATAACACGATCACGCGCTGCGCACAATATAGGGCTGTATAGGGCTGCACGCGCTTAACACGATTATAGGGCTGCGGGCAAACGGGCAAATAACACGATTAACACGATCGGGCAAACACGCGCTTATAGGGCTGCGGGCAAATAACACGATTATAGGGCTGTAACACGATTAACACGATTAACACGATCGGGCAAACGGGCAAACACGCGCTCGGGCAAAGCGCACAAGCGCACAATATAGGGCTGGCGCACAAGCGCACAATAACACGATCGGGCAAAGCGCACAACACGCGCTCGGGCAAACACGCGCTTATAGGGCTGCGGGCAAACACGCGCTTAACACGATCACGCGCTCACGCGCTCGGGCAAACGGGCAAATATAGGGCTGGCGCACAATAACACGATTAACACGATTATAGGGCTGTATAGGGCTGGCGCACAATAACACGATCGGGCAAATAACACGATCACGCGCTGCGCACAAGCGCACAAGCGCACAATAACACGATTATAGGGCTGCGGGCAAACACGCGCTGCGCACAAGCGCACAATATAGGGCTGCACGCGCTGCGCACAATAACACGATGCGCACAATATAGGGCTGTAACACGATTATAGGGCTGCGGGCAAAGCGCACAACACGCGCTTATAGGGCTGTAACACGATCACGCGCTGCGCACAAGCGCACAACGGGCAAACACGCGCTCACGCGCTTAACACGATTAACACGATCACGCGCT'
#k = 12
# 
#print(FrequentKmers(txt_in,k))


