#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 14:13:49 2018

@author: jennifermichaud
"""

def ShortFreqArray(text, k): 
   """Inputs a string of text and an integer kmer length k.  Generates a 
   frequency array  as a dictionary where the key is the kmer and the entry 
   is the count of that kmer in the text.  Avoids generation of every possible 
   kmer of length k in the frequency array.""" 
   kmer = ""
   Freq_array = {}
   for i in range(len(text)+1 - k):
        kmer = text[i:i + k]
        if kmer in Freq_array:
            Freq_array[kmer] += 1
        else:
            Freq_array[kmer] = 1
   return Freq_array