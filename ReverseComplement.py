#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 10:41:48 2018

@author: jennifermichaud
"""

def ReverseComplement(text):
    """ Inputs a DNA string.  Outputs the reverse complement of the DNA strand."""
    text = text.upper()
    compd = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    revtext = text[::-1]
    revcomp = []
    for let in revtext:
        revcomp.append(compd[let])
    return ''.join(revcomp)

