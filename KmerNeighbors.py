#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 13:26:24 2018

@author: jennifermichaud
"""

def Neighbors(pattern, d):
    """Inputs a DNA pattern as a text string and a number of allowed mismatches, d. 
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