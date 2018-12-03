#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 16:46:23 2018

@author: jennifermichaud
"""

def RemoveListDuplicates(listoflists):
    """Inputs a list of list of numbers and outputs a list with duplicates
    removed.
    """
    removedups = []
    for slist in listoflists:
        removedups.append(" ".join(list(map(str, slist))))
    removedups = set(removedups)
    removedups = list(removedups)
    correctedlist = []
    for dlist in removedups:
        temp = list(map(int, dlist.split(' ')))
        correctedlist.append(temp)
    return correctedlist