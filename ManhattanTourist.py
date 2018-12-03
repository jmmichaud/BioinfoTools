#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 19:33:47 2018

@author: jennifermichaud
"""

import numpy as np

def Converttonumpy(returnsepstrings):
    listofstrings = returnsepstrings.split("\n")
    listlist = []
    for lst in listofstrings:
        lst = list(map(int, lst.split(' ')))
        listlist.append(lst)
    return np.array(listlist)


def ManhattanTourist(n, m, down, right):
    """ The Mattahattan Tourist problem describes an optimal one-way path
    (with no overlap or backtracking) that a tourist in Manhattan should take to 
    visit as many tourist attractions as possible. Such calculations are useful
    to a variety of computational problems. This function calculates the 
    length of the longest path possible when the dimensions are given ("the street 
    grid") and the number of attractions on each vertical and horizontal block. 
    The function inputs grid dimensions n (vertical) and m (horizontal), as well 
    as a down(n x (m+1)) and right(m x (n+1)) matrices as numpy arrays.  
    Creates a directed acyclic graph represented by a matrix of 
    n+1(rows) by m+1(columns) of nodes denoting longest path lengths 
    using the weight values in the down and right matrices. 
    Outputs the length of longest path.
    """
    pathlengthm = np.zeros((n+1, m+1))#initialize a matrix to calc the longest path lengths
    currlength = 0
    for i in range(1,n+1): #fill in path lengths for first column using down matrix
        currlength += down[i-1,0] 
        pathlengthm[i, 0] = currlength  
    currlength = 0
    for j in range(1, m+1): #fill in path lengths for first row using right matrix 
        currlength += right[0, j-1] 
        pathlengthm[0, j] = currlength
    for x in range(1,n+1):# for remainder of matrix evaluate left and up calcs to determine max
        for y in range(1,m+1):# each up or left path includes previous value and weight
            uppath = pathlengthm[x-1,y] + down[x-1,y]
            leftpath = pathlengthm[x,y-1] + right[x,y-1]
            maxpath = max(uppath, leftpath)
            pathlengthm[x,y] = maxpath
    return int(pathlengthm[n,m])
            
        
##TEST DATA 
    
#n = 4     
#m = 4    
#down_in = """1 0 2 4 3
#4 6 5 2 1
#4 4 5 2 1
#5 6 8 5 3"""
#right_in = """3 2 4 0
#3 2 4 2
#0 7 3 3
#3 3 0 2
#1 3 2 2"""
#    
##exp result = 34
    
#n = 12     
#m = 5   
#down_in = """3 0 2 2 0 4
#0 2 2 3 1 3
#0 4 4 1 3 3
#4 3 2 4 2 1
#3 0 2 1 0 4
#3 1 2 1 4 2
#1 2 4 2 1 1
#1 1 1 0 1 4
#0 0 2 4 0 4
#4 4 1 3 4 2
#4 2 2 0 2 1
#0 0 1 4 1 2"""
#right_in = """0 4 2 1 1
#1 0 3 3 3
#2 3 0 0 1
#2 0 3 2 4
#4 3 3 0 2
#4 4 0 1 2
#4 2 3 0 0
#3 0 2 1 4
#1 2 3 4 3
#0 2 3 3 3
#3 3 2 4 4
#1 2 1 1 4
#2 3 4 3 4"""

##expected result = 47
#
#
#down = Converttonumpy(down_in)    
#      
#right = Converttonumpy(right_in)  
#
#print(ManhattanTourist(n, m, down, right))