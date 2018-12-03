#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 16:20:30 2018

@author: jennifermichaud
"""
import math
import numpy as np

    

def OutputLCSBacktrack(seq1,seq2, backtrack):
    """Inputs two sequences as strings and a backtrack matrix.
    Outputs two strings that represent the alignment between seq1 and seq2 
    respectively incorporating longest common substring between them using
    bactrack matrix.
    """
    i = len(seq1)
    j= len(seq2)
    if i==0 or j==0:
        return
    subseq1 = []
    subseq2 = []
    while i > 0 and j > 0:
        if backtrack[i,j]== "↘": #matches add the base to each subsequence and move to diagonal up/left position
            subseq1.append(seq1[i-1])
            subseq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif backtrack[i,j]== "↓": #append subseq1 with seq1 base, '-' to subseq2 and move to up position
            subseq1.append(seq1[i-1])
            subseq2.append("-")
            i -= 1
        elif backtrack[i,j]== "→": #append subseq2 with seq2 base, '-' to subseq1 and move to left position
            subseq1.append("-")
            subseq2.append(seq2[j-1])
            j -= 1 
    if i == 1:
        subseq1.append(seq1[i-1])
        subseq2.append("-")
    elif j == 1:
        subseq1.append("-")
        subseq2.append(seq2[j-1])        
    subseq1 = subseq1[::-1]
    subseq2 = subseq2[::-1]
    return subseq1, subseq2


def LocalLCSBacktrack(seq1,seq2, backtrack, max1, max2):
    """Inputs two sequences as strings, a backtrack matrix, and the position of
    the max score.
    Outputs alignment incorporating longest common substring between them using 
    bactrack matrix.
    """
    i = max1
    j = max2
    if i==0 or j==0:
        return
    subseq1 = []
    subseq2 = []
    stop = 0
    while stop == 0:
        if backtrack[i,j]== "↘":
            subseq1.append(seq1[i-1])
            subseq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif backtrack[i,j]== "↓":
            subseq1.append(seq1[i-1])
            subseq2.append("-")
            i -= 1
        elif backtrack[i,j]== "→":
            subseq1.append("-")
            subseq2.append(seq2[j-1])
            j -= 1 
        elif backtrack[i,j]== ">":
            stop = 1
    subseq1, subseq2 = subseq1[::-1], subseq2[::-1]
    return subseq1, subseq2


def LocalAlignment(seq1, seq2, scoringmatrix, sigma = -5):
    """Inputs two strings of peptide sequences (seq1, seq2), a scoringmatrix as a nested
    dictionary, and a default indel penalty. The scoring matrix dictionary allows 
    for specific penalties/ rewards between matches and mismatches.
    Outputs the maximum alignment score of these strings as well as two strings 
    that represent the alignment between seq1 and seq2 respectively incorporating 
    longest common substring between them using a bactrack matrix. This local 
    aligment algorithm allows for low scoring pairings to be skipped allowing for
    higher score between shorter subregions of the sequences.
    """
    v = len(seq1)
    w= len(seq2)
    maxscore = 0
    backtrack = np.empty([v+1,w+1], dtype = str ) #initialize bactrack matrix
    pathlengthm = np.zeros((v+1, w+1))#initialize a matrix to calc the longest path lenths
    for n in range(1,v+1):#initalize 1st row and 1st column with indel penalty score as they represent insertions and deletions
        pathlengthm[n,0] = sigma + pathlengthm[n-1,0]
    for m in range(1,w+1):
        pathlengthm[0,m] = sigma + pathlengthm[0,m-1]
    for i in range(1,v+1): # for remainder of matrix evaluate max score from precursor nodes.  
        for j in range(1,w+1):
            pathlengthm[i,j] = max(0,pathlengthm[i-1,j] + sigma, pathlengthm[i,j-1] + sigma, pathlengthm[i-1,j-1]+scoringmatrix[seq1[i-1]][seq2[j-1]])
            if  pathlengthm[i,j] == pathlengthm[i - 1,j]+ sigma:
                backtrack[i,j] = "↓"
            elif pathlengthm[i,j] == pathlengthm[i,j-1]+ sigma:
                backtrack[i,j] = "→"
            elif pathlengthm[i,j] == pathlengthm[i-1,j-1]+scoringmatrix[seq1[i-1]][seq2[j-1]]:
                backtrack[i,j] = "↘"
            elif pathlengthm[i,j] == 0:  #indicates that the max score is low so should be skipped in alignment,  makes for a "shortcut"
                backtrack[i,j] = ">"
            if pathlengthm[i,j] > maxscore:
                maxscore = pathlengthm[i,j]
                max1, max2 = i, j         
    align1, align2 = LocalLCSBacktrack(seq1,seq2, backtrack, max1, max2)
    return str(int(maxscore)) +"\n"+ ''.join(align1) +"\n" + ''.join(align2)

##TEST DATA
    
#PAM250 = {'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3}, 'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2, 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3, 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0}, 'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2, 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5}, 'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1, 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7}, 'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2, 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1}, 'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2, 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0, 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0}, 'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2, 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1, 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4}, 'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2, 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2, 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2}, 'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2, 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3, 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1}, 'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2, 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2}, 'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2, 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0, 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4}, 'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2, 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6, 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5}, 'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1, 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3}, 'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2, 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0, 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4}, 'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3}, 'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5, 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6, 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0}, 'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1, 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2}, 'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1, 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}}

#seq1 = "MEANLY"
#seq2 = "PENALTY"

#print(LocalAlignment(seq1, seq2, PAM250))


def LocalForFitting(seq1, seq2, backtrack, max1, max2):
    """Inputs two sequences as strings, a backtrack matrix, and the position of
    the max score (two coordinates representing row and column.
    Outputs  two strings representing the aligment of seq1 and seq2 respectively
    incorporating best alignment of full seq2 against portion
    of seq1."""
    i = max1
    j = max2
    if i==0 or j==0:
        return
    subseq1 = []
    subseq2 = []
    while j>0:
        if backtrack[i,j]== "↘":
            subseq1.append(seq1[i-1])
            subseq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif backtrack[i,j]== "↓":
            subseq1.append(seq1[i-1])
            subseq2.append("-")
            i -= 1
        elif backtrack[i,j]== "→":
            subseq1.append("-")
            subseq2.append(seq2[j-1])
            j -= 1 
        elif backtrack[i,j]== ">":
            subseq1.append(seq1[i-1])
            subseq2.append(seq2[j-1])
            i -= 1
            j -= 1
    subseq1, subseq2 = subseq1[::-1], subseq2[::-1]
    return subseq1, subseq2


def FittingAlignment(seq1, seq2):
    """Inputs two strings.  Finds a fitting alignment of a shorter string(seq2)
    against a longer string(seq1). Matches get a score of +1; mismatches and 
    indel penalites are -1.
    Outputs the maximum alignment score of seq1/seq2 as well as two strings
    representing the aligment of seq1 and seq2 respectively.
    """
    v = len(seq1)
    w = len(seq2)
    maxscore = -math.inf
    backtrack = np.empty([v+1,w+1], dtype = str ) #initialize bactrack matrix
    pathlengthm = np.zeros((v+1, w+1))#initialize a matrix to calc the longest path lenths
    for m in range(1,w+1):#initalize 1st row with indel penalty score. To give "free taxi ride" to begining of second seq leave first columns as zeroes.
        pathlengthm[0,m] = -1 + pathlengthm[0,m-1]
    for i in range(1,v+1): # for remainder of matrix evaluate max score from precursor nodes.  
        for j in range(1,w+1):
            match = -1
            if seq1[i-1] == seq2[j-1]:
                match = 1
            pathlengthm[i,j] = max(pathlengthm[i-1,j] -1, pathlengthm[i,j-1] -1, pathlengthm[i-1,j-1]+match)
            if  pathlengthm[i,j] == pathlengthm[i - 1,j]-1:
                backtrack[i,j] = "↓"
            elif pathlengthm[i,j] == pathlengthm[i,j-1]-1:
                backtrack[i,j] = "→"
            elif pathlengthm[i,j] == pathlengthm[i-1,j-1]+match:
                backtrack[i,j] = "↘"
            elif pathlengthm[i,j] == 0:
                backtrack[i,j] = ">"
    for k in range(0,v+1): #select maxscore and its coordinates from last column to be sure all seq2 included
        if pathlengthm[k,-1] > maxscore:
            maxscore = pathlengthm[k,w]
            max1, max2 = k, w                  
    align1, align2 = LocalForFitting(seq1,seq2, backtrack, max1, max2) #create alignment by backtracking bactracking matrix  
    return str(int(maxscore)) +"\n"+ ''.join(align1) +"\n" + ''.join(align2)

##TEST DATA
    
#seq1 = "GTAGGCTTAAGGTTA"
#seq2 = "TAGATA"

#seq1 = "GTTGGATTACGAATCGATATCTGTTTG"
#seq2 = "ACGTCG"    

#print(FittingAlignment(seq1, seq2)) 


def OverlapAlignment(seq1, seq2):
    """Inputs two sequences and generates the highest scoring overlap alignment 
    between the end of seq1 and beginning of seq2
    using alteration of global/local alignment. For scoring, matches count +1. 
    Indels and mismatches are -2.
    Outputs two strings that reperesent the overlap alignment between seq1 and
    seq2 respectively.
    In the longest path matrix, the first column is set to zero so that a "free 
    taxi ride" is provided to begining of seq2.  Similarly the first row  is set
    to zero to give a "free taxi ride" from end of prefix set maxscore and 
    coordinates according to best max score of last row which represents the 
    end of the prefix of seq2.
    """
    v = len(seq1)
    w= len(seq2)
    maxscore = -math.inf
    backtrack = np.empty([v+1,w+1], dtype = str ) #initialize bactrack matrix
    pathlengthm = np.zeros((v+1, w+1)) #initialize a matrix to calc the longest path lenths
    for i in range(1,v+1):# leave first row and first column as zeroes. for remainder of matrix evaluate max score from precursor nodes. 
        for j in range(1,w+1):
            match = -2
            if seq1[i-1] == seq2[j-1]:
                match = 1
            pathlengthm[i,j] = max(pathlengthm[i-1,j] -2, pathlengthm[i,j-1] -2, pathlengthm[i-1,j-1]+match)
            if  pathlengthm[i,j] == pathlengthm[i - 1,j]-2:
                backtrack[i,j] = "↓"
            elif pathlengthm[i,j] == pathlengthm[i,j-1]-2:
                backtrack[i,j] = "→"
            elif pathlengthm[i,j] == pathlengthm[i-1,j-1]+match:
                backtrack[i,j] = "↘"
    k = w 
    while k >0: #select maxscore and its coordinates from last row to get longest lengh of prefix of seq2
        if pathlengthm[v,k] > maxscore:
            maxscore = pathlengthm[v,k]
            max1, max2 = v, k
        k -= 1      
    align1, align2 = LocalForFitting(seq1,seq2, backtrack, max1, max2)#create alignment by backtracking bactracking matrix  
    outlist = [str(int(maxscore)), ''.join(align1),  ''.join(align2)]
    return outlist[0] +"\n"+ outlist[1]+"\n"+ outlist[2]

##TEST DATA
    
#seq1 = "PAWHEAE"
#seq2 = "HEAGAWGHEE"

#print(OverlapAlignment(seq1, seq2))


def GlobalAlignmentAffine(seq1, seq2, scoringmatrix, opening = -11, extension = -1):
    """Inputs two strings of peptide sequences (seq1, seq2), a scoringmatrix as 
    a nested dictionary, and default penalties for gap opening(sigma) and gap 
    extension(epsilon). The scoring matrix dictionary allows 
    for specific penalties/ rewards between matches and mismatches. This algorithm
    calculates three longest paths using insertions, matches/mismathes, or deletions
    with respect to seq1.  The matches/mismatches lenth matrix uses values from
    the other 2 length matrices to determine their values. The resulting matrices
    are used to generate one backtrack matrix. This approach allows for comparison 
    extension versus introducing insertions/delection (opening) making for stronger
    alignments.
    Outputs the maximum alignment score of these strings as well as two strings 
    that represent the aligment of the seq1 and seq2 respectively. 
    """
    v = len(seq1)
    w = len(seq2)
    backtrack = np.empty([v+1,w+1], dtype = str ) #initialize bactrack matrix
    lower = np.zeros((v+1, w+1)) #initialize 3 matrices to calc the longest path lenths, lower is to represent insertions
    middle = np.zeros((v+1, w+1))#match/mismatches
    upper = np.zeros((v+1, w+1)) #deletions with respect to seq1
    lower[1,0] = opening  #the begining of an alignment starts with an opening
    middle[1,0] = opening 
    upper[1,0] = -math.inf  #the first row of lower and first column of upper are set to -inf as they are not valid paths
    lower[0,1] = -math.inf  
    middle[0,1] = opening 
    upper[0,1] = opening 
    for n in range(2,v+1): #Extend first rows and columns with extension penalty
        lower[n,0] = extension + lower[n-1,0]
        middle[n,0] = extension + middle[n-1,0]
        upper[n,0] = extension + upper[n-1,0]
    for m in range(2,w+1):
        lower[0,m] = extension + lower[0,m-1]
        middle[0,m] = extension + middle[0,m-1]
        upper[0,m] = extension + upper[0,m-1]
    for i in range(1,v+1): # for remainder of matrix evaluate max score from precursor nodes. 
        for j in range(1, w+1):
            lower[i,j] = max(lower[i-1,j] + extension, middle[i-1,j] + opening)
            upper[i,j] = max(upper[i,j-1] + extension, middle[i,j-1] + opening)
            middle[i,j] = max(lower[i,j], upper[i,j], middle[i-1,j-1]+scoringmatrix[seq1[i-1]][seq2[j-1]])
            if middle[i,j] == lower[i,j]:
                backtrack[i,j] = "↓"
            elif middle[i,j] == upper[i,j]:
                backtrack[i,j] = "→"
            elif middle[i,j] ==  middle[i-1,j-1]+scoringmatrix[seq1[i-1]][seq2[j-1]]:
                backtrack[i,j] = "↘"         
    align1, align2 = OutputLCSBacktrack(seq1,seq2, backtrack) #create alignment by backtracking bactracking matrix   
    return str(int(middle[-1][-1])) +"\n"+ ''.join(align1) +"\n"+ ''.join(align2)

##TEST DATA
    
#Blosum62_dict = {'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2}, 'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2}, 'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2}, 'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3}, 'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3}, 'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3}, 'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1}, 'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2}, 'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2}, 'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1}, 'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1}, 'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2}, 'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1}, 'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3}, 'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2}, 'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2}, 'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2}, 'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2}, 'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1}, 'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7}}
    
#seq1 = "PRTEINS"
#seq2 = "PRTWPSEIN"
#    
#print(GlobalAlignmentAffine(seq1, seq2, Blosum62_dict))

#seq1 = "YHFDVPDCWAHRYWVENPQAIAQMEQICFNWFPSMMMKQPHVFKVDHHMSCRWLPIRGKKCSSCCTRMRVRTVWE"
#seq2 = "YHEDVAHEDAIAQMVNTFGFVWQICLNQFPSMMMKIYWIAVLSAHVADRKTWSKHMSCRWLPIISATCARMRVRTVWE"
#    
#print(GlobalAlignmentAffine(seq1, seq2, Blosum62_dict))

#seq1 = "KLQDIFFADFDCACSMPAAKTETPEYPCGPNCIHDQEAARSVSSHRRSETRIECRSNRHERQYWAHGGIPCINPYRRNQEEFIN"
#seq2 = "KKTDSRFADFDCYPWVYMRSPPEYPCGPNNPMVDMFWLIHDQAAARQMFVYVASHRETRIELKSFRHERQYWAHGGIPCINPYRRNQEEFIN"
#
#print(GlobalAlignmentAffine(seq1, seq2, Blosum62_dict))


   