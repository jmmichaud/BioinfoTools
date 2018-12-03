#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 13:40:39 2018

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


def HammingDistance(text1, text2):
    """Computes the hamming distance between
    two equal length dna strings.  Each mismatch 
    is a score of 1"""
    hamming = 0
    for i in range(len(text2)):
        if text1[i] != text2[i]:
            hamming += 1
    return hamming


def MotifsScore(motifs, consensus):
    """Inputs a list of DNA strings (dnalist) and a string consensus sequence (consensus)
    between them.
    Outputs the culmulative distance between the consensus string and each DNA
    strand.
    """   
    score = 0
    for motif in motifs:
        score += HammingDistance(motif, consensus)
    return score


def WeightedRandom(profile, selecteddna, k): 
    """Inputs a probability matrix (profile), a dna string (selecteddna), and a 
    kmer length (k). 
    Outputs an integer that is the index of the staring position of the kmer
    within the dna string randomly selected but weighted by probalities generated
    by calculating the probability of base of the kmer staring at that position of 
    the dna string given the probability matrix.  Randomness avoids forcing
    """
    prob = []
    probcalc = 1
    for n in range(len(selecteddna)-k+1): #test each kmer along the dna string
        test = selecteddna[n:n+k]
        for m in range(k): #calculate probility for each kmer and add to prob list
            if test[m] == 'A':
                probcalc *= profile[0, m]
            elif test[m] == 'C':
                probcalc *= profile[1, m]  
            elif test[m] == 'G':
                probcalc *= profile[2, m]   
            elif test[m] == 'T':
                probcalc *= profile[3, m]     
        prob.append(probcalc)
        probcalc = 1
    prob = np.array(prob)
    if sum(prob) != 1: #normalize probabilities to sum to 1
        prob = prob/np.sum(prob)
    return np.random.choice(len(prob), p = prob)


def GibbsMedianCore(dnalist, k, t, N): 
    """Inputs a list of DNA strings (dnalist), a kmer length (k), and DNA list 
    length (t), and a number of iterations (N).
    Ouputs a list of best median motifs for each of the DNA strings in the given 
    list and the score of the motifs. This efficient algorithm is based upon 
    Gibbs Sampling (Bayesian Inference) combining the advantages of weighted 
    probability and randomness.
    """
    motifs = []
    for i in range(t): #generate list of random kmers for each string in dna
        rando = np.random.choice((len(dnalist[0])-k+1))
        motifs.append(dnalist[i][rando:rando+k])
    bestmotifs = copy.deepcopy(motifs)
    for j in range(N):
        ind = np.random.choice(t) # choose index of random dna in list
        del motifs[ind]#remove the kmer of selected random index # from motifs
        pcount = LaplaceCount(motifs, k)  
        profile = Profile(pcount, motifs, k)
        pos = WeightedRandom(profile, dnalist[ind], k)
        motifs.insert(ind, dnalist[ind][pos:pos+k])  #based upon profile weighted random position insert new kmer into motifs
        consensus = Consensus(profile, k)
        motifscore = MotifsScore(motifs, consensus)
        bestmotifscore = MotifsScore(bestmotifs, consensus)
        if motifscore < bestmotifscore:
            bestmotifs = copy.deepcopy(motifs)
        motifs = copy.deepcopy(bestmotifs)
    scoremotifs = [bestmotifs, bestmotifscore]
    return scoremotifs   # outputs both the best motifs and the score


def GibbsMedianMotifs(dnalist, k, t, N): 
    """Inputs a list of DNA strings (dnalist), a kmer length (k), and DNA list 
    length (t), and a number of iterations (N).
    Ouputs a list of best median motifs for each of the DNA strings in the given 
    list minimizing the hamming distance between them. This efficient algorithm 
    is based upon Gibbs Sampling (Bayesian Inference) combining the advantages 
    of weighted probability and randomness.  It gives the GibbsMedianMotifs() 25 random
    starts and returns the best scoring collection of median strings.
    """    
    bestscore = math.inf
    bestmotifcollection = []
    scorelist = []
    bestscorepos = 0
    count = 0
    while count < 25:
        bestlist = GibbsMedianCore(dnalist, k, t, N) #Gives motif list and score
        bestmotifcollection.append(bestlist[0])
        scorelist.append(bestlist[1])
        count += 1
    for i in range(len(scorelist)): #select lowest scoring list
        if scorelist[i] < bestscore:
            bestscore = scorelist[i]
            bestscorepos = i
    return bestmotifcollection[bestscorepos]


#TEST DATA
#k_in = 8
#t_in = 5
#N_in = 100
#inputdna = """CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA
#GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
#TAGTACCGAGACCGAAAGAAGTATACAGGCGT
#TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
#AATCCACCAGCTCCACGTGCAATGTTGGCCTA"""
##Correct_out = TCTCGGGG
##CCAAGGTG
##TACAGGCG
##TTCAGGTG
##TCCACGTG
#
#dnalist_in = inputdna.split("\n")
#print(' '.join(GibbsMedianMotifs(dnalist_in, k_in, t_in, N_in)))


#NOTE
#Algorithm can be further optimized by increasing or decreasing either the number of starts
#or N, the number of iterations.  From testing it seems increasing the number
#of starts has the most impact to accuracy.
