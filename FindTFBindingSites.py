#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 13:32:02 2018

@author: jennifermichaud
"""
"""  It uses Gibb sampling techniques.
"""
import numpy as np
import math
import copy

#FUNCTIONS
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
        rando = np.random.choice((len(dnalist[0])-(k+1)))
        motifs.append(dnalist[i][rando:rando+k])
    bestmotifs = copy.deepcopy(motifs)
    for j in range(N):
        ind = np.random.choice(t) # choose random kmer in list
        del motifs[ind]#remove the random kmer from motifs
        pcount = LaplaceCount(motifs, k)  
        profile = Profile(pcount, motifs, k)
        pos = WeightedRandom(profile, dnalist[ind], k)
        motifs.insert(i, dnalist[ind][pos:pos+k])  #based upon profile weighted random position insert new kmer into motifs
        consensus = Consensus(profile, k)
        motifscore = MotifsScore(motifs, consensus)
        bestmotifscore = MotifsScore(bestmotifs, consensus)
        if motifscore < bestmotifscore:
            bestmotifs = copy.deepcopy(motifs)
        motifs = copy.deepcopy(bestmotifs)
    scoremotifs = [bestmotifs, bestmotifscore, consensus]
    return scoremotifs   # this outputs both the best motifs and the score


def ModGibbsMedianMotifs(dnalist, k, t, N, M):
    """Inputs a list of DNA strings (dnalist), a kmer length (k), and DNA list 
    length (t), a number of iterations (N), and a number of starts (M).
    Ouputs a list of best median motifs for each of the DNA strings in the given 
    list minimizing the hamming distance between them. This efficient algorithm 
    is based upon Gibbs Sampling (Bayesian Inference) combining the advantages 
    of weighted probability and randomness.  This modified version of 
    GibbsMedianMotifs gives the GibbsMedianMotifs() M random
    starts and returns the best scoring collection of median strings, the consensus
    sequence between them, and the score of the best scoring collection.
    """  
    bestscore = math.inf
    bestmotifcollection = []
    scorelist = []
    consensuslist = []
    bestscorepos = 0
    count = 0
    while count < M:
        bestlist = GibbsMedianCore(dnalist, k, t, N)
        bestmotifcollection.append(bestlist[0])
        scorelist.append(bestlist[1])
        consensuslist.append(bestlist[2])
        count += 1
    for i in range(len(scorelist)):
        if scorelist[i] < bestscore:
            bestscore = scorelist[i]
            bestscorepos = i
            bestconsensus = consensuslist[i]
    motifsconsensus = [bestmotifcollection[bestscorepos], bestconsensus, bestscore]
    return motifsconsensus

## Find transcription binding sites from upstream sites of several genes. 

# Load DNA sequences.  In this case, they are a text file '\n' seperated capitalized
#sequences and which are stripped to make into a list of the sequences.    
#This example takes upstream regions 10 genes with expression changes during
# hypoxia from Mycobacterium tuberculosis in effort to find transcription factor 
# binding sites for dormancy survival regulator (DosR), a transcription factor
# related to dormancy. Following example attempts to identify  DosR transcription factor binding sites 
#from 10 upstream sites out of 25 of identified DosR genes from Mycobacterium tuberculosis.  
# Use algorithm with k= 20, 2000 runs N = 200. Other kmer lengths can also be tested.

f = open('/Users/.../DosR.txt', 'r') 
dnatxt = f.read()
f.close()

# put data into list if each DNA fragment is seperated by '\n'
dnalist = [y for y in (x.strip() for x in dnatxt.splitlines()) if y]

# Set algorithm parameters
#  k = motif length, t = number of dna sequences, N = number of iterations of 
# core algorithm optimization, M = number of starts of core algorithm.  t should
# not be altered.
k = 20
N = 1000
M = 30

t = len(dnalist)

#Generate data
motiflist = ModGibbsMedianMotifs(dnalist, k, t, N, M)
print('Median motifs: ')
print('\n'.join(motiflist[0]))  #  Motiflist
print('Consensus sequence: ' + motiflist[1]) #consensus sequence
print('Motif score: ' + str(motiflist[2])) #Motif score

# Below is the output from 3 different runs using different values of N (iterations)
# and M(number of starts). From this data N=1000, M=30 results in a score of 0 
# finding exact sequences with no mismatches from the consensus sequence in all
# tested regions.  From this data 'GTTTCCATCGATTAGGAGGC' very likely captures
# a DosR transcription factor binding site in Mycobacterium tuberculosis


##output k = 20, N=20, M= 20
#Median motifs: 
#CCCTAGCCCTGGCCACGATG
#CCCTAGCCCTGGCCACGATG
#CCCTAGCCCTGGCCACGATG
#CCCTAGCCCTGGCCACGATG
#CTTCGGCCCCACCCACGAGG
#CGCTAACCCTGGCTTCGATG
#CCCCAGCGAAGGAGACGGCG
#CTTCGGCCCCACCCACGAGG
#CGCTAACCCTGGCTTCGATG
#CTTCGGCCCCACCCACGAGG
#Consensus sequence: CCCTAGCCCTGGCCACGATG
#Motif score: 40

##output k = 20, N=1000, M= 25
#Median motifs: 
#TGCGCAGCGACGGAAACTGC
#TGCGCAGCGACGGAAACTGC
#TGCGCAGCGACGGAAACTGC
#TGCGCAGCGACGGAAACTGC
#TGCGCAGCGACGGAAACTGC
#TGCGCAGCGACGGAAACTGC
#TGCGCAGCGACGGAAACTGC
#TGCGCAGCGACGGAAACTGC
#GCCCCAGCGAAGGAGACGGC
#TGCGCAGCGACGGAAACTGC
#Consensus sequence: TGCGCAGCGACGGAAACTGC
#Motif score: 6

##output k = 20, N=1000, M= 30
#Median motifs: 
#GTTTCCATCGATTAGGAGGC
#GTTTCCATCGATTAGGAGGC
#GTTTCCATCGATTAGGAGGC
#GTTTCCATCGATTAGGAGGC
#GTTTCCATCGATTAGGAGGC
#GTTTCCATCGATTAGGAGGC
#GTTTCCATCGATTAGGAGGC
#GTTTCCATCGATTAGGAGGC
#GTTTCCATCGATTAGGAGGC
#GTTTCCATCGATTAGGAGGC
#Consensus sequence: GTTTCCATCGATTAGGAGGC
#Motif score: 0



