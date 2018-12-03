#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 16:47:26 2018

@author: jennifermichaud
"""
import math
import numpy as np
import random
import copy

def NumberToPattern(kmerdex, k): 
    """Inputs a integer dnadex that is a numerical representation of a DNA 
    sequence, and an interger kmer length k.
    Outputs the DNA sequence encodes by the quarternary code).
    """
    pos = k
    divfour = int(kmerdex / 4)
    remfour = int(kmerdex % 4)
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
    
def HammingDistance(text1,text2):
    """Computes the hamming distance between
    two equal length dna strings.  Each mismatch 
    is a score of 1.
    """
    hamming = 0
    for i in range(len(text2)):
        if text1[i] != text2[i]:
            hamming += 1
    return hamming

def DistanceBetweenPatternAndStrings(pattern, dnalist):
    """Inputs a string pattern and a list of strings. Outputs an integer score
    that is the culumaltive sum of the minimum achievable hamming distance between 
    the pattern and each string.
    """
    k = len(pattern)
    distance = 0
    for dstring in dnalist:
        hammingdistance = math.inf
        for i in range(len(dstring)-k+1):
            test = dstring[i:i+k]
            testham = HammingDistance(pattern, test)
            if hammingdistance > testham:
                hammingdistance = testham
        distance += hammingdistance    
    return distance


def BruteMedianMotifs(dnalist, k):
    """Inputs a list of DNA strings (dnalist) and a kmer or motif length (k).
    Outputs a list of motifs that represent the lowest achievable distance between
    the motif and the DNA fragments in the list. This brute force algorithm
    generates all possible motifs of length k and then tests those motifs against
    the DNA list."""
    distance = math.inf
    medianlist = []
    for i in range((4**k)-1):
        pattern = NumberToPattern(i, k) #generates all possible k-length kmers
        patdist = DistanceBetweenPatternAndStrings(pattern, dnalist) #calcuate score of motif
        if distance > patdist: #if distance is lower than previous lowest score adds to list
            distance = patdist
            medianlist = []
            medianlist.append(pattern)
        elif distance == patdist: #if distance matches lowest score adds to list
            medianlist.append(pattern)
            
    return medianlist


#TEST DATA
#k_in = 6
#inputdna = """TAAGACCATTTTCTAACTTCGGGATTTCCGAAACGGAATTTC
#CCCGCTCCGAACCGTTATTAAGACTGCAGGAGACATTGTGGT
#GAAGCTTGCATGTGACAAGCGAGTCGACAGGGACATTAAGAC
#TCAATCCGGGCGTGCTGACCTAGCGTAGGTAGACACTAACAC
#CCGGGTCCGAGCAAGCCAGCCCTAAGTTTCTAACACGGCTCA
#CACGTTTCCCGGTAATACATCGAACGCCTATTCCCGGGTGCG
#CAGGCGGCCCCACACAGAATAGAGGCTATACAGATATAACAC
#AGGCTCGGCCGTTAAAACTTAAGCTGACATGACATTGACCCA
#ACGGAACCTCTTTCTCATTACTATTATGTAGCTGCGTAATAC
#GGACTTAGAGTAGTGTGCTAAAACTGGCCCTCCTTGGAGATG"""
#dnalist_in = inputdna.split()   
##Correct_out = TAACAC TAAGAC
#
#print(' '.join(BruteMedianMotifs(dnalist_in, k_in)))
  

def MostProbableKmer(dna, k, profile):
    """ Inputs a DNA string (dna), a kmer length (k), and a probality matrix (profile)
    that conveys the probability of finding a particular DNA base at a given 
    position.  The probability matrix is a numpy array where row 0 through 3 
    represent the probablity of finding an A, C, G, or T, respectively, and 
    each column corresponds to a position in the k-lenthg motif.   
    Outputs the most probable motif or kmer using the probability matrix.
    """
    prob_list = []
    prob_calc = 1
    for i in range(len(dna)-k +1):#for length of dna create test kmer
        test = dna[i:i+k]
        for x in range(k): #calculate probability score for each test kmer
            if test[x] == 'A':
                prob_calc *= profile[0, x]
            elif test[x] == 'C':
                prob_calc *= profile[1, x]
            elif test[x] == 'G':
                prob_calc *= profile[2, x]
            elif test[x] == 'T':
                prob_calc *= profile[3, x]
        prob_list.append(prob_calc)
        prob_calc = 1
    highestprob = max(prob_list) #Find highest probability score
    position = prob_list.index(highestprob) #Find position of highest score
    prob_kmer = dna[position:position +k] #Using poistion generate highest probable kmer
    return prob_kmer


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
    for n in range(4): #generate probability for each base 0-4 represents each base/row
        for l in range(k): #for each column/position 0-k
            prob_ACGT[n, l] = np.divide(count_ACGT[n,l],len(motif_list)) #prob is count of base divided by length of DNA list
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
    for i in range(k): #for each column postition find highest probility between the rows and the postion of the highest score is used to determine highest probable base.
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


def MotifsScore(dnalist, consensus):
    """Inputs a list of DNA strings (dnalist) and a string consensus sequence (consensus)
    between them.
    Outputs the culmulative distance between the consensus string and each DNA
    strand.
    """
    score = 0
    for motif in dnalist:
        score += HammingDistance(motif, consensus)
    return score
                    
        
def GreedyMedianMotifs(dna, k, t):
    """Inputs a list of DNA strings, a kmer length (k), and a length of dnalist (t).
    Outputs the best median strings between the DNA strands.  Utilizing matrices
    that generate the counts and probability of each base at a given position and
    a generated consensus sequences from the generated probility matrices, finds 
    motifs that minimize a score based upon hamming distance between a 
    computed consensus sequence and DNA strands.  
    """
    bestmotifs = []
    motifs = []
    bestscore = math.inf
    for i in range(len(dna[0])-k+1): #choose base motif from 1st string of DNA and begins motif list
        motif = dna[0][i:i+k]
        motifs.append(motif)
        for j in range(1,t):#make profile(probability matrix) of current list use to find most probable kmer in current dna string and adds it to motif list
            count = LaplaceCount(motifs, k)
            profile = Profile(count, motifs, k)
            mpkmer = MostProbableKmer(dna[j], k, profile)
            motifs.append(mpkmer)
        consensus = Consensus(profile, k)
        score = MotifsScore(motifs, consensus) #score motif list using profile to make a consensus seq
        if score < bestscore: #supplants current best motif is the score is lower than other motif lists
            bestmotifs = motifs
            bestscore = score
        motifs = []
    return bestmotifs            
            

#TEST DATA  
#k_in = 3
#t_in = 5
#inputdna = """GGCGTTCAGGCA
#AAGAATCAGTCA
#CAAGGAGTTCGC
#CACGTCAATCAC
#CAATAATATTCG"""
##Correct_out = TCA
##TCA
##TCG
##TCA
##TCG
#dnalist_in = inputdna.split("\n")
#
#print('\n'.join(GreedyMedianMotifs(dnalist_in, k_in, t_in)))


def ProfileMostProbableMotifs(profile, dnalist, k, t):
    """Inputs a probability matrix (profile), a list of DNA strings (dnalist),
    a kmer length (k), and DNA list length(t). 
    Outputs a list of most probable kmers from each strand of DNA given the probability
    profile.
    """
    motifs = []
    for j in range(0,t):#generate most probable kmer list using current profile
            mpkmer = MostProbableKmer(dnalist[j], k, profile)
            motifs.append(mpkmer)
    return motifs


def RandomizedMotifSearch(dnalist, k, t):
    """Inputs a probability matrix (profile), a list of DNA strings (dnalist),
    a kmer length (k), and DNA list length(t). 
    Outputs a list of the best median motifs between the DNA strands and the score
    of motifs.  The motifs minimize a score based upon hamming distance, a random 
    selection of inital motifs, and an optimization of a motif list and profile matrix. 
    A randomselection of kmers from each of the DNA strings is used to generate 
    a probability matrix. The probability matrix is used to generate new candidate 
    best motifs until the score is no longer optimized.
    """
    motifs = []
    for i in range(t): #generate list of random kmers for each string in dna
        rando = random.randint(0,(len(dnalist[0])-(k+1)))
        motifs.append(dnalist[i][rando:rando+k])
    bestmotifs = motifs
    while True: #until score is no longer optimized
        pcount = LaplaceCount(motifs, k)
        profile = Profile(pcount, motifs, k)
        consensus = Consensus(profile, k)
        motifs = ProfileMostProbableMotifs(profile, dnalist, k, t) #motifs are populated as dictated by probility matrix
        motifscore = MotifsScore(motifs, consensus)
        bestmotifscore = MotifsScore(bestmotifs, consensus)
        if motifscore < bestmotifscore: #if new motif list minimizes score it becomes best motifs
            bestmotifs = motifs
        else:
            scoremotifs = [bestmotifs, motifscore]
            return scoremotifs   # this outputs both the best motifs and the score

def RandomizedMedianMotifs(dna, k, t, M=1000): 
    """Inputs a probability matrix (profile), a list of DNA strings (dnalist),
    a kmer length (k), and DNA list length(t). 
    Outputs the best median motifs between the DNA strands that minimize a score
    based upon hamming distance, a random selection of inital motifs, and an 
    optimization of a motif list and profile matrix. A random
    selection of kmers from each of the DNA strings is used to generate a probability
    matrix. The probability matrix is used to generate new candidate best motifs
    until the score is no longer optimized. This is run a M times (1000 default) 
    and outputs the lowest scoring motif list.
    """
    bestscore = math.inf
    bestmotifcollection = []
    scorelist = []
    bestscorepos = 0
    count = 0
    while count < M:
        bestlist = RandomizedMotifSearch(dna, k, t)
        bestmotifcollection.append(bestlist[0])
        scorelist.append(bestlist[1])
        count += 1
    for i in range(len(scorelist)):
        if scorelist[i] < bestscore:
            bestscore = scorelist[i]
            bestscorepos = i
    return bestmotifcollection[bestscorepos]
        

##TEST DATA    
    
#k_in = 8
#t_in = 5
#inputdna = """CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
#GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
#TAGTACCGAGACCGAAAGAAGTATACAGGCGT
#TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
#AATCCACCAGCTCCACGTGCAATGTTGGCCTA"""
##corrext out=
##TCTCGGGG
##CCAAGGTG
##TACAGGCG
##TTCAGGTG
##TCCACGTG
#dnalist_in = inputdna.split("\n")
#
#
#print('\n'.join(RandomizedMedianMotifs(dnalist_in, k_in, t_in)))   


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
    