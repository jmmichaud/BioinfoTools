#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 13:41:22 2018

@author: jennifermichaud
"""

import copy
import collections


def GenerateCycloSpectra(peptidemasses):
    """Inputs a list of integer masses (min 2 fragments) of a cyclic peptide.
    Outputs the expanded spectra (sums of adjacent fragments) as a list of masses.
    """
    expandedspectra = []
    L = len(peptidemasses) 
    for peptide in peptidemasses:
        expandedspectra.append(peptide)
    expandedspectra.append(sum(peptidemasses))
    i = 2
    while i < L:
        for n in range(L):
            mass = peptidemasses[n:n+i]
            if len(mass) < i:
                mass = peptidemasses[n:]
                frontneeded = i - len(mass)
                mass.extend(peptidemasses[:frontneeded])
            masssum = sum(mass)              
            expandedspectra.append(masssum)
        i += 1
    expandedspectra.sort()   
    return expandedspectra


def GenerateLinearSpectra(peptidemasses):
    """Inputs a list of integer masses (min 2 fragments) of a linear peptide.
    Outpus the expanded spectra (sums of adjacent fragments) as a list of masses"""
    expandedspectra = []
    L = len(peptidemasses) 
    for peptide in peptidemasses:
        expandedspectra.append(peptide)
    expandedspectra.append(sum(peptidemasses))
    i = 2
    while i < L:
        for n in range(L):
            mass = peptidemasses[n:n+i]
            if len(mass) == i:
                masssum = sum(mass)              
                expandedspectra.append(masssum)
        i += 1
    expandedspectra.sort()      
    return expandedspectra


def PeptideCycloScoring(peptideseq, spectrum):
    """Inputs a cyclic peptide (peptideseq) as a list of integers and mass 
    spectra (spectrum) as a list of integers.  
    Outputs the score of the peptide spectrum generated from peptideseq against 
    the inputted spectrum.
    Each match of the peptide generated spectrum to the inputted spectrum yields
    one point. 
    """
    if peptideseq == []:
        return 0
    else:
        peptidemasses = GenerateCycloSpectra(peptideseq)
        score = 0
        specdup = copy.deepcopy(spectrum)
        for spec in peptidemasses:
            if spec in specdup:
                score += 1
                specdup.remove(spec) #remove mass from spectrum to avoid duplicates
        return score


def PeptideLinearScoring(peptideseq, spectrum):
    """Inputs a linear peptide as a list of integers (peptideseq)
    and mass spectra as a list of integers (spectrum).  
    Outputs the score of the peptide generated spectrum against the 
    inputted spectrum.
    Each match of the peptide generated spectrum to the inputted spectrum yields
    one point. 
    """
    if peptideseq == []:
        return 0
    else:
        peptidemasses = GenerateLinearSpectra(peptideseq)
        score = 0
        specdup = copy.deepcopy(spectrum)
        for spec in peptidemasses:
            if spec in specdup:
                score += 1
                specdup.remove(spec) #remove mass from spectrum to avoid duplicates
        return score
       
        
def TrimPeptideLeaderboard(leaderboard, spectrum, N):
    """Inputs a list of lists of candidate peptides for an unknown peptide
    as a leaderboard for the spectrum of the unknown peptide. Outputs candidate 
    peptides with N number of top scores.
    """  
    toppeptides = []
    scoredict = {}
    scorelist = []
    scoreset = []
    for peptide in leaderboard:
        peptidejoin = ' '.join(list(map(str, peptide)))
        peptidescore = PeptideLinearScoring(peptide, spectrum)
        scoredict[peptidejoin] = peptidescore 
        scorelist.append(peptidescore)
        scoreset.append(peptidescore)
    scoredict = collections.OrderedDict(scoredict)
    scoredict = sorted(scoredict.items(), key=lambda t: t[1], reverse = True)
    scorelist.sort(reverse= True)
    scoreset = list(set(scoreset))
    scoreset.sort(reverse= True)
    countdict = {}
    for sco in scorelist:
        if sco in countdict:
            countdict[sco] += 1
        else:
            countdict[sco] = 1       
    Z = 0
    C = 0
    for sc in scoreset:
        Z += countdict[sc]
        if Z >= N:
            C = Z
            break
    if C != Z:
        C = Z
    L = 0       
    for pep in scoredict:
        if L < C:
            toppeptide = pep[0]
            if ' ' in toppeptide:
                toppeptide = list(map(int, toppeptide.split(" ")))
            else:
                toppeptide = [int(toppeptide)]
            toppeptides.append(toppeptide)
            L += 1
    return toppeptides        

    
def CorrectedAAList(spectrum, M):
    """Inputs a mass spectrum as a list of integers and an integer (M).
    Generates a list of masses that represent the likely ammino acids present in
    the spectra.  Uses all non-zero products when elements of spectrum are 
    subtracted from one another and a dictionary to count the number of times a 
    mass appears in the amended mass list.  
    Ouputs a list of  the top M masses (with ties) between 57-200 that appear a 
    minimum of M times in the convoluted mass list. 
    """
    spectrum.sort(reverse = True)
    convoluted = [] #generate convoluted mass list
    for i in range(0,len(spectrum)-1):
        for j in range(i+1,len(spectrum)):
            newmass = spectrum[i] - spectrum[j]
            if newmass != 0:
                convoluted.append(newmass)
    convoluted.sort() 
    cl = copy.deepcopy(convoluted) # Remove masses not between 57 - 200
    for mass in cl:
        if mass < 57 or mass > 200:
            convoluted.remove(mass)
    convolutedcountd = {} #create dictionary, convolutedcountd to count number of each score
    for num in convoluted:
        if num in convolutedcountd:
            convolutedcountd[num] += 1
        else:
            convolutedcountd[num] = 1 
    scorelist = []
    for key in convolutedcountd:
        scorelist.append(convolutedcountd[key])
    scoredict = collections.OrderedDict(convolutedcountd) #Use set of masses and countdict to determine the number,C, of topscores to take
    scoredict = sorted(scoredict.items(), key=lambda t: t[1], reverse = True)
    scorelist.sort(reverse= True)
    scoreset = list(set(scorelist))
    scoreset.sort(reverse= True)
    countdict = {}  #create another dictionary 
    for sco in scorelist:
        if sco in countdict:
            countdict[sco] += 1
        else:
            countdict[sco] = 1    
    Z = 0
    C = 0
    for sc in scoreset:#Find number of items you need from dictionary (C) for top N scores including ties
        Z += countdict[sc]
        if Z >= M:
            C = Z
            break
    if C != Z:
        C = Z     
    L = 0 
    AAmasslist = []      
    for AAmass in scoredict:  #generate masses 
        if L < C:
            AAmasslist.append(AAmass[0])
            L += 1       
    AAmasslist.sort()
    return AAmasslist

            
def CyclopeptideSequencing(spectrum, M, N):
    """Inputs two integers (M and N) and a mass spectrum of an unknown 
    cyclopeptide as list of integer masses. 
    Outputs all lead peptides of a maximal score using N as a cutoff for the Trim
    function including ties as a list of lists where  each peptide is a list 
    of integer masses. 
    Uses an AA list generated from processing of the 
    spectrum and returning M topscoring masses between 57 - 200 including ties.
    Allows for non-standard amino acids and takes advantage of both using a leader 
    board to keep N top scoring peptides and scoring algorithsm that compare 
    the spectrum to each generated spectra of peptides from the leaderboard to 
    generate a leader peptide list where each peptide is maximal achieved score.  
    """  
    leaderboard = []
    leaderpeptides = []
    parentmass = max(spectrum)
    looper = 0
    AAmasslist = CorrectedAAList(spectrum, M)
    if 0 in spectrum: #remove 0 from spectrum
        spectrum.remove(0)
    monocount = 0
    monopeptides = []
    for spec in spectrum: #Build AA list
        if spec <= 200:
            monocount += 1
            leaderboard.append([spec])
            monopeptides.append(spec)
    for mass in AAmasslist: #Start peptide list
        if mass not in monopeptides:
            leaderboard.append([mass])
#    b=0  # a counter that prints to keep track of loop progress
    while leaderboard != []: #generate strings consistent with spectra, add all possibilitiess from monopeptide/AA list and then remove combinations whose mass isn't consistent with the spectra
#        print(b)
#        b += 1
        duppeps = copy.deepcopy(leaderboard)
        leaderboard = []
        for AA in AAmasslist: # generate extended peptides
            for i in range(len(duppeps)):
                templist = copy.deepcopy(duppeps[i])
                templist.append(AA)
                leaderboard.append(templist)
        leadups = copy.deepcopy(leaderboard)
        if len(leadups[0]) >= monocount:
            looper += 1
        for dup in leadups:
            if sum(dup) == parentmass:
                if leaderpeptides == []:
                    leaderpeptides.append(copy.deepcopy(dup)) 
                dupscore =  PeptideCycloScoring(dup, spectrum)
                lpscore =  PeptideCycloScoring(leaderpeptides[0], spectrum)
                if dupscore > lpscore:
                    leaderpeptides = []
                    leaderpeptides.append(copy.deepcopy(dup))
                elif dupscore == lpscore:
                    leaderpeptides.append(copy.deepcopy(dup))
            elif sum(dup) > parentmass: #this will result in the leader board being eventually empty.
                    leaderboard.remove(dup)
        leaderboard = TrimPeptideLeaderboard(leaderboard, spectrum, N) # Trim list.  Keep values of top N scores
    return leaderpeptides

def ConvertPeptideMassToAA(listofpeptidelists):
    """Inputslist of peptide list of list where each sublist is a list of interger
    masses representing a peptide. 
    Outputs a list of peptides using one letter aminoacid code
    of corresponded peptides in AA code. Nonstandard or unknown amino acids 
    will be shown as a number representing its integer mass enclosed in '()'.
    """
    MasstoAAdict = {0: '', 57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: '(I/K)', 114: 'N', 115: 'D', 128: '(K/Q)', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}
    newpeptides = []
    for lst in listofpeptidelists:
        templist = []
        for mass in lst:
            if mass in MasstoAAdict:
                templist.append(MasstoAAdict[mass])
            else:
                templist.append('(' + str(mass) + ')')
        temppep = ''.join(templist)
        newpeptides.append(temppep)
    return newpeptides

## TEST DATA
    
#Spectrum_in = '57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493'    
#M_in = 20
#N_in = 60


Spectrum_in = '328 887 229 540 655 128 584 128 688 360 1143 889 832 532 129 483 1187 872 1013 1129 815 912 370 1116 412 372 559 584 418 172 584 916 1116 1059 1187 660 698 461 115 200 313 603 128 1002 1116 1015 428 469 1044 242 1115 647 1116 332 71 186 57 589 256 1001 1244 1015 556 712 114 532 185 229 315 1173 57 816 1017 303 456 784 826 1072 475 103 788 788 300 431 759 557 1145 460 660 232 101 229 546 227 889 485 227 769 1015 1141 355 641 712 357 813 597 355 874 231 1012 775 685 243 884 931 929 128 357 456 99 1130 660 687 1017 944 941 704 761 783 1058 429 0 887 988'

M_in = 17
N_in = 359
 
Spectrum_in = list(map(int, Spectrum_in.split(" ")))  

petideseq = CyclopeptideSequencing(Spectrum_in, M_in, N_in) 
print(petideseq)
print(ConvertPeptideMassToAA(petideseq))



