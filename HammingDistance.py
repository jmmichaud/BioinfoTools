#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 13:07:44 2018

@author: jennifermichaud
"""

def HammingDistance(text1,text2):
    """Inputs two equal length DNA strings. Outputs the hamming distance between
    them as an integer.  Each mismatch is a score of 1."""
    hamming = 0
    for i in range(len(text1)):
        if text1[i] != text2[i]:
            hamming += 1
    return hamming


#TEST DATA   
#dna1 = 'AAATAGCTGCGTCTCAAGTTCGCGATTCGGAAACCCGGACGTACTCATAGGCCGCTAACCAGATATACGCCTAAGATTTAAGATATGCGAGGCCCCGCAATAAATGTCAGCCTGTCGCAATGAAGCCCGACCCTACGTATTGACTTGTTGGTATCTGTATTATAATATAGCCGATCGGCGCCCGTTCACGTAGTAGGGAACTTCTCACGGGGATAGGACCAGCCTTGCGTTCATGTTCACAGCCCCAACAGTTCCCCACAGCGGGGGTTTGGGCCTAACCCGGATTTCCGCCCTTCCAACCAAGACTCCGGACACTGTGTCGGCGGTCGTGCTATAGGTCAACTCCCCATAAGTAGGTCACTAAAGCACGCAGACTATTTCACGTTTGTCCGACGCGCCCCTACTATCCTGGACACGTACGGTTTCTCTGATTAAGCATATCCTCCACTAATGTCACTGAGTTTTGTTCACGGTGATAACGCACGAGTTGAGCCCGCGCCTAGGCATACCAGCAGATTCCCATTGAAAAGGGGTGTGTACAATGAAGTACGCGCTCGTGAGGGATGACTTTCGGGCGTCCCGCATGCCACGTCTTTATTTCGAGTCTCCGGCTTAAAATCCCACAATAACGCGACCCTGTGATCCATGTTTGGTCCCATAAAGGGAGCTCTCACTAAAGCAACATAGAGATGTCGGCCGACTCTTTATAGAAACGGTGCCCATTCAAGCGGGCTATCCTTATACAGTGAAGCGCTGTGCCTGTCAGGTGGAAGTTGTTCTACGAAATAGCCAGGTAGCTAACTGGAACCAAGAACTCACGGACAATCCTGTTCCGCATTCGAGTAGTGCGACCCGCTCCAGTAGTAAGAGATCTAATCCTGATGACAGTTAGTTATCTGGATTCGTCCCCCTGGCAGATAAGCGAGTGCGAGAAGTCTCGCTCCCAGTCAACGGGCCGCCTCGGCCGTATTATTAAGTGGAGTACGATCACGTAGGATGTACACGGCGTACGACCGGAGCTCGGGTTTTAGTGGAAAATCAACCAGACCCGTGCAGTAATACCCAGGTTGACTGTTCCTTGATGGCGGTGAATCGGGCTCGCTAGTCCCGGATACGCGTCATCCTCCTCGAAAGCTTAAGCCGCAGGCACACTCGTACGACGTGACTACGTCGCTGGAACGAGGTTCAG' 
#dna2 = 'AGTAAAATTGCTGCGCCTACCGATCACAGCGACGTACCACGCGCGAAGGATTTTGCTTCTATCCACTCCACTTTAACCATCTGTGTGGTTCGTAACGAAAATCCACGAAGAAGGCTGCTCAGTTAAACGGCCTGCATGTCGGCAGGACCGTTTGGTTTCGCAAGGTTATCACGGTCCTCTCCGAACAATGTAGTCTTTCAATAGTCAGCTGATGTAGGTCCTCTCGAGGTGTTCGGACGGCGAAACTCCAACAGTCAACGTGCTAGGACGTATGCCCCTGGCCCTCACGTTGAGCCTGACGCAGCTGGATTGGAGGACGGGCCGAAGTCCAGGAAGTACAAGGGTGCCGACTCAATATAACCAGAAAGAAATCTGAGCAAAACTCGCGTAGTAGTATTCGGGGGTCTTACGAGGCCTCTGGTCAACCGCACCGTTGGAGTTGCATAGTGCCAATAGGTAGCGAGGAGGGCATGGTCGTCAAACCGGCCAACGGATATAACAATCAGAGATTCTAACTAGTTCCCCGTAGTAGGCATGACGATTGTCACTTCGCCTAAGACACACTATTAAATTGGTAGACACAAGCCAGGTCCAGAGGATCACCTATGAGCAGTGATCCTGGTCCTCACTACGGCATAATTTGTCCAACAAAAAGGAATTTGAGCCAATGTTACCGCCGCATGAAGAGATGTTATATCTGATGTAGACTGTCCGTCCACACTAGATAGAATGAACGTATAACACCAACGTTGATGAACTCGTAGACACACATATCGTATGAGGCGAACGACCACGCCGGATTGATTGGCACACGAAAGTGCACGCTTGTCCTCCCCCGTGCAATAATTGACAATCGTTAAGTTGCGCCGTATACCTTATAAAATTGGTCCTACGAGATCGGTTTTGAGGCATGGGCATGTATACCAATTCATAAAGTACCCTGTCGTACCTGGTTACTTCGTGGCGCCAACTCAATAGGAGGGAACACCGAGGGGCGGACGGCAGCTGAGTTTCGTGAAAGTTATACGTTGAAATAGGTACTTCCGATAGGAGTGGGTGGCTAAGACCTCGAGGCATCGGATGTTTGTCTTCAGCGCGGACTCGTATCCACCAGGAGAGGTATGTCTGATGTTTTGTAGTGTACTTGTGCAAACAGTTCTCGGGATGGGACCGTTTGTAGTCGGCGTAA'
##corr_answer = 893
#
#print(HammingDistance(dna1,dna2))