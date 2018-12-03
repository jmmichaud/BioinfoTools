#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 13:10:12 2018

@author: jennifermichaud
"""

import copy
import random
import numpy as np


class PairedReadDNASeqData(object):
    """Inputs  paired reads as a 2-row numpy array where row 0 
    holds the first read list and row 1 the second read list and reads in 
    same column are paired and at read pair distance (d).
    """
    def __init__(self, rawreads, d):
        self.rawreads = rawreads
        self.d = d
    
    def GetParams(self):
        """ Determines number of reads, length of reads, and if .  If read length is uniform
        sets reads as kmers and k as the kmer lenght.  All of these parameters
        are then usuable by other functions.
        """
        self.count1 = len(self.rawreads[0,])
        self.count2 = len(self.rawreads[1,])
        self.pairedunequal = 0
        if self.count1 == self.count2:
            self.readcount = self.count1
        elif self.count2 < self.count1:
            self.readcount = self.count2
            np.delete(self.rawreads, [self.count1-1], axis=1) #removes columns that lack paired reads
            self.pairedunequal += self.count1 - self.count2
        elif self.count1 < self.count2:
            self.readcount = self.count1
            np.delete(self.rawreads, [self.count2-1], axis=1) #removes columns that lack paired reads
            self.pairedunequal += self.count2 - self.count1
        readlengthlist1 = [len(read) for read in self.rawreads[0,]]
        readlengthlist2 = [len(read) for read in self.rawreads[1,]]
        self.readlengthdict = {}  #generate dictionary for how often a read lenth is repeated
        for rlen in readlengthlist1: 
            if rlen in self.readlengthdict:
                self.readlengthdict[rlen] += 1
            else:
                self.readlengthdict[rlen] = 1
        for rlen in readlengthlist2: 
            if rlen in self.readlengthdict:
                self.readlengthdict[rlen] += 1
            else:
                self.readlengthdict[rlen] = 1
        for key in self.readlengthdict:
            if len(self.readlengthdict) == 1:
                self.readlength = key
                self.uniform = "Y"
            else:
                self.readlength = "Variable"
                self.uniform = "N"
        if self.uniform == "Y":  #Room to ammend code to allow exceptions for non uniform read lengths
            self.kmers = self.rawreads
            self.k = self.readlength
        if self.pairedunequal > 0:
            print("\n" + self.pairedunequal + " unpaired reads were removed.")
        print("\nNumber of reads: " + str(self.readcount))
        print("Readlength is uniform?: " + self.uniform)
        print("Readlength: " + str(self.readlength))
        
        return self.readcount, self.readlength, self.uniform, self.readlengthdict, self.kmers, self.k
    
    def DeBruijnPaired(self):
        """Inputs paired DNA strings of same kmer length, as a 2-row numpy array 
        where row 0 holds the first kmer list and row 1 the second kmer list and 
        reads in same column are paired.
        Outputs DeBruijn graph as a dictionary of k-1 pairs represented as strings
        "first-k-1-mer second-k-1-mer" 
        """
        try:
            self.kmers
            self.readcount
        except NameError:
            self.GetParams()
        self.dbdict = {} 
        for i in range(self.readcount): #generate left and right edges from each kmer pair
            kmer1 = self.kmers[0,i]
            kmer2 = self.kmers[1,i]
            leftedge = kmer1[:-1] + " " + kmer2[:-1] #concatenate edges to allow use in dictionary
            rightedge = kmer1[1:] + " " + kmer2[1:]
            if leftedge not in self.dbdict:
                self.dbdict.setdefault(leftedge, []).append(rightedge) #adds edge to dictionary with a list as entry, appends right edge
            else:
                self.dbdict[leftedge].append(rightedge)
        return self.dbdict
    
    def InOutDict(self):
        """Inputs a DeBruijn graph as a dictionary in the format:
        {left edge: [right edge1, right edge2...],...}.
        Ouputs the number in degree and out degree of all nodes in a DB graph
        as a dictionary in the format:{node1 : [#indegree, #outdegree]...}. 
        """
        try:
            self.dbdict
        except NameError:
            self.DeBruijnPaired()  
        self.inoutdict = {}
        for key in self.dbdict: 
            if key not in self.inoutdict: #add new dictionary entry
                self.inoutdict[key] = [0, 0]
                for item in self.dbdict[key]:
                    self.inoutdict[key][1] += 1 #add outdegrees
                    if item not in self.inoutdict: #add new dictionary entry
                        self.inoutdict[item] = [0, 0]
                        self.inoutdict[item][0] = 1 #add indegree
                    else:  #for existing entry
                        self.inoutdict[item][0] += 1 #add indegree
            else:  #for existing entry
                for item in self.dbdict[key]:
                    self.inoutdict[key][1] += 1  #add outdegrees
                    if item not in self.inoutdict: #add new dictionary entry
                        self.inoutdict[item] = [0, 0]
                        self.inoutdict[item][0] = 1 #add indegree
                    else: #for existing entry
                        self.inoutdict[item][0] += 1 #add indegree
        return self.inoutdict
    
    
    def StartEndNodes(self):
        """Inputs a dictionary that enumerates the number of indegrees and outdegrees
        for all nodes in a DeBruin Graph in the format: {node1 : [#indegree, #outdegree]...}.
        Outputs the starting node for a semi-balanced path. A semi-balanced path is a path
        where all nodes have the same number of indegrees as outdegrees except for 2
        nodes that signify the begining and end of the path.  Start node has one less 
        indegree and the end node has one less outdegree.
        Outputs two variables start node of a Eularian path as a string and a string
        value that indicates whether the path is balanced. 'B' means fully balanced 
        (circular), 'S' semi-balanced only two nodes that are unbalanced.
        """
        try:
            self.inoutdict
        except NameError:
            self.InOutDict()
        self.startnodes = []
        self.endnodes = []
        self.balanced = 'U' # 'BU' balanced(unbranched), 'BB' balanced(branched), 'SU' semi-balanced (unbranched), ''SB' semi-balanced (branched), 'U' unbalanced
        unbcount = 0
        degreemax = 'N'
        branched = 'N'
        for key in self.inoutdict:
            outminusmin = self.inoutdict[key][1] - self.inoutdict[key][0]
            if outminusmin == 1:
                self.startnodes.append(key)
                unbcount += 1
            elif outminusmin == -1:
                unbcount += 1
                self.endnodes.append(key)
            elif abs(outminusmin) >1:
                unbcount += abs(outminusmin) 
                degreemax = 'Y'
            elif outminusmin == 0:
                if self.inoutdict[key][1] + self.inoutdict[key][0] > 2:
                    branched = 'Y'     
        if unbcount == 0:
            if branched == 'Y':
                self.balanced = 'BB'
            else:
                self.balanced = 'BU'
        elif unbcount <= 2:
            if branched == 'Y':
                self.balanced = 'SB'
            else:
                self.balanced = 'SU'
        elif degreemax == 'Y':
            self.balanced = 'U'
        return self.startnodes, self.endnodes, self.balanced

    
    def EulerianCycle(self):
        """Inputs a DeBruijn graph as a dictionary in the format:
        {left edge: [right edge1, right edge2...],...}.
        Outputs a Eulerian cycle that visits each edge once using DB Graph for
        unbranched balanced and semibalanced DeBruijn Graphs as a list of 
        ordered edges. 
        """
        try:
            self.dbdict
            self.inoutdict
            self.startnodes
            self.balanced
        except NameError:
            self.StartEndNodes()
        dict2 = copy.deepcopy(self.dbdict)
        self.eulpath = []
        temppath = []
        startposition = self.startnodes[0]
        if self.balanced == 'BU':    
            startposition = random.choice(list(dict2))
        elif self.balanced == 'U':
            raise Exception('The DeBruijn Graph indicates an unbalanced path. Program terminated.')
        curredge = startposition   
        while dict2 != {}:  #deletes edges from dictionary as they are used
            if curredge in dict2: #if leftedge in dictionary, looks up entry
                self.eulpath.append(curredge)
                randpos = np.random.randint(0,len(dict2[curredge])) #choose random rightedge if more than one
                position = dict2[curredge][randpos]
                dict2[curredge].remove(dict2[curredge][randpos]) #remove edge from dictionary
                if dict2[curredge] == []: #delete left edge (key) if no right edges remain
                    del dict2[curredge]
                curredge = position #set position to new left edge
            elif self.eulpath != []: #if edges still exist and eulpath not empty, add edge to temp list and remove from eulpath
                temppath.append(curredge) #adds dead end location to temppath
                curredge = self.eulpath.pop(-1) #remove last location from eulerian path and set to current location
        self.eulpath.append(curredge)
        while temppath != []:
            self.eulpath.append(temppath.pop(-1))
        return self.eulpath #outputs the list of  ordered kmer strings and the balanced parameter
    
    def GenomePathReconstruction(self, patterns):
        """Inputs a sequence of k-mers Pattern(1), … ,Pattern(n) such that the last k - 1 
        symbols of Patterni are equal to the first k-1 symbolsb of Patterni+1 for 1 ≤ i ≤ n-1.
        outputs a string Text of length k+n-1 such that the i-th k-mer in Text is 
        equal to Patterni  (for 1 ≤ i ≤ n)"""
        try:
            self.balanced
        except NameError:
            self.StartEndNodes()
        sequence = []
        sequence.append(patterns[0])
        km = len(patterns[0])
        for i in range(1,len(patterns)):
            sequence.append(patterns[i][-1])
        self.assembled = ''.join(sequence)
        if self.balanced == 'BU':
            self.assembled = self.assembled[:-km]
        return self.assembled
    
    def SequencefromPairedReads(self):
        """Inputs paired reads as a 2-row numpy array where the columns correspond
        to paired reads, a kmer length (k), and a distance between paired reads (d).  
        Outputs a reconstructed DNA sequence using caluclation of a DeBruijn graph
        and Eularian path. Returns error message if no overlap exists.
        """
        try:
            self.kmers
            self.k
            self.d
            self.eulpath
            self.balanced
        except NameError:
            self.EulerianCycle()
        def GenomePathReconstruction(patterns, bal):
            sequence = []
            sequence.append(patterns[0])
            km = len(patterns[0])
            for i in range(1,len(patterns)):
                sequence.append(patterns[i][-1])
            self.assembled = ''.join(sequence)
            if self.balanced == 'BU':
                self.assembled = self.assembled[:-km]
            return self.assembled
        pathlist = []
        firstpatterns = []
        secondpatterns = []
        for string in self.eulpath:
            pathlist.append(string.split(" "))
        for pair in pathlist:
            firstpatterns.append(pair[0])
            secondpatterns.append(pair[1])
        prefixstring = GenomePathReconstruction(firstpatterns, self.balanced)
        suffixstring = GenomePathReconstruction(secondpatterns, self.balanced)
        overlap = 1
        if prefixstring[self.k+self.d:]!= suffixstring[:-self.k-self.d]: #check if overlap betweeen paired strings matches
            overlap = 0
        if not overlap:
            return "there is no string spelled by the gapped patterns"
        self.assembly = prefixstring + suffixstring[self.k+self.d:]
        return self.assembly 
        
    def PairedNonBranchingContigs(self):
        """ Inputs a DeBruijn graph as a dictionary in the format:
        {left edge: [right edge1, right edge2...],...}, a dictionary that enumerates 
        the number of indegrees and outdegrees for all nodes in a DeBruin Graph in 
        the format: {node1 : [#indegree, #outdegree]...} and the output of the function
        StartEndNodes() that contains a list of start nodes. Takes unbalanced graphs
        and generates the longest contigs that can be assembled as nonbranching 
        contigs. Generally uses look up of next nodes and if balanced (1 indegree, 
        1 out) adds to the contig.  End or branched nodes (exists as greater than 
        1 in- and out- degrees) terminate contigs and if a branched node new contigs
        are started.  Avoids assembly of contigs that have more than one possible
        order.
        Output is a list of non branching contigs that include potenial overlaps 
        between contigs of k-1 and a numpy array of contigs that serves as a 
        mini-scaffold demonstrating the overlaps between contigs. Each row in the 
        array represents a branched contig, the columns [0,...n] order the contigs 
        from beginning to end each as a list where adjacent columns represent adjacent 
        contigs. Branched contigs where multiple possibilities exist are listed in 
        the same list/column.  For example, the array [[ABCDE], [CDEFGHI, CDEQPRGHI], 
        [GHIJKLM]] represents two possible contigs: 'ABCDEFGHIJKLM' and
        'ABCDEQPRGHIJKLM'. [[list(['CAAATGCAT']), list(['CATACGCTCA', 'CATCA']), 
        list(['TCACCCAG', 'TCAT'])]]. For more complex branching columns may not
        always represent adjacent contigs.   
        """
        try:
            self.dbdict
            self.inoutdict
            self.startnodes
            self.balanced
        except NameError:
            self.StartEndNodes()
        def TrimNodes(remainingdegrees, remainingpaths, node):
            """Removes nodes and node counts from appropriate dictionaries"""
            remainingdegrees[node][1] -= 1
            if remainingdegrees[node][0] == 0 and remainingdegrees[node][1] == 0: #if no remaining in- or out- degrees delete nodes
                del remainingpaths[node] 
                del remainingdegrees[node]
            return remainingdegrees, remainingpaths
        self.contiglist1 = []
        self.contiglist2 = []
        lastedges = []
        remainingpaths = copy.deepcopy(self.dbdict)  #duplicated of dbdict and inout dict are made so nodes can be deleted as used in contigs.
        remainingdegrees = copy.deepcopy(self.inoutdict)
        numofrows = len(self.startnodes)*2
        numofcolumns = 1
        for key in self.inoutdict: #create an array as large as potential contig fragments, a column allotted for every
            if self.inoutdict[key][1] > 1:
                numofcolumns += 1
        self.contigarray = np.empty((numofrows, numofcolumns), dtype = object)
        for x in range(numofrows):
            for y in range(numofcolumns):
                self.contigarray[x,y] = []
        starts = copy.deepcopy(self.startnodes)
        i = 0
        j=0
        node = random.choice(starts)  #initalize a temp path list with a start node
        starts.remove(node) #remove start node from list
        nonbranchingpaths1 = []
        nonbranchingpaths2 = []
        kmerlist = node.split(' ')
        nonbranchingpaths1.append(kmerlist[0])
        nonbranchingpaths2.append(kmerlist[1])
        while remainingpaths != {}:
            if remainingdegrees[node][1] > 1: #if more than one outdegree or terminal end reached
                pass
            elif node in remainingpaths and remainingpaths[node] != []: #if an internal node already used skips it
                nextnode = random.choice(remainingpaths[node]) #add next node by random choice if muliple outgoing nodes exist
                remainingpaths[node].remove(nextnode) #remove selected outgoing node from db graph list
                if remainingdegrees[node][0] !=0: #if not a start node subtracts 1 from indegrees
                    remainingdegrees[node][0] -= 1 
                remainingdegrees, remainingpaths = TrimNodes(remainingdegrees, remainingpaths, node)
                node = str(nextnode)
                kmerlist = node.split(' ')
                nonbranchingpaths1.append(kmerlist[0][-1])
                nonbranchingpaths2.append(kmerlist[1][-1])
                while remainingdegrees[node][1] ==1 and remainingdegrees[node][0]== 1: #while encountering edges with 1 in- and 1 out- degree           
                    nextnode = random.choice(remainingpaths[node])
                    remainingpaths[node].remove(nextnode)
                    if remainingdegrees[node][0] !=0: #if not a start node subtracts 1 from indegrees
                        remainingdegrees[node][0] -= 1 
                    remainingdegrees, remainingpaths = TrimNodes(remainingdegrees, remainingpaths, node)
                    node = str(nextnode)
                    kmerlist = node.split(' ')
                    nonbranchingpaths1.append(kmerlist[0][-1])
                    nonbranchingpaths2.append(kmerlist[1][-1])
            if remainingdegrees[node][0] !=0: #if not a start node subtracts 1 from indegrees
                remainingdegrees[node][0] -= 1 
            if node not in lastedges and node in remainingpaths:
                if remainingdegrees[node][1] !=0: 
                    lastedges.append(node) #keep last edges to start next contig
                elif remainingdegrees[node][1] == 0: #delete end node once reached
                    del remainingdegrees[node]
            if lastedges != []:        
                node = lastedges[0] #usefirst entry in list
            self.contiglist1.append(''.join(nonbranchingpaths1))
            self.contiglist2.append(''.join(nonbranchingpaths2))
            self.contigarray[i,j].append(''.join(nonbranchingpaths1))
            self.contigarray[i+1,j].append(''.join(nonbranchingpaths2))
            try:  #initializes the previous last choice variable without overwriting if it already exists
                prevlast
            except NameError:
                prevlast = ""
            if remainingdegrees[node][1] > 0: #if there are outdegree for the node make a new contig with it as the start
                if prevlast != node:
                    j += 1
                nonbranchingpaths1 = []
                nonbranchingpaths2 = []
                kmerlist = node.split(' ')
                nonbranchingpaths1.append(kmerlist[0])
                nonbranchingpaths2.append(kmerlist[1])
                nextnode = random.choice(remainingpaths[node])
                remainingpaths[node].remove(nextnode)
                remainingdegrees, remainingpaths = TrimNodes(remainingdegrees, remainingpaths, node)
                if remainingdegrees[node][1] == 0:
                    lastedges.remove(node)   
                prevlast = str(node) #needed only for making numpy map
                node = str(nextnode)
                kmerlist = node.split(' ')
                nonbranchingpaths1.append(kmerlist[0][-1])
                nonbranchingpaths2.append(kmerlist[1][-1])
            elif remainingdegrees[node][1] == 0: #if endnode is reached
                if remainingdegrees[node][0] == 1: #remove endnode 
                    del remainingdegrees[node]
                if starts != []: #start new contig with a new start node
                    i += 2
                    j = 0
                    node = random.choice(starts)
                    starts.remove(node)
                    nonbranchingpaths1 = []
                    nonbranchingpaths2 = []
                    kmerlist = node.split(' ')
                    nonbranchingpaths1.append(kmerlist[0])
                    nonbranchingpaths2.append(kmerlist[1])
                    nextnode = random.choice(remainingpaths[node])
                    remainingpaths[node].remove(nextnode)
                    remainingdegrees, remainingpaths = TrimNodes(remainingdegrees, remainingpaths, node)
                    node = str(nextnode)
                    kmerlist = node.split(' ')
                    nonbranchingpaths1.append(kmerlist[0][-1])
                    nonbranchingpaths2.append(kmerlist[1][-1])
            contig1 = ''.join(nonbranchingpaths1)    
            contig2 = ''.join(nonbranchingpaths2) 
            if remainingpaths == {} and contig1 not in self.contiglist1:  #adds last fragment to list
                self.contiglist1.append(contig1)
                self.contigarray[i,j].append(contig1)
            if remainingpaths == {} and contig2 not in self.contiglist2:  #adds last fragment to list
                self.contiglist1.append(contig1)
                self.contigarray[i+1,j].append(contig2)
            print(remainingpaths)
            print(remainingdegrees)
            print(self.contigarray)
            print("node "+ str(node))
            print("nonbranchingpaths "+ str(nonbranchingpaths1), str(nonbranchingpaths2))
            print("last edges "+ str(lastedges))
        return self.contiglist1, self.contiglist2, self.contigarray  


    def AssemblyPairedNonBranchingContigs(self):
        """Inputs a contig array from non branching contigs, a kmerlength (k), and a paired read 
        distance (d).  Condenses paired reads by the k+d overlap between them. 
        Outputs a numpy array that contains non-branching contigs where each row in the 
        array represents a branched contig, the columns [0,...n] order the contigs 
        from beginning to end each as a list where adjacent columns represent adjacent 
        contigs. Branched contigs where multiple possibilities exist are listed in 
        the same list/column. Where a paired read could not be resolved (the overlap)
        between them does not match) "n/a" is stored.
        """
        try:
            self.contigarray
            self.k
            self.d
        except NameError:
            self.PairedNonBranchingContigs()
        r, c = self.contigarray.shape
        self.assembledcontigs = np.empty((int(r/2), c), dtype = object) #create empty array with a list for each entry and half the rows of contigarray
        for x in range(int(r/2)):
            for y in range(c):
                self.assembledcontigs[x,y] = []
        i = 0
        l = 0
        while i < r:   #combine contigs from corresponding rows using an overlap
            for j in range(c):
                prefixlist = self.contigarray[i,j]
                suffixlist = self.contigarray[i+1,j]
                overlap = 1
                for m in range(len(prefixlist)):
                    prefix = prefixlist[m]
                    suffix = suffixlist[m]
                    if prefix[self.k+self.d:]!= suffix[:-self.k-self.d]: #check if overlap betweeen paired strings matches
                        overlap = 0
                    if overlap:
                        self.assembledcontigs[l,j].append(prefix + suffix[self.k+self.d:])
                    else:
                        self.assembledcontigs[l,j].append("n/a")     
            i += 2
            l += 1
        return self.assembledcontigs

def ReadPairsConvert(pairedreads):
    """Inputs a list of paired reads in the format "kmer1|kmer2..." 
    and outputs paired reads in format of np array. Row 1 is kmer1
    Row2 is kmer2. Column indicates pairs"""
    pread_array = np.empty((2, len(pairedreads)),dtype=object)
    k = int((len(pairedreads[0])-3)/2)
    for i in range(len(pairedreads)): 
        pr1 = pairedreads[i][0:k+1]
        pr2 = pairedreads[i][k+2:]
        pread_array[0,i] = pr1
        pread_array[1,i] = pr2
    return pread_array   



##### TEST DATA ###################################

pr_in = """GAGA|TTGA
TCGT|GATG
CGTG|ATGT
TGGT|TGAG
GTGA|TGTT
GTGG|GTGA
TGAG|GTTG
GGTC|GAGA
GTCG|AGAT"""
pr_in = pr_in.split("\n")

k_in = 4
d_in = 2
seq_in = ReadPairsConvert(pr_in)
## correct assembly 'GTGAGATGTTGA'
    

##################################################

# Code to perform paired read assembly   

# Initalize reads and check parameters
seqdata =  PairedReadDNASeqData(seq_in, d_in) 
seqdata.GetParams()

# Generate DeBruijn Graph, In and Outdegree parameters including if the Debruijn
# is balanced,
seqdata.DeBruijnPaired()
seqdata.InOutDict()
seqdata.StartEndNodes()

#Determine type of function to use based upon how balanced the DeBruijn Graph is

#If data is unbranched run standard Eularian Function
if seqdata.balanced == 'BU' or seqdata.balanced == 'SU':
    seqdata.EulerianCycle()
    seqdata.SequencefromPairedReads()
    
# If data is branched run function to assemble into all possible unbranched contigs   
# Alteration of non-braching contigs function to allow for a non empty remaining paths could allow its addtion here    
elif seqdata.balanced == 'BB' or seqdata.balanced == 'SB' or seqdata.balanced == 'U':
    seqdata.PairedNonBranchingContigs()
    seqdata.AssemblyPairedNonBranchingContigs()
    

print("\n\nAssembly Data")
print("===========================")
if seqdata.balanced == 'BU':
    print('The read data suggests a circular DNA strand.')
    print('\nAn intact sequence was able to be assembled from the provided reads.')
    print('\nAssembled sequence: ' + str(seqdata.assembly))
elif seqdata.balanced == 'BB':
    print('The read data suggests a circular DNA strand.')
    print('Non-branching contigs were assembled from the provided reads. \nThe contigs in the first and last contigs likely have a redundant region of overlap.')
    print('\nAssembled contigs(numpy array): ' + str(seqdata.assembledcontigs))
    print('\nEach row indicates an independant contig. \nAdjacent columns indicate contigs with overlapping edges.')
elif seqdata.balanced == 'SU':
    print('The read data suggests a linear DNA strand.')
    print('\nAn intact sequence was able to be assembled from the provided reads.')
    print('\nAssembled sequence: ' + str(seqdata.assembly))
elif seqdata.balanced == 'SB':
    print('The read data suggests a linear DNA strand.')
    print('Non-branching contigs were assembled from the provided reads.')
    print('\nAssembled contigs(numpy array): ' + str(seqdata.assembledcontig))
    print('\nEach row indicates an independant contig. \nAdjacent columns indicate contigs with overlapping edges.') 
elif seqdata.balanced == 'U':
    print('\n\nAn unbalanced or semi-balanced DeBruijn Graph was not found from the read data.')
    print('Contigs were assembled from available data.')
    print('\nAssembled contigs(numpy array): ' + str(seqdata.assembledcontig))
    print('\nEach row indicates an independant contig. \nAdjacent columns indicate contigs with overlapping edges.') 

    