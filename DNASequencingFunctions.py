#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 11:49:59 2018

@author: jennifermichaud
"""

import numpy as np
import random
import copy


def ReverseComplement(text):
    """ Inputs a DNA string.  Outputs the reverse complement of the DNA strand."""
    text = text.upper()
    compd = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    revtext = text[::-1]
    revcomp = []
    for let in revtext:
        revcomp.append(compd[let])
    return ''.join(revcomp)


def KmerComposition(dnastring, k):
    """Inputs a DNA string (dnastring) and a length(k).
    Outputs a list of DNA strings that deconvolute the string into a list of 
    k-lenth strings where each string has a k-1 base overlap with at least one 
    other DNA string and they are ordered in the list according to their overlap 
    such that the DNA string [i][1:] matches DNA string [i+1][:k].
    """
    kmers = []
    for n in range(len(dnastring)-k +1):
        kmer = dnastring[n:n+k]
        kmers.append(kmer)
    return kmers


def DeBruijnG(patterns):
    """Inputs a collection of k-length kmers (patterns) as list of strings.
    Outputs the overlap graph in the form of an adjacency list: 
    left edge -> right edge1, right edge2 . Adjacency is determined by the 
    overlap of a suffix(length k-1) of one pattern with a prefix(k-1) of another 
    pattern. Each pattern edge k-1 gets a node to assign overlaps for a resulting 
    k-2 overlaps between k-1 nodes. DeBruijn Graphs generate graphs of edges and 
    allow for the use of Eulerian paths (visit each edge once) which are easier 
    to solve than Hamiltonian paths (visit each vertex once).
    """
    dbmap = {} #map of nodes to neighbors
    nodes = {} #maps k-1 mers to nodes
    for kmer in patterns:  #generates nodes
        kmL, kmR = kmer[:-1], kmer[1:]
        nodeL, nodeR = None, None
        if kmL in nodes:
            nodeL = nodes[kmL]
        else:
            nodeL = nodes[kmL] = kmL
        if kmR in nodes:
            nodeR = nodes[kmR]
        else:
            nodeR = nodes[kmR] = kmR
        dbmap.setdefault(nodeL, []).append(nodeR)
    outputlist = []
    for edge in dbmap:
        outputlist.append(edge + " -> " + ",".join(dbmap[edge]))
    return "\n".join(outputlist)


def DeBruijn(patterns):
    """Inputs a collection of k-length kmers (patterns) as list of strings.
    Outputs the overlap graph in the form of an dictionary: 
    {left edge: [right edge1, right edge2...],...}. Adjacency is determined by the 
    overlap of a suffix(length k-1) of one pattern with a prefix(k-1) of another 
    pattern. Each pattern edge k-1 gets a node to assign overlaps for a resulting 
    k-2 overlaps between k-1 nodes. DeBruijn Graphs generate graphs of edges and 
    allow for the use of Eulerian paths (visit each edge once) which are easier 
    to solve than Hamiltonian paths (visit each vertex once).
    """
    dbmap = {} #map of nodes to neighbors
    for kmer in patterns:  #generates nodes
        leftkmer, rightkmer = kmer[:-1], kmer[1:]
        if leftkmer not in dbmap:
            dbmap.setdefault(leftkmer, []).append(rightkmer) #adds edge to dictionary with a list as entry, appends right edge
        else:
            dbmap[leftkmer].append(rightkmer) #appends rightedge to node
    return dbmap

            
##TEST DATA
#inputdna = """GAGG
#CAGG
#GGGG
#GGGA
#CAGG
#AGGG
#GGAG"""
#dna_in= inputdna.split("\n")
#
#
#print(DeBruijnG(dna_in))
#print(DeBruijn(dna_in))


def DeconvoluteDBgraph(dbgraph):
    """Inputs a DeBruijn Graph in the format : left edge -> right edge1, right edge2 
    separated by '\n'.
    Ouputs the DB graph in dictionary format:
    {left edge: [right edge1, right edge2...],...}.
    """
    graphlist = dbgraph.split("\n")
    dbdict = {}
    for line in graphlist:
        allspace =  line
        allspace = allspace.replace(',', ' ')
        edgelist = allspace.split(' ')
        dbdictlist = []
        for i in range(2,len(edgelist)):
            dbdictlist.append(edgelist[i])
        dbdict[edgelist[0]] = dbdictlist
    return dbdict


##TEST DATA
#graph_in ="""1 -> 2
#2 -> 3
#3 -> 4,5
#6 -> 7
#7 -> 6"""  
#    
#print(DeconvoluteDBgraph(graph_in))

    
def InOutDict(dbdict):
    """Inputs a DeBruijn graph as a dictionary in the format:
    {left edge: [right edge1, right edge2...],...}.
    Ouputs the number in degree and out degree of all nodes in a DB graph
    as a dictionary in the format:{node1 : [#indegree, #outdegree]...}. 
    """
    inoutdict = {}
    for key in dbdict: 
        if key not in inoutdict: #add new dictionary entry
            inoutdict[key] = [0, 0]
            for item in dbdict[key]:
                inoutdict[key][1] += 1 #add outdegrees
                if item not in inoutdict: #add new dictionary entry
                    inoutdict[item] = [0, 0]
                    inoutdict[item][0] = 1 #add indegree
                else:  #for existing entry
                   inoutdict[item][0] += 1 #add indegree
        else:  #for existing entry
           for item in dbdict[key]:
               inoutdict[key][1] += 1  #add outdegrees
               if item not in inoutdict: #add new dictionary entry
                   inoutdict[item] = [0, 0]
                   inoutdict[item][0] = 1 #add indegree
               else: #for existing entry
                   inoutdict[item][0] += 1 #add indegree
    return inoutdict


def StartNode(inoutdict):
    """Inputs a dictionary that enumerates the number of indegrees and outdegrees
    for all nodes in a DeBruin Graph in the format: {node1 : [#indegree, #outdegree]...}.
    Outputs the starting node for a semi-balanced path. A semi-balanced path is a path
    where all nodes have the same number of indegrees as outdegrees except for 2
    nodes that signify the begining and end of the path.  Start node has one less 
    indegree and the end node has one less outdegree.
    Outputs two variables start node of a Eularian path as a string and a string
    value that indicates whether the path is balanced. 'B' means fully balanced 
    (circular), 'S' semi-balanced only two nodes that are unbalanced
    """
    startnode = ""
    balanced = 'U' # 'B' balanced, 'S' semi-balanced, 'U' unbalanced
    unbcount = 0
    degreemax = 'N'
    for key in inoutdict:
        outminusmin = inoutdict[key][1] - inoutdict[key][0]
        if outminusmin == 1:
            startnode = key
            unbcount += 1
        elif outminusmin == -1:
            unbcount += 1
        elif abs(outminusmin) >1:
            unbcount += abs(outminusmin) 
            degreemax = 'Y'
    if unbcount == 0:
        balanced = 'B'
    elif unbcount <= 2:
        balanced = 'S'
    elif degreemax == 'Y':
        balanced = 'U'
    return startnode, balanced


def StartEndNodes(inoutdict):
    """Inputs a dictionary that enumerates the number of indegrees and outdegrees
    for all nodes in a DeBruin Graph in the format: {node1 : [#indegree, #outdegree]...}.
    Outputs the starting node for a semi-balanced path. A semi-balanced path is a path
    where all nodes have the same number of indegrees as outdegrees except for 2
    nodes that signify the begining and end of the path.  Start node has one less 
    indegree and the end node has one less outdegree.
    Outputs two variables start node of a Eularian path as a string and a string
    value that indicates whether the path is balanced. 'B' means fully balanced 
    (circular), 'S' semi-balanced only two nodes that are unbalanced
    """
    startnodes = []
    endnodes = []
    balanced = 'U' # 'B' balanceD, 'S' semi-balanced, 'U' unbalanced
    unbcount = 0
    degreemax = 'N'
    for key in inoutdict:
        outminusmin = inoutdict[key][1] - inoutdict[key][0]
        if outminusmin == 1:
            startnodes.append(key)
            unbcount += 1
        elif outminusmin == -1:
            unbcount += 1
            endnodes.append(key)
        elif abs(outminusmin) >1:
            unbcount += abs(outminusmin) 
            degreemax = 'Y'
    if unbcount == 0:
        balanced = 'B'
    elif unbcount <= 2:
        balanced = 'S'
    elif degreemax == 'Y':
        balanced = 'U'
    return startnodes, endnodes, balanced


def EulerianCycle(dbdict):
    """Inputs a DeBruijn graph as a dictionary in the format:
    {left edge: [right edge1, right edge2...],...}.
    Outputs a Eulerian cycle that visits each edge once using DB Graph for
    balanced and semibalanced DeBruijn Graphs. Can alter code to either output
    as a list of edges in the order of the path and the balanced paramter or 
    as string with the edges separated by '->'.
    """
    inoutdict = InOutDict(dbdict)
    dict2 = copy.deepcopy(dbdict)
    eulpath = []
    temppath = []
    startposition, balanced = StartNode(inoutdict)
    if balanced == 'B':    
        startposition = random.choice(list(dict2))
    elif balanced == 'U':
        raise Exception('The DeBruijn Graph indicates an unbalanced path. Program terminated.')
    curredge = startposition   
    while dict2 != {}:  #deletes edges from dictionary as they are used
        if curredge in dict2: #if leftedge in dictionary, looks up entry
            eulpath.append(curredge)
            randpos = np.random.randint(0,len(dict2[curredge])) #choose random rightedge if more than one
            position = dict2[curredge][randpos]
            dict2[curredge].remove(dict2[curredge][randpos]) #remove edge from dictionary
            if dict2[curredge] == []: #delete left edge (key) if no right edges remain
                del dict2[curredge]
            curredge = position #set position to new left edge
        elif eulpath != []: #if edges still exist and eulpath not empty, add edge to temp list and remove from eulpath
            temppath.append(curredge) #adds dead end location to temppath
            curredge = eulpath.pop(-1) #remove last location from eulerian path and set to current location
    eulpath.append(curredge)
    while temppath != []:
        eulpath.append(temppath.pop(-1))
    return eulpath, balanced #outputs the list of  ordered kmer strings and the balanced parameter


##TEST DATA
#kmers_in = """AAAT
#AATG
#ACCC
#ACGC
#ATAC
#ATCA
#ATGC
#CAAA
#CACC
#CATA
#CATC
#CCAG
#CCCA
#CGCT
#CTCA
#GCAT
#GCTC
#TACG
#TCAC
#TCAT
#TGCA"""
#patterns_in  =  kmers_in.split("\n") 
#
#dbgraph_in = DeBruijn(patterns_in)
#print(dbgraph_in)
#print(EulerianCycle(dbgraph_in))
##Output CAA->AAA->AAT->ATG->TGC->GCA->CAT->ATC->TCA->CAT->ATA->TAC->ACG->CGC->
## GCT->CTC->TCA->CAC->ACC->CCC->CCA->CAG


def AssemblyNonBranchingContigs(dbdict, inoutdict, startendnodes):
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
    startnodes, endnodes, balanced = startendnodes
    contiglist = []
    lastedges = []
    remainingpaths = copy.deepcopy(dbdict)  #duplicated of dbdict and inout dict are made so nodes can be deleted as used in contigs.
    remainingdegrees = copy.deepcopy(inoutdict)
    numofrows = len(startnodes)
    numofcolumns = 1
    for key in inoutdict: #create an array as large as potential contig fragments, a column allotted for every
        if inoutdict[key][1] > 1:
            numofcolumns += 1
    contigarray = np.empty((numofrows, numofcolumns), dtype = object)
    for x in range(numofrows):
        for y in range(numofcolumns):
            contigarray[x,y] = []
    starts = copy.deepcopy(startnodes)
    i = 0
    j=0
    node = random.choice(starts)  #initalize a temp path list with a start node
    starts.remove(node) #remove start node from list
    nonbranchingpaths = []
    nonbranchingpaths.append(node)
    while remainingpaths != {}:
        if node in remainingpaths: #if an internal node already used skips it
            nextnode = random.choice(remainingpaths[node]) #add next node by random choice if muliple outgoing nodes exist
            remainingpaths[node].remove(nextnode) #remove selected outgoing node from db graph list
            remainingdegrees[node][1] -= 1
            if remainingdegrees[node][0] !=0: #if not a start node subtracts 1 from indegrees
                remainingdegrees[node][0] -= 1 
            if remainingdegrees[node][0] == 0 and remainingdegrees[node][1] == 0: #if no remaining in- or out- degrees delete nodes
                del remainingpaths[node] 
                del remainingdegrees[node]
            node = str(nextnode)
            nonbranchingpaths.append(node[-1])                              
            while remainingdegrees[node][1] ==1 and remainingdegrees[node][0]== 1: #while encountering edges with 1 in- and 1 out- degree           
                nextnode = random.choice(remainingpaths[node])
                remainingpaths[node].remove(nextnode)
                remainingdegrees[node][1] -= 1
                if remainingdegrees[node][0] !=0: #if not a start node subtracts 1 from indegrees
                    remainingdegrees[node][0] -= 1 
                if remainingdegrees[node][0] == 0 and remainingdegrees[node][1] == 0: #if no remaining in- or out- degrees  delete nodes
                    del remainingpaths[node] 
                    del remainingdegrees[node]           
                node = str(nextnode)
                nonbranchingpaths.append(node[-1])
        if node not in lastedges and node in remainingpaths:
            if remainingdegrees[node][1] !=0: 
                lastedges.append(node) #keep last edges to start next contig
            elif remainingdegrees[node][1] == 0: #delete end node once reached
                del remainingdegrees[node]
        if lastedges != []:        
            node = lastedges[0] #usefirst entry in list
        contiglist.append(''.join(nonbranchingpaths))
        contigarray[i,j].append(''.join(nonbranchingpaths))
        try:  #initializes the previous last choice variable without overwriting if it already exists
            prevlast
        except NameError:
            prevlast = ""
        if remainingdegrees[node][1] > 0: #if there are outdegree for the node make a new contig with it as the start
            if prevlast != node:
                j += 1
            nonbranchingpaths = []
            nonbranchingpaths.append(node)
            nextnode = random.choice(remainingpaths[node])
            remainingpaths[node].remove(nextnode)
            remainingdegrees[node][1] -= 1
            if remainingdegrees[node][0] !=0: #if not a start node subtracts 1 from indegrees
                remainingdegrees[node][0] -= 1 
            if remainingdegrees[node][0] == 0 and remainingdegrees[node][1] == 0: #if no remaining in- or out- degrees degrees delete nodes
                del remainingpaths[node] 
                del remainingdegrees[node]
                lastedges.remove(node)
            prevlast = str(node) #needed only for making numpy map
            node = str(nextnode)
            nonbranchingpaths.append(node[-1])
        elif remainingdegrees[node][1] == 0: #if endnode is reached
            if remainingdegrees[node][0] == 1: #remove endnode 
                del remainingdegrees[node]
            if starts != []: #start new contig with a new start node
                i += 1
                j = 0
                node = random.choice(starts)
                starts.remove(node)
                nonbranchingpaths = []
                nonbranchingpaths.append(node)
                nextnode = random.choice(remainingpaths[node])
                remainingpaths[node].remove(nextnode)
                remainingdegrees[node][1] -= 1
                if remainingdegrees[node][0] !=0: #if not a start node subtracts 1 from indegrees
                    remainingdegrees[node][0] -= 1 
                if remainingdegrees[node][0] == 0 and remainingdegrees[node][1] == 0: #if no remaining in- or out- degrees
                    del remainingpaths[node] #remove
                    del remainingdegrees[node]
                node = str(nextnode)
                nonbranchingpaths.append(node[-1])
        contig = ''.join(nonbranchingpaths)        
        if remainingpaths == {} and contig not in contiglist:  #adds last fragment to list
            contiglist.append(contig)
            contigarray[i,j].append(contig)
    return contiglist, contigarray
                

##TEST DATA
#kmers_in = """AAAT
#AATG
#ACCC
#ACGC
#ATAC
#ATCA
#ATGC
#CAAA
#CACC
#CATA
#CATC
#CCAG
#CCCA
#CGCT
#CTCA
#GCAT
#GCTC
#TACG
#TCAC
#TCAT
#TGCA"""
#patterns_in  =  kmers_in.split("\n") 
#
#
#dbdict_in = DeBruijn(patterns_in)
##print(dbdict_in)
#inoutdict_in = InOutDict(dbdict_in)
##print(inoutdict_in)
#
#startendnodes_in = StartEndNodes(inoutdict_in)
##print(NonBranchingContigs1(dbdict_in, inoutdict_in, 5))
#print(AssemblyNonBranchingContigs(dbdict_in, inoutdict_in, startendnodes_in))
### correct_output = (['CAAATGCAT', 'CATACGCTCA', 'CATCA', 'TCAT', 'TCACCCAG'], 
###array([[list(['CAAATGCAT']), list(['CATACGCTCA', 'CATCA']),list(['TCAT', 'TCACCCAG'])]], dtype=object))


def GenomePathReconstruction(genomepathlist, balanced):
    """inputs a sequence of k-mers Pattern(1), … ,Pattern(n) such that the last k - 1 
symbols of Patterni are equal to the first k-1 symbolsb of Patterni+1 for 1 ≤ i ≤ n-1.
outputs a string Text of length k+n-1 such that the i-th k-mer in Text is 
equal to Patterni  (for 1 ≤ i ≤ n)"""
    sequence = []
    sequence.append(genomepathlist[0])
    k = len(genomepathlist[0])
    for i in range(1,len(genomepathlist)):
        sequence.append(genomepathlist[i][-1])
    assembled = ''.join(sequence)
    if balanced == 'B':
        assembled = assembled[:-k]
    return assembled


##TEST DATA
##balanced (circular) example
#kmers_in = """ATTAC
#TACAG
#GATTA
#ACAGA
#CAGAT
#TTACA
#AGATT"""
#correct_out 'GATTACA'
    
#print(dbgraph_in)
#eulerp, balanced = EulerianCycle(dbgraph_in)
#print(eulerp)
#print(balanced)
#print(GenomePathReconstruction(eulerp, balanced))
#correct_out 'GATTACA'


def MakeReadPair(dnastring, k, d):
    """Inputs a DNA string (dnastring), a kmer length (k), and a read-pair 
    distance (d). Generates a set of k-length kmer read pairs with d distance 
    between them. Ouputs a list in lexicographic order with each pair represented 
    as "kmer1|kmer2". 
    """
    pairlist = []
    for i in range(len(dnastring)-2*k-d+1):
        kmer1 = dnastring[i:i+k]
        kmer2 = dnastring[i+k+d:i+2*k+d]
        pairlist.append(kmer1 + "|" + kmer2)
    return pairlist
    
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
            

def StringSpelledByGappedPatterns(pairedreads, k, d, balanced = 'S'):
    """Inputs  paired reads as a 2-row numpy array where row 0 
    holds the first kmer list and row 1 the second kmer list and reads in 
    same column are paired, a kmer length (k), and read pair distance (d). 
    Outputs the string made by overlap of prefix strings (from kmer1s)
    concatenated with suffixstrings (from kmer2s).  Returns error message
    if no overlap exists"""
    firstpatterns = list(pairedreads[0,])
    secondpatterns = list(pairedreads[1,])
    prefixstring = GenomePathReconstruction(firstpatterns, balanced)
    print(prefixstring)
    suffixstring = GenomePathReconstruction(secondpatterns, balanced)
    print(suffixstring)
    for i in range(k+d+1,len(prefixstring)):
        if prefixstring[i] != suffixstring[i-k-d]:
            return "there is no string spelled by the gapped patterns"
    return prefixstring + suffixstring[-k-d:]

##TEST DATA
#print(MakeReadPair("TAATGCCATGGGATGTT", 3, 2))   
#pairedreads_in = ReadPairsConvert(MakeReadPair("TAATGCCATGGGATGTT", 3, 2))
#print(ReadPairsConvert(MakeReadPair("TAATGCCATGGGATGTT", 3, 2)))
#print(StringSpelledByGappedPatterns(pairedreads_in, 3, 2))



def DeBruijnPaired(pairedreads):
    """Inputs paired DNA strings of same kmer length, as a 2-row numpy array 
    where row 0 holds the first kmer list and row 1 the second kmer list and 
    reads in same column are paired.
    Outputs DeBruijn graph as a dictionary of k-1 pairs represented as strings
    "first-k-1-mer second-k-1-mer" 
    """
    dbmap = {} 
    for i in range(len(pairedreads[0,])): #generate left and right edges from each kmer pair
        kmer1 = pairedreads[0,i]
        kmer2 = pairedreads[1,i]
        leftedge = kmer1[:-1] + " " + kmer2[:-1] #concatenate edges to allow use in dictionary
        rightedge = kmer1[1:] + " " + kmer2[1:]
        if leftedge not in dbmap:
            dbmap.setdefault(leftedge, []).append(rightedge) #adds edge to dictionary with a list as entry, appends right edge
        else:
            dbmap[leftedge].append(leftedge)
    return dbmap
                   

def EulerianPaired(pairedreads):
    """Inputs paired reads as a 2-row numpy array where the columns correspond
    to paired reads. Generates a cycle or path that visits each paired edge once 
    using paired DB Graphs using 2 db dictionaries of each paired unit.
    Outputs paired eularian path as a 2-row numpy array where columns
    correspond to paired path nodes.  The ony difference between this and the 
    non-paired version is the DeBruijnPaired() function.
    """
    dbdict = DeBruijnPaired(pairedreads)
    inoutdict = InOutDict(dbdict)
    dict2 = copy.deepcopy(dbdict)
    startposition, balanced = StartNode(inoutdict)
    curredge = startposition
    eulpath = []
    temppath = []
    if balanced == 'B':    
        startposition = random.choice(list(dict2))
    elif balanced == 'U':
        raise Exception('The DeBruijn Graph indicates an unbalanced path. Program terminated.')
    curredge = startposition
    while dict2 != {}:
        if curredge in dict2: #if leftedge in dictionary looks up entry
            eulpath.append(curredge)
            randpos = np.random.randint(0,len(dict2[curredge])) #choose random rightedge if more than one
            position = dict2[curredge][randpos]
            dict2[curredge].remove(dict2[curredge][randpos]) #remove edge from dictionary
            if dict2[curredge] == []: #delete left edge (key) if no right edges remain
                del dict2[curredge]
            curredge = position #set position to left edge
        elif eulpath != []: #if edges still exist and eulpath not empty, add edge to temp list and remove from eulpath
            temppath.append(curredge) #adds dead end location to temppath
            curredge = eulpath.pop(-1) 
        else: #if edges still exist but dead end reached 
            temppath.append(curredge) #adds dead end location to temppath
            curredge = eulpath.pop(-1) #remove last location from eulerian path and set to current location 
    eulpath.append(curredge)
    while temppath != []:
        eulpath.append(temppath.pop(-1))
    return eulpath, balanced 


def SequencefromPairedReads(pairedreads, k, d):
    """Inputs paired reads as a 2-row numpy array where the columns correspond
    to paired reads, a kmer length (k), and a distance between paired reads (d).  
    Outputs a reconstructed DNA sequence using caluclation of a DeBruijn graph
    and Eularian path. Returns error message if no overlap exists.
    """
    eulpath, balanced = EulerianPaired(pairedreads)
    pathlist = []
    FirstPatterns = []
    SecondPatterns = []
    for string in eulpath:
        pathlist.append(string.split(" "))
    for pair in pathlist:
        FirstPatterns.append(pair[0])
        SecondPatterns.append(pair[1])
    PrefixString = GenomePathReconstruction(FirstPatterns, balanced)
    SuffixString = GenomePathReconstruction(SecondPatterns, balanced)
    for i in range(k+d+1,len(PrefixString)):
        if PrefixString[i] != SuffixString[i-k-d]:
            return "there is no string spelled by the gapped patterns"
    return PrefixString + SuffixString[-k-d:] 
    
def PairedNonBranchingContigs(dbdict, inoutdict, startendnodes):
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
    startnodes, endnodes, balanced = startendnodes
    def TrimNodes(remainingdegrees, remainingpaths, node):
        """Removes nodes and node counts from appropriate dictionaries"""
        remainingdegrees[node][1] -= 1
        if remainingdegrees[node][0] !=0: #if not a start node subtracts 1 from indegrees
            remainingdegrees[node][0] -= 1 
        if remainingdegrees[node][0] == 0 and remainingdegrees[node][1] == 0: #if no remaining in- or out- degrees delete nodes
            del remainingpaths[node] 
            del remainingdegrees[node]
        return remainingdegrees, remainingpaths
    contiglist1 = []
    contiglist2 = []
    lastedges = []
    remainingpaths = copy.deepcopy(dbdict)  #duplicated of dbdict and inout dict are made so nodes can be deleted as used in contigs.
    remainingdegrees = copy.deepcopy(inoutdict)
    numofrows = len(startnodes)*2
    numofcolumns = 1
    for key in inoutdict: #create an array as large as potential contig fragments, a column allotted for every
        if inoutdict[key][1] > 1:
            numofcolumns += 1
    contigarray = np.empty((numofrows, numofcolumns), dtype = object)
    for x in range(numofrows):
        for y in range(numofcolumns):
            contigarray[x,y] = []
    starts = copy.deepcopy(startnodes)
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
        if node in remainingpaths: #if an internal node already used skips it
            nextnode = random.choice(remainingpaths[node]) #add next node by random choice if muliple outgoing nodes exist
            remainingpaths[node].remove(nextnode) #remove selected outgoing node from db graph list
            remainingdegrees, remainingpaths = TrimNodes(remainingdegrees, remainingpaths, node)
            node = str(nextnode)
            kmerlist = node.split(' ')
            nonbranchingpaths1.append(kmerlist[0][-1])
            nonbranchingpaths2.append(kmerlist[1][-1])
            while remainingdegrees[node][1] ==1 and remainingdegrees[node][0]== 1: #while encountering edges with 1 in- and 1 out- degree           
                nextnode = random.choice(remainingpaths[node])
                remainingpaths[node].remove(nextnode)
                remainingdegrees, remainingpaths = TrimNodes(remainingdegrees, remainingpaths, node)
                node = str(nextnode)
                kmerlist = node.split(' ')
                nonbranchingpaths1.append(kmerlist[0][-1])
                nonbranchingpaths2.append(kmerlist[1][-1])
        if node not in lastedges and node in remainingpaths:
            if remainingdegrees[node][1] !=0: 
                lastedges.append(node) #keep last edges to start next contig
            elif remainingdegrees[node][1] == 0: #delete end node once reached
                del remainingdegrees[node]
        if lastedges != []:        
            node = lastedges[0] #usefirst entry in list
        contiglist1.append(''.join(nonbranchingpaths1))
        contiglist2.append(''.join(nonbranchingpaths2))
        contigarray[i,j].append(''.join(nonbranchingpaths1))
        contigarray[i+1,j].append(''.join(nonbranchingpaths2))
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
            try: #if node has been fully visited delete from last edges
                remainingpaths[node]
            except NameError:
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
        if remainingpaths == {} and contig1 not in contiglist1:  #adds last fragment to list
            contiglist1.append(contig1)
            contigarray[i,j].append(contig1)
        if remainingpaths == {} and contig2 not in contiglist2:  #adds last fragment to list
            contiglist1.append(contig1)
            contigarray[i+1,j].append(contig2)
    return contiglist1, contiglist2, contigarray



def AssemblyPairedNonBranchingContigs(contigarray, k, d):
    """Inputs a contig array from non branching contigs, a kmerlength (k), and a paired read 
    distance (d).  Condenses paired reads by the k+d overlap between them. 
    Outputs a numpy array that contains non-branching contigs where each row in the 
    array represents a branched contig, the columns [0,...n] order the contigs 
    from beginning to end each as a list where adjacent columns represent adjacent 
    contigs. Branched contigs where multiple possibilities exist are listed in 
    the same list/column.
    """
    r, c = contigarray.shape
    assembledcontigs = np.empty((int(r/2), c), dtype = object) #create empty array with a list for each entry and half the rows of contigarray
    for x in range(int(r/2)):
        for y in range(c):
            assembledcontigs[x,y] = []
    i = 0
    l = 0
    while i < r:   #combine contigs from corresponding rows using an overlap
        for j in range(c):
            prefixlist = contigarray[i,j]
            suffixlist = contigarray[i+1,j]
            overlap = 1
            for m in range(len(prefixlist)):
                prefix = prefixlist[m]
                suffix = suffixlist[m]
                if prefix[k+d:]!= suffix[:-k-d]: #check if overlap betweeen paired strings matches
                        overlap = 0
                if overlap:
                    assembledcontigs[l,j].append(prefix + suffix[k+d:])
                else:
                    assembledcontigs[l,j].append("n/a")     
        i += 2
        l += 1
    return assembledcontigs

        
        
    
    
##TEST DATA     
#pr_in = """GAGA|TTGA
#TCGT|GATG
#CGTG|ATGT
#TGGT|TGAG
#GTGA|TGTT
#GTGG|GTGA
#TGAG|GTTG
#GGTC|GAGA
#GTCG|AGAT"""
#pr_in = pr_in.split("\n")
#
#k_in = 4
#d_in = 2
#pr_in = ReadPairsConvert(pr_in)
##Correct out = GTGGTCGTGAGATGTTGA
    
#pr_in = """ACC|ATA
#ACT|ATT
#ATA|TGA
#ATT|TGA
#CAC|GAT
#CCG|TAC
#CGA|ACT
#CTG|AGC
#CTG|TTC
#GAA|CTT
#GAT|CTG
#GAT|CTG
#TAC|GAT
#TCT|AAG
#TGA|GCT
#TGA|TCT
#TTC|GAA"""
#pr_in = pr_in.split("\n")

#k_in = 3
#d_in = 1
#pr_in = ReadPairsConvert(pr_in)
#
#print(pr_in)
#dbgraph = DeBruijnPaired(pr_in)
#print(dbgraph)
#inoutdict = InOutDict(dbgraph)
#print(InOutDict(dbgraph))
#print(StartEndNodes(inoutdict))

#
#print(EulerianPaired(pr_in))
#print(SequencefromPairedReads(pr_in, k_in, d_in))
#contig1, contig2, contigarray = PairedNonBranchingContigsPairedReads(dbgraph, inoutdict, StartEndNodes(inoutdict))
#print(PairedNonBranchingContigs(dbgraph, inoutdict, StartEndNodes(inoutdict)))
#print(AssemblyPairedNonBranchingContigs(contigarray, k_in, d_in))

