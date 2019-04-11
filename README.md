BioinfoTools
========

This collection of algorithms, scripts, and reference material serve as an organized repository to allow future access, use, and application of its contents. Python algorithms were written with the guidance of online resources such as the Bioinformatics Specialization on Coursera (Pavel Pevzner, Philip Compeau), Rosalind, the various other MOOC courses.  Statistics and probability reference guides largely represent material derived from the HarvardX Data Science courses on edX (Rafael Irizarry). 

## Table of Contents

1. [Bioinformatic Algorithms](#bioinformatic-algorithms)
    1. [Finding kmer Motifs in a DNA Sequence](#finding-kmer-motifs-in-a-dna-sequence) 
    2. [Finding kmer Motifs between DNA Sequences](#finding-kmer-motifs-between-dna-sequences)
    3. [DNA Sequencing](#dna-sequencing)
    4. [Peptide Sequencing](#peptide-sequencing)
    5. [Sequence Alignment](#sequence-alignment)
2. [Bioinformatic Tools](#bioinformatic-tools)
    1. [Multiple Sequence Alignment](#bioinformatic-tools)
    2. [Clustering and PCA](#bioinformatic-tools)
3. [Statistics and Probability Reference](#statistics-and-probability-reference)
    1. [Discrete Probability Reference](#statistics-and-probability-reference)
    2. [Distribution Functions Reference](#statistics-and-probability-reference)
    3. [Inference and Modeling Reference](#statistics-and-probability-reference)
4. [Basic Functions](#basic-functions)


## **Bioinformatic Algorithms**
--------


### Finding kmer Motifs in a DNA Sequence
---------

Kmers are short regions of DNA that have utility in finding unique patterns in DNA.  They can indicate transcription factor binding sites, polymerase binding sites, and a variety of regulatory elements. This collection of functions and programs find kmers within a genome or region of DNA.

**[Brute Force Motif Enumeration](https://github.com/jmmichaud/BioinfoTools/blob/master/BruteMotifEnumeration.py)** - PYTHON 3.6. Inputs a list of DNA fragments as a list of strings (dna), a kmerlength (k), and an allowed number of mismatches (d). Outputs a list of highest occurring motif or motifs of length k found between the DNA fragments with at most d mismatches for each occurrence. It employs an algorithm to determine all motifs of length k with at most d mismatches for each DNA fragment and then finds common motifs between the motif collections for each DNA fragment. This brute force algorithm is thorough but computationally expensive. It is only practical for short DNA lengths with low values of allowed mismatches.  

**[Find DnaA Box Within a Genome](https://github.com/jmmichaud/BioinfoTools/blob/master/FindDNAboxfromGenome.py)** - PYTHON 3.6. A small program that inputs a genome and finds candidates for DnaA box. G-C skew is used to identify region that contains origin of replication (ORI).  A 500 bp window identified by the skew or a user defined values are used to find highest occurring 9-mer motifs with at most 1 mismatch. The algorithm employs quaternary encoding and frequency arrays to improve computational efficiency allowing for location of larger motifs and more allowed mismatches.  The code can easily be altered to find motifs of different length and mismatches in any defined region of DNA within a genome.

**[Finding Kmers in Clumps](https://github.com/jmmichaud/BioinfoTools/blob/master/FindingKmersInClumps.py)** - PYTHON 3.6. Specific motifs are generally grouped or found in clumps or short regions within a genome. This collection of functions inputs a DNA sequence or genome as a string, a kmer length (k), a window or clump length (L), and a minimum number of occurrences (t). Outputs all kmers that occur within L-length windows that occur a minimum of t times within a L length window within the sequence.  Utilizes an efficient frequency array dictionary of the first L-length window as a basis and alters the dictionary based upon sliding one base at a time and changing the frequency of only the first and last kmers in the window. This algorithm can also be altered to give kmers with and number of occurrences and positions within the genome.


### Finding kmer Motifs between DNA Sequences
---------

Kmer motifs can also be shared between different homologs within different species and so it is useful to look at several DNA regions simultaneously.  It is also useful to look at multiple site within the same genome as distinct to identify kmers from different regions in the DNA.

**[Median Motifs Functions](https://github.com/jmmichaud/BioinfoTools/blob/master/MedianMotifFunctions.py)** - PYTHON 3.6. A collection of functions that demonstrate the development of an efficient method to identify motifs in different sequences of DNA.  The functions generally input a list of DNA strings and a kmer or motif length to find within the DNA strings.  The found kmers are sequences that generate a minimum score based on hamming distance. The approaches include brute force (exhaustive search), greedy (using probability matrices), randomization (random selection of intial kmers and probability), and Gibbs (Gibbs Sampling).

**[Gibbs Median Motifs](https://github.com/jmmichaud/BioinfoTools/blob/master/GibbsMedianMotifs.py)** - PYTHON 3.6.  A subset of functions from 'Median Motif Functions' that provide an efficient and optimized algorithm to find shared motifs between discrete DNA sequences.  The main function inputs a list of DNA strings, a kmer length, and DNA list length, and a number of iterations (N).  It outputs a list of best median motifs for each of the DNA strings in the given list minimizing the hamming distance between them.  It gives the GibbsMedianMotifs() 25 random starts and returns the best scoring collection of median strings. This algorithm can be further optimized by increasing or decreasing either the number of starts or N, the number of iterations.  This efficient algorithm is based upon Gibbs Sampling (Bayesian Inference) combining the advantages of weighted probability and randomness.  Using weighted random selection improves the accuracy and efficiency of this approach over a completely random selection approach.

**[Find Transcription Factor Binding Sites](https://github.com/jmmichaud/BioinfoTools/blob/master/FindTFBindingSites.py)** - PYTHON 3.6.  This code utilizes a modified version of the Gibbs Median Motifs algorithm to find transcription factor binding sites.  Using a list of DNA strings (dnalist), a kmer length (k), and DNA list length (t), a number of iterations (N), and a number of starts (M) it finds median strings in each DNA region.  This shows and example to find transcription factor binding sites in Mycobacterium tuberculosis.


### DNA Sequencing
---------

**[DNA Sequencing Functions](https://github.com/jmmichaud/BioinfoTools/blob/master/DNASequencingFunctions.py)** - PYTHON 3.6.  A collection of functions that input single or paired reads that assemble reads as intact sequences or non-branched contigs (largest runs of DNA that can be assembled without including paths with multiple possible variations). Assembly is guided through the generation of Debruijn Graphs and Eulerian Paths.

**[DNA Sequencing Single Reads](https://github.com/jmmichaud/BioinfoTools/blob/master/DNASequencingSingleReads.py)** - PYTHON 3.6.  Defines a class of single read DNA sequencing data and associated functions that assemble single read data into an intact assembly or non-branching contigs as data dictates.  Assembly is guided through the generation of Debruijn Graphs and Eulerian Paths.

**[DNA Sequencing Paired Reads](https://github.com/jmmichaud/BioinfoTools/blob/master/DNASequencingPairedReads.py)** - PYTHON 3.6.  Defines a class of paired read DNA sequencing data and associated functions that assemble paired read data into an intact assembly or non-branching contigs as data dictates.  Assembly is guided through the generation of Debruijn Graphs and Eulerian Paths.

### Peptide Sequencing
---------

**[Peptide Sequencing Functions](https://github.com/jmmichaud/BioinfoTools/blob/master/PeptideSequencingFunctions.py)** - PYTHON 3.6.  A collection of functions that input mass spectrometry measurements of unknown peptides and generate linear and cyclopeptide sequences.  The approaches generally entail building peptides by mass one amino acid at a time and trimming to eliminate possibilities not compatible with provided spectra at each stage to minimize run time.  Most function utilize normalized integer mass data.  An attempt is made to sequence a cyclic antibiotic using real mass spectrometry data.

**[Cyclopeptide Sequencing](https://github.com/jmmichaud/BioinfoTools/blob/master/CyclopeptideSequencing.py)** - PYTHON 3.6.  - A subset of 'Peptide Sequencing Functions' that optimize cyclopeptide sequencing.  The optimized Cyclopeptide Sequencing function, inputs two integers (M and N) and a mass spectrum of an unknown cyclopeptide as list of integer masses.  Outputs all lead peptides (as a list of masses and a list of in one-letter AA code) of a maximal score using N as a cutoff for a trim function including ties. Uses an AA list generated from processing of the spectrum and returning M top scoring masses between 57 - 200 including ties allowing for non-standard amino acids. It takes advantage of both using a leader board to keep N top scoring peptides and scoring algorithms that compare the spectrum to each generated spectra of peptides from the leaderboard to generate a leader peptide list where each peptide is maximal achieved score. 

### Sequence Alignment
---------

**[Sequence Alignment Functions](https://github.com/jmmichaud/BioinfoTools/blob/master/SequenceAlignmentFunctions.py)** - PYTHON 3.6. A collection of functions that input two either DNA or peptide sequences and generate alignments that maximize a score or path length using a Manhattan Tourist algorithm approach.  This method generally employs the creation of longest path matrices that use precursor nodes to calculate a maximum score and this is used to generate a "backtrack" matrix that uses symbols and allow for reverse construction of alignment sequences. Different approaches to matrix initialization enable both local and global alignments as well as alignments between sequences of differing length, overlapping sequences, and subsequences within provided sequences. The incorporation and comparison of both opening (introducing an insertion/deletion) and extension (continuing an insertion or deletion) allow for strong global alignments.

**[Sequence Alignment Functions Expanded](https://github.com/jmmichaud/BioinfoTools/blob/master/SequenceAlignmentFunctionsExpanded.py)** - PYTHON 3.6.  A larger set of functions that includes all of the functions in 'Sequence Alignment Functions' and several more that follow the development of the global alignment algorithm as well as a function to calculate the number of changes needed to make one sequence into another provided sequence.

## **Bioinformatic Tools**
--------

**[Multiple Sequence Alignment](https://github.com/jmmichaud/BioinfoTools/blob/master/MultipleSeqAlignment.md)** - R.  This R notebook demonstrates several multiple sequence alignments on 16S sequences using the DECIPHER package in R.  DECIPHER like many other aligners maximizes a score that combines structural and evolutionary alignment. It performs iterative sequence alignment of multiple sequences by first aligning two sequences and adding in subsequent sequences one at a time. There are other functions to refine alignments such as aligning DNA to RNA, alignments on very large datasets, staggered alignment that ignore non-homologous regions, and adjustments that shift gaps to allow better alignment. It reads in fasta file of several 16 sequences and uses AlignSeqs() and AdjustALignment() with default settings to align them.  AdjustAlign() allows for gaps to be shifted to improve alignments.  This was useful to align 16S sequences of different lengths and different coverage.

**[Clustering and PCA](https://github.com/jmmichaud/BioinfoTools/blob/master/ClusteringPCA.md)** - R.  This R notebook describes how to perform hierarchal clustering to create dendrograms and how to perform simple principal component analysis (PCA) in R using multiple methods.  Clustering is applied to differential gene expression (DGE) files.  PCA is performed on storm data with multiple parameters. 


## **Statistics and Probability Reference**
--------

The R notebooks in this section are reference material for statistic and probability concepts, their application, and the functions in R useful to their application. Most examples derive from Harvard Data Science course series on edX (PH125). The examples are to serve as a reference only and examples are either taken directly or adapted from examples in these courses.

 **[Discrete Probability Reference](https://github.com/jmmichaud/BioinfoTools/blob/master/DiscreteProbabilityReference.md)** - R. This is a brief guide to discrete probability in R. It provides background on the statical methods covered, pertinent R functions, and examples of their application.  It covers permutations, combinations, and the use of monte carlo simulations to make predictions. 

 **[Distribution Functions Reference](https://github.com/jmmichaud/BioinfoTools/blob/master/DistributionFunctionsReference.md)** - R. This is a brief guide to distribution functions and their application to continuous distributions in R. It provides background on the statical methods covered, pertinent R functions, and examples of their application.  It covers culmulative distribution, the central limit theorem, and law of large numbers applied to determining probability of events and sample modeling. The examples examine primarily leveraging a normal distribution even for variables that are not normally distributed and explain the 4 base R functions that are useful for a variety of distributions which are also referenced. The examples include modeling the behavior of a roulette wheel to evaluate potential profit and examining mortgage interest and foreclosure rates to increase the probability of profit.

 **[Inference and Modeling Reference](https://github.com/jmmichaud/BioinfoTools/blob/master/InferenceModelingReference.md)** - R. This is a brief guide to inference and modeling. It provides background on the statical methods covered, pertinent R functions, and examples of their application.  It covers p-value, spread, margin of error, and confidence intervals. These parameters are applied to look at samples to estimate values of the population the derive from and to evaluate the reliability and significance estimates. How to perform monte carlo simulations to make estimates about a population are also examined. The examples included apply these concepts to polling statistics.

## Basic Functions
--------

A repository of useful functions that are utilized in a variety of computation tasks.


**[Count DNA nucleotides](https://github.com/jmmichaud/BioinfoTools/blob/master/CountDNAnulcleotides.py)** - PYTHON 3.6.  Inputs a string of DNA.  Outputs the number of each nucleotide
    ACGT respectively as a space separated string of integers.
    
**[Reverse Complement](https://github.com/jmmichaud/BioinfoTools/blob/master/ReverseComplement.py)** - PYTHON 3.6.  Inputs a DNA string.  Outputs the reverse complement of the DNA strand.

**[Hamming Distance](https://github.com/jmmichaud/BioinfoTools/blob/master/HammingDistance.py)** - PYTHON 3.6.  Inputs two equal length DNA strings. Outputs the hamming distance between
    them as an integer.  Each mismatch is a score of 1.
    
**[GCSkew Minimum and Maximum](https://github.com/jmmichaud/blob/master/BioinfoTools/GCSkewMinMax.py)** - PYTHON 3.6.  Inputs DNA as a text string.  Functions can output either the positions of minimum or maximum GC SKew score. The GC Skew is calculated as the running difference between G and C bases with each G base contributing +1 and each C, -1. 

**[Kmer Neighbors](https://github.com/jmmichaud/BioinfoTools/blob/master/KmerNeighbors.py)** - PYTHON 3.6.  Inputs a DNA pattern as a text string and a number of allowed mismatches, d.  Outputs all permutations of the pattern than have no more than d mismatches. 

**[Frequent Kmers](https://github.com/jmmichaud/BioinfoTools/blob/master/FrequentKmers.py)** - PYTHON 3.6. Inputs a text string and an integer kmer length, k. Outputs the most most frequent kmers of length k within a given text using an algorithm to count the occurrences of kmers or patterns within the text. Will return ties (all kmers that are repeated the discovered maximum occurrences.) Output is a list of strings.

**[Quarternary Encoding](https://github.com/jmmichaud/BioinfoTools/blob/master/QuarternaryEncoding.py)** - PYTHON 3.6.  Two functions that either input a DNA sequence and output an integer representation using quarternary encoding or input the integer representation and a kmer or sequence length (k) and outputs the DNA sequence.  The quartnary encoding works the same as binary except each DNA base represents a digit, 'A' = 0, 'C' = 1, 'G'  = 2, 'T' = 3.

**[Short Frequency Array](https://github.com/jmmichaud/BioinfoTools/blob/master/ShortFrequencyArray.py)** - PYTHON 3.6.  Inputs a string of text and an integer kmer length (k).  Generates a frequency array as a dictionary where the key is the kmer and the entry is the count of that kmer in the text.  Avoids generation of every possible kmer of length k in the frequency array and only creates entries for kmers that exist.

**[Count Profile Consensus](https://github.com/jmmichaud/BioinfoTools/blob/master/CountProfileConsensus.py)** - PYTHON 3.6. Base functions that input a list of DNA text strings and a kmer or motif length. Together they compute counts of each DNA base at a given position of k, the probability of a given base at a given position of k, and a k-length consensus sequence that represents the most probably shared motif between the DNA strings using a probability matrix.  These functions are useful for algorithms that find a consensus sequence and closest shared motifs or median strings of a group of DNA sequences and multiple sequence alignments.

**[DNA Deconvolution Reconstruction](https://github.com/jmmichaud/BioinfoTools/blob/master/DNADeconvolutionReconstruction.py)** - PYTHON 3.6. Two opposing functions where one take a DNA string and deconvolute it into its component kmers as a list of k-length strings such that each string has a k-1 base overlap with at least one other DNA string and they are ordered in the list according to their overlap such that the DNA string [i][1:] matches DNA string [i+1][:k]. The second function takes such as list and reconstructs the DNA string.

**[Dynamic Programing Change](https://github.com/jmmichaud/BioinfoTools/blob/master/DynamicProgramingChange.py)** - PYTHON 3.6. A set of functions that input a total money amount and a list of coins and generate the smallest amount of coins possible to sum to the total. The outputs are the minimum number of coins required and a list of the coins used to reach that total.  It is a simple demonstration of the use of dynamic programming to improve computational efficiency.

**[Manhattan Tourist](https://github.com/jmmichaud/BioinfoTools/blob/master/ManhattanTourist.py)** - PYTHON 3.6. The Manhattan Tourist problem describes an optimal one-way path (with no overlap or backtracking) that a tourist in Manhattan should take to visit as many tourist attractions as possible. Such calculations are useful to a variety of computational problems. This function calculates the length of the longest path possible when the dimensions are given ("the street grid") and the number of attractions on each vertical and horizontal block. The function inputs grid dimensions n (vertical) and m (horizontal), as well as down(n x (m+1)) and right(m x (n+1)) matrices as numpy arrays. Creates a directed acyclic graph represented by a matrix of n+1(rows) by m+1(columns) of nodes denoting longest path lengths using the weight values in the down and right matrices. Outputs the length of longest path.

**[Remove List Duplicates](https://github.com/jmmichaud/BioinfoTools/blob/master/RemoveListDuplicates.py)** - PYTHON 3.6.  Short function to remove duplicate lists from a list of lists.

**[Useful DNA RNA Protein Functions](https://github.com/jmmichaud/BioinfoTools/blob/master/UsefulDNARNAProteinFunctions.py)** - PYTHON 3.6.  A short list of function that perform conversions from DNA to RNA and vice versa, protein translation from RNA, and a function that takes a short DNA sequence and find encoding substrings within a provided protein sequence.
    






