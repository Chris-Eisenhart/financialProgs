#!/usr/bin/env python2.7
# reverseORFModule REDO
# Chris Eisenhart  ceisenha@ucsc.edu/ceisenhart@soe.ucsc.edu
# 11.15.2014-12.16.2014 
"""
This module contains functions intended for use with the 
python program reverseORF.  The following functions are
contained within this module. Note that these functions 
are all intended for use with Sulfolobus solfatariucs. 
Much of the code is generalized and will work with other 
species, but as a whole this code is designed for 
Sulfolobus solfataricus, and any other species information
will likely result in erroneous output. 

reverseComplement(string) - Returns the reverse compliment of 
	the string DNA. The string is inverted. 

readCodonTable (file) - Returns a collections.Counter() holding
	the codon values from input file. 

scanSequence (string) - Returns an integer that is 
	the length of the longest ORF seen in the input string. 

weightedRandom(collections.Counter) -  Returns a randomly 
	selected key from the hash table. The key is selected
	on a weighted random based on its value. For example
	a key with value 10 will be returned more than a 
	key wil value 1. 

getModel1Result(fastaFile) - Returns the probability of an
	ORF 388 or longer occurring. Calculated using GC content

getModel2Result(fastaFile) - Returns the probability of an ORF
	388 or longer occurring. Calculated using codon counts	
	from a similar organism.  

generateSequence3 (collections.Counter, int) - Returns a 
	randomly generated string. The string is constrained so 
	that it begins with a start codon, and all other 
	start/stop codons are removed.  

generateSequence4 (collections.Counter, string) - Returns a string 
	generated from a template string. For each AA on the
	template string a random codon which encodes for that
	AA is appended to the result string. 

getORFs (int, int, collections.Counter, collections.Counter,
			 int, file, file) - Generate a 
	set number of sequences (default 10,000), each 
	sequence is scanned for the longest ORF. The 
	number of times each ORF occures is stored in a hash.
	The hash is returned as a collections.Counter. 

printHistogram (collections.Counter, string, int, bool)-
	Prints a histogram to sys.stdout. The histogram has
	several statistics in the header.  

countKmers (collections.Counter, string, int) - Kmers of 
	the specified size are counted and stored in 
	a collections.Counter object. The updated object
	is returned. 
"""
from __future__ import print_function  
import sys, operator, fileinput, collections, urllib
import string, os.path, re, argparse, bisect
import gzip, math, random, itertools, fastaFastqParser

def reverseComplement(dna):
    """
    INPUT: 
	dna - A string.
    OUTPUT:	
	result - A string. 
    Returns a string that is the reverse-complement of  
    string "dna". 
    """
    complement_table = string.maketrans("ACGT", "TGCA")
    result = dna[::-1].translate(complement_table)
    return result 

def readCodonTable(URL):
    """
    INPUT:
	inputFile - A file like object. 
    OUTPUT:
	countsTable - A collections.Counter object that holds 
	    DNA codons as keys and an integer value that corresponds	
	    to the codons frequency. For example, countsTable['ATG'] = 10986
	aaTable- A dict object that holds all DNA codons as keys
	    and the corresponding single letter amino acid abreviation
	    as values. For example aaTable['ATG'] = M 
    Takes in a URL string and opens it as a file. The
    codon table is loaded into a collections.Counter()
    object and a dict object. Together these hold the codon, and 
    corresponding amino acid and counts.  The hash is returned after
    all inputFile lines are read. The codon table must
    be in the format seen at http://www.kazusa.or.jp/codon
    /cgi-bin/showcodon.cgi?species=273057.  The format is 
    codon, amino acid, fraction, frequency per thousand, 
    number. The number element is assumed to be in the 
    format '( #number)'.
    """
    # Open the URL 
    page = urllib.urlopen(URL)
    table = False
    # A trigger for identifying the table in the HTML code. The table will
    # start after a <PRE> tag, this boolean will be used to identify if 
    # The program is before or after the <PRE> tag. 
    countsTable = collections.Counter()
    # This will store the codon as a string and its observed
    # counts as an integer. For example countsTable['ATG'] = 10986
    aaTable = dict()
    # This will store the codon as a string and its amino acid as a string
    # For example aaTable['ATG'] = M
    for line in page:
	if line.startswith( "<PRE>"):
        # The table starts and table is switched to True
	    table = True
	    continue
	if line.startswith("</PRE>"): table = False
        # The table ends and table is switched back to False
        if table: 
        # At this point only lines of the table are being 
        # considered
	    splitLine = line.split()
	    codon = None
	    aa = ""
	    count = 0
	    itemNumber = 0
	    for item in splitLine:
	    # After splitting the elements are grouped by the 
	    # five columns and the '(' in the number element 
	    # ( number = '( ####)'). This gives a total of six 
	    # elements per group.  The first is the codon, the 
	    # second is the amino acid, and the 6th is the 
	    # raw counts. The raw counts has an erroneous ')'
	    # appended, so the last character is removed. 
	        itemNumber += 1
    		complement_table = string.maketrans("U", "T")
	        if itemNumber is 1: codon = item.translate(complement_table)
	        if itemNumber is 2: aa = item
	        if itemNumber is 6:
	        # The end of a group of elements, the completed 
	        # key value pairs are inserted into countsTable and
	        # aaTable 
		    count = int(item[:-1])
		    countsTable[codon] = count
		    aaTable[codon] = aa
	            itemNumber = 0
    return countsTable, aaTable
	
def scanSequence(sequence):
    """
    INPUT:
	sequence - A string. 
    OUTPUT:
        len(result) - An integer. 
    Takes in a DNA sequence as a string. The reverse compliment
    of the sequence is determined and the sequence is reversed. 
    All three reading frames are scanned over and the longest 
    ORF is identified. The length of this ORF is returned. 
    """
    # I worked on this function with Mary Lenore Pafford. 
    result = 1
    # The default ORF length, this value will be replaced
    # by ORF lengths from the sequence. 
    frame = 0
    # This will keep track of the current frame in the sequence, 
    # it will be reset after reaching 3.
    ORFs = collections.Counter()
    # Holds the frame number and the length of the current ORF
    for i in range(3):
        ORFs[i] = None
    # Initialize the three frames to none, signifying no ORF.
    stopCodons = {"TGA","TAG","TAA"}
    # Iterate over the sequence one time using an integer to 
    # keep track of the frame. If a start sequence is seen 
    # in a frame then that frame becomes a key is ORFs and the
    # length of the sequence is the value. If a stop codon is 
    # seen in a frame then the length of the ORF is compared
    # to the current longest ORF, if longer it replaces it.  
    for start in range(len(sequence)):
	if frame == 3: frame = 0
	if (sequence[start:start+3] == "ATG" and ORFs[frame]
					 is None):
	    ORFs[frame] = 0
	if (sequence[start:start+3] in stopCodons): 
	    if ORFs[frame] is not None: ORFs[frame] += 1
	    if ORFs[frame] > result: result = ORFs[frame] 
	    ORFs[frame] = None
	if ORFs[frame] is not None: 
	    ORFs[frame] += 1
	frame += 1
    # This addresses the ORFs that had start codons, but 
    # no stop codons. They have not been compared to the
    # result yet, so it is done here. If the ORF length is 
    # longer than the result it replaces it. 
    for DNAFrame, ORFLength in ORFs.iteritems():
	if ORFLength > result: result = ORFLength
    return result
	  

#Stack Overflow user:eumiro
def weightedRandom(weights):
    """
    INPUT:
	weights - A collections.Counter object. 
    OUTPUT:
	key - A string.
    Computes a weighted random variable using the 
    input hash table. The random variable is returned. 
    """
    number = random.random() * sum(weights.values())
    # Choose a random number from the cumulative total
    # codon counts. 
    for key, value in weights.iteritems():
    # Find and return the key that corresponds to the 
    # random value. 
	if number < value:
	    break
	number -= value
    return key

def getModel1Result(dnaFile):
    """
    INPUT:
    	dnaFile - An opened file like object. Must be in fasta/.fa format. 
    OUTPUT:
	result - A float that corresponds to the probability of an ORF of
	    388 occuring given the current model. 
	totalNucleotides - An integer, the total number of nucleotides seen. 
	ATG - A float, the probability of a start codon. 
	TAG + TGA + TAA - A float, the probability of a stop codon. 
	weights['C'] + weights['G'] - A float, the percent of the genome that 
	    is either G or C.  
    Calculates the probability of an ORF of length 388
    occurring in Sulfolobus solfataricus. This model uses 
    the GC content to do a direct calculation. 
    """
    weights = collections.Counter()
    totalNucleotides = 0 
    for fasta in fastaFastqParser.read_fasta(dnaFile, None, False, False): 
	# Count the kmers in the input fasta file, store them 
	# in weights. The codons are string keys, and the count is
	# value stored as an integer. 
	totalNucleotides += len(fasta.sequence)
        countKmers(weights, fasta.sequence, 1)
    number = sum(weights.values())
    for key, value in weights.iteritems():
	# Assign probabilities based on observed count
	weights[key] = float(value)/float(number)
    ATG= weights['A']*weights['G']*weights['T']
    TAG= weights['T']*weights['A']*weights['G']
    TGA= weights['T']*weights['G']*weights['A']
    TAA= weights['T']*weights['A']*weights['A']
    result = float(ATG*((1-TAG-TGA-TAA)**387))
    result = result * 2 * (number-387)
    return (result, 2 * totalNucleotides, ATG, (TAG + TGA +TAA), (weights['C']+weights['G']))
    
def getModel2Result (dnaFile):
    """
    INPUT:
    	dnaFile - An opened file like object. Must be in fasta/.fa format. 
    OUTPUT:
	result - A float that corresponds to the probability of an ORF of
	    388 occuring given the current model. 
	totalKmers - An integer, the total number of kmers seen. 
	weights['ATG'] - A float, the probability of a start codon. 
	weights['TAG'] + weights['TGA'] + weights['TAA'] - A float, the probability of a stop codon. 
    Calculates the probabillity of an ORF of length 388
    occuring in Sulfolobus solfataricus. This model uses 
    codon counts to do a direct calculation. 
    """
    weights = collections.Counter()
    totalKmers = 0
    for fasta in fastaFastqParser.read_fasta(dnaFile, None, False, False): 
	# Count the kmers in the input fasta file, store them 
	# in weights. The codons are string keys, and the count is
	# value stored as an integer. 
	totalKmers += len(fasta.sequence) - 2
        countKmers(weights, fasta.sequence, 3)
    number = sum(weights.values())
    # number is the sum of the kmer counts, the total number
    # of kmers seen. 
    for key, value in weights.iteritems():
	# Assign probabilities based on observed count
	weights[key] = float(value)/float(number)
    result = float(float(weights["ATG"])*((1-float(weights["TAG"])-float(weights["TGA"])-float(weights["TAA"]))**387)) 
    result = result * 2 * (number-387)
    return (result, 2 * totalKmers, weights['ATG'], (weights['TAG'] + weights['TGA'] +weights['TAA']))

 
def generateSequence3(kmerWeights, length):
    """
    INPUT:
	kmerWeights - A collections.Counter object. 
	length - An integer.
    OUTPUT:
	seq - A string.
    Generates a sequence of length 'length'. This sequence
    is generated by appending complete codons, therefore the
    sequence length must be divisible by 3. The sequence
    generated will be an ORF of length length, start and stop
    codons will not be generated with the exception of the
    first codon.  
    """
    seq = "ATG"
    kmerWeights["TAG"] = 0
    kmerWeights["TGA"] = 0
    kmerWeights["TAA"] = 0
    # Set the start and stop codon counts to 0, 
    # these codons will not be generated. 
    for i in range((length/3)-1):
        # The length is in nucleotides, it is converted
        # to amino acid length. A weighted random 
        # codon is appended to the sequence each iteration
	seq += weightedRandom(kmerWeights)
    # Reverse compliment the seqeunce. 
    result = reverseComplement(seq)
    return result

def generateSequence4(aaTable, codonTable, refSeq):
    """
    INPUT:
	aaTable - A collections.Counter object. 
	refSeq - A string.
    OUTPUT:
	seq - A string.
    Generates a sequence of length 3 * refSeq. This sequence
    is generated by appending complete codons, therefore the
    sequence length must be divisible by 3. The sequence
    generated will be the reverse compliment of a gene that 
    codes for protein sequence refSeq. 
    """
    seq = ""
    for char in refSeq:
    # Iterates over the reference sequence provided. Each
    # character represents an amino acid. For each amino acid
    # all the potential codons are determined and inserted
    # into a hash table, then one is chosen at random. This 
    # random codon is appended to the sequence. 
	codons = collections.Counter()
	# codons["M"] will return the counts of methionine. 
	for AAcodon, aminoAcid in aaTable.iteritems():
	    # Find all codons for the amino acid and load them into 
	    # the hash table.
	    if char == aminoAcid:
		# Find the correct count for the amino acid
	        for codon, count in codonTable.iteritems():
		    if codon == AAcodon: codons[AAcodon] = count
	seq += weightedRandom(codons)    
    # The sequence is reverse complimented and returned.
    result = reverseComplement(seq)
    return result


def getORFs(seqNumber, seqLength, codonTable, aaTable, model, dnaFile, proteinFile):
    """
    INPUT:
	seqNumber - An integer.
	seqLength - An integer.
	codonTable - A collections.Counter object.
	aaTable - A collections.Counter object. 
	model - An integer. 
    	dnaFile - An opened file like object. Must be in fasta/.fa format. 
	proteinFile - An opened file like object. Must be in fasta/.fa format. 
    OUTPUT:
	ORFs - A collections.Counter object.
    Generates seqNumber of sequences. The sequences
    are generated using one of four models, specified by
    the model option. 
    Each sequence is scanned for ORF and the 
    longest ORF is returned. The count of these
    returned ORFs is sorted in ORFs. ORFs 
    is returned using the ORF length as a key, and
    its occurrences as the value. 
    """
    ORFs = collections.Counter()
    # This will hold the ORF length as a key andt he number
    # of times an ORF of that length is seen as a value. 
    # The key is a string and the value is an int
    if model is None: model = 4
    if model is 4:
    # For model 4 a protein reference sequence is needed
        refSeq = ""
        for fasta in fastaFastqParser.read_fasta(proteinFile, None, False, False): 
    	    refSeq += fasta.sequence
    # Handle models 1 and 2, these require no simulation
    if model == 1:result = getModel1Result(dnaFile)
    if model == 2:result = getModel2Result(dnaFile)
    # generates a sequence using the specified model then
    # scans the sequence for an ORF, returning the length of 
    # the longest ORF as a integer. This integer becomes
    # a string key in ORFs, each time this int is returned
    # its value is incremented. 
    for i in range(seqNumber):
	if model == 3: 
	    sequence = generateSequence3(codonTable, seqLength)
	    ORFs[scanSequence(sequence)] += 1
	if model == 4: 
	    sequence = generateSequence4(aaTable, codonTable, refSeq)
	    ORFs[scanSequence(sequence)] += 1
    if model == 1: print ("Exp. number of ORFs 388 or longer"
			" based on gc content " + str(result))
    if model == 2: print ("Exp. number  of ORFs 388 or longer"
			" based on direct counts " + str(result))
    return ORFs

def printHistogram(ORFs, dnaFile, sequences, rOutput):
    """
    INPUT:
	ORFs - A collections.Counter object. 
    	dnaFile - An opened file like object. Must be in fasta/.fa format. 
	gC - A float. 
	bases - An integer.
	kmers - An integer. 
	sequences - An integer. 
	rOutput - Boolean.
    OUTPUT:
	None
    Prints a histogram to sys.stdout.  The histogram contains 
    basic statistics about Sulfolobus solfataricus. 
    The input hash table is sorted then printed.
    The hash table is printed to four columns (tab separated);
    the first is ORF length, the second is the number of 
    occurrences, the third is a calculated probability, and 
    the fourth is the probability of seeing
    an ORF that length or greater. If  rOutput is
    selected column headers are printed and the
    histogram header is omitted, additionally the columns are
    comma separated. 
    """
    dnaFile1 = gzip.GzipFile(dnaFile, "r")
    model1Result = getModel1Result(dnaFile1)
    dnaFile2 = gzip.GzipFile(dnaFile, "r")
    model2Result = getModel2Result(dnaFile2)
    if not rOutput:
        print ("# Reading from " + dnaFile)
        print ("# The double-stranded genome contains " + str(model1Result[1])  
		+ " bases, and " +str(model2Result[1]) + " 3-mers.")
        print ("# probablity of G+C is " + str(model1Result[4]))
        print ("# probablitiy of ATG start codon is " + str(model1Result[2]) + " based on"
		"GC content, " + str(model2Result[2]) + " on direct count")
        print ("# probablitiy of stop codon is " + str(model1Result[3]) +" based on"
		"GC content, "+ str(model2Result[3]) + " on direct count")
        print ("# Expected number of 388 or longer orfs based"
		" on GC frequency is " + str(model1Result[0]))
        print ("# Expected number of 388 or longer orfs based on"
		" 3-mer frequency is " + str(model2Result[0]))
    sortedORFs = sorted(ORFs.iteritems(), key = operator.itemgetter(0))
    start = 1
    if rOutput: print ("key,value,probability,pValue")
    # Print the histogram and calculate stats. The 
    # probability is calculated by dividing the 
    # number of times it was seen by the number of sequences
    # processed. The P-value is created by recursively 
    # subtracting probabilities from 1. 
    for key, value in sortedORFs:
	if rOutput:
	    print (str(key) + "," + str(value) + "," + 
		    str(float(value)/sequences) + "," +
		    str(start))
	else:
	    print (str(key) + "\t" + str(value) + "\t" + 
		    str(float(value)/sequences) + "\t" +
		    str(start))
	start = start-(float(value)/sequences)

def countKmers(counts, sequence, kmerSize):
    """
    INPUT:
	counts - A collections.Counter object
	sequence - A string. 
	kmerSize - An integer. 
    OUTPUT:
	counts - A collections.Counter object. 
    Counts kmers of kmerSize in the sequence. Each
    time a kmer is seen the count is incremented.
    """
    for start in range(len(sequence)):
        counts[sequence[start:start+kmerSize]] += 1
    return counts
