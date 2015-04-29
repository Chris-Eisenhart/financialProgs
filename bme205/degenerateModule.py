#!/usr/bin/env python2.7SD
# degenerateModule.py
# Chris Eisenhart  ceisenha@ucsc.edu/ceisenhart@soe.ucsc.edu
# 12.12.2014 - 12.16.2014
"""
getCodons(degenCodon):
	Takes in a degenerate codon as a string and returns all
	possible nucleotide codons that the degenerate codon
	could produce. 
codonListToAAList(codonList, aaTable):
	Takes in a list of nucleotide codons and a dict object 
	that maps nucleotides to amino acids. The nucleotides are
	translated and the occurence of each amino acid is counted. A 
	collections.Counter that maps amino acids to counts is returned.
imbalance(aaList):
	Takes in a collections.Counter object that maps amino acids
	to counts and returns a float. 
processAllCodons(URL):
	Takes in a URL as a string. Generates all possible degenerate 
	codons, for each degenerate codon all possible amino acids 
	are determined. For each amino acid set the 'best' degenerate 
	codon is determined by taking the degenerate codon with the
	lowest imbalance score and the highest frequency. For each
	amino acid set a line is printed, the line has four columns
	the first is the amino acid set, the second is the 'best' 
	degenerate codon, the third is the minimal imbalance score 
	for the amino acid list, and the fourth is the frequency 
	associated with the 'best' degenerate codon.
readCodonTable(URL):
	Takes in a URL as a string. Returns two hash tables, the
	first is a collections.Counter object that maps nucleotide
	codons to their count. The second is a dict object that maps
	nucleotide codons to their correpsonding amino acid single
	letter abreviation. 
"""
from __future__ import print_function  
import sys, operator, fileinput, collections, urllib
import string, os.path, re, argparse, bisect
import gzip, math, random, itertools

# Alphabet corresponds to all possible degenerate nucleotides.
# degenNucTable has a single entry for each degenerate nucleotide,
# mapping the degenerate nucleotide to the basic nucelotides (ATCG)
# it represents. For example degenNucTable['M'] = AC.  
alphabet = "ACGTRYKMSWBDHVN"
degenNucTable = dict()
degenNucTable['A'] = ['A']
degenNucTable['C'] = ['C']
degenNucTable['G'] = ['G']
degenNucTable['T'] = ['T']
degenNucTable['R'] = ['GA']
degenNucTable['Y'] = ['TC']
degenNucTable['K'] = ['GT']
degenNucTable['M'] = ['AC']
degenNucTable['S'] = ['GC']
degenNucTable['W'] = ['AT']
degenNucTable['B'] = ['GTC']
degenNucTable['D'] = ['GAT']
degenNucTable['H'] = ['ACT']
degenNucTable['V'] = ['GCA']
degenNucTable['N'] = ['AGCT']

def getCodons(degenCodon):
    """
    INPUT:
	degenCodon - A string; it should be three characters long 
	    with characters in the accepted extended alphabet. 
    OUTPUT:
	result - A list of strings, each string is three characters long
	    and contains only the characters A, C, T and G. 
    """
    # Result will hold all possible nucleotide codons as strings. 
    result = []
    # Turn degenCodon into a list to allow direct indexing 
    degenCodon = list(degenCodon)
    # unprocessedResult is a list of tuples. Each tuple has three elements, each element corresponds to a 
    # nucleotide. For example degenCodon = [ [(A, T, G)], [(G, C, T)] ... ]
    unprocessedResult = itertools.product("".join(degenNucTable[degenCodon[0]]), "".join(degenNucTable[degenCodon[1]]),
			 "".join(degenNucTable[degenCodon[2]]))
    # Remove tuple punctuation, appending the codons as strings to result. 
    for item in unprocessedResult:
	result.append("".join(item))
    return result
    
def codonListToAAList(codonList, aaTable):
    """
    INPUT:
	codonList - A collections.Counter object. The keys are nucleotide
	    codons and the values are the corresponding count. For example
	    codonList["ATG"] = 10948. 
	aaTable - A collections.Counter object. The keys are nucleotides and
	    the values are amino acid single letter abreviations. For example
	    aaTable["ATG"] = 'M'
    OUTPUT:
	result - A collections.Counter object. The keys are amino acid 
	    single letter abreviations and the values are the number 
	    of times the amino acid occured in codonList. For example
	    result["M"] = 10. 
    For each codon in codonList the corresponding amino acid is determined.
    The number of times each amino acid occures is counted.
    """
    # result will hold the amino acids and their counts. For example 
    # result['M'] = 10, result['W'] = 4. 
    result = collections.Counter()
    for codon in codonList:
	for codonToAA, AA in aaTable.iteritems(): 
	    if codon == codonToAA:
		result[AA] += 1
    return result

def imbalance(aaList):
    """
    INPUT:
	aaList - A collections.Counter object, the values are expected
	    to be integer. 
    OUTPUT:
	result - A float, the imbalance score as defined by Kevin Karplus.
    AaList is a collections.Counter object that stores the amino acid and 
    the number of times it was seen. For example aaList['M'] = 10.  
    """
    counts = aaList.values()
    return (max(counts) - min(counts))/float(sum(counts))

def processAllCodons(URL):
    """
    INPUT:
	none
    OUTPUT:
	none
    Generates all possible degenerate codons. For each codon three values are calculated
    all possible amino acids generated, the average frequency, and the imbalance score. 
    For each set of amino acids there are several degenerate codons, the best degenerate
    codon for each amino acid set is determined and printed with its corresponding values. 
    These rows are printed to sys.stdout and generate a histogram. The columns are 
    amino acid list, degenerate codon, imbalance score, and frequency. 
    """
    # allDegenCodons is a list of all possible degenerate codons, for example
    # NNN and AAA will both be present in this list.
    allDegenCodons = itertools.product(alphabet, repeat = 3)
    # Tables is a tuple with two elements, both elements are hash tabls. The first
    # is a collections.Counter object which holds the three letter codon in 
    # basic nucleotides (A, C, T, G). For each codon there is a number associated
    # this number corresponds to the number of times the codon is seen. For example
    # tables[0]["ATG"]= 103209.  The second table uses amino acid single letter 
    # abreviations as keys, and the nucleotide associated with the amino acid as a value.
    # For example tables[1]["M"] = ATG.  
    tables = readCodonTable(URL)
    # This will hold the results histogram. There are four columns, the first is the
    # list of amino acids, the second is the degenerate codon responsible for that list
    # the third is the imbalance score, the fourth is the frequency of the degenerate codon. 
    resultTable = []
    # Used for normalizing, this is the sum of all values in the codon-count table. 
    totalValues = sum(tables[0].values())
    # freqTable will hold the codon and its corresponding probability. Go through 
    # the codon count table normalizing each value. 
    freqTable = collections.Counter()
    for key, value in tables[0].iteritems():
	freqTable[key] = float(value/float(totalValues))
    
    # For each degenerate codon three histogram columns values must be made. 
    # All possible amino acids encoded by the degenerate codon, the imbalance
    # of the degenerate codon and the average frequency of the degenerate codon. 
    for degenCodon in allDegenCodons:	
	# Get all possible nucleotide codons the degenerate codon produces
	codonList = getCodons("".join(degenCodon))
	# Translate these codons into amino acids and count the occurences.
	# For example if ten codons coded for M then AAList["M"]=10
	AAList = codonListToAAList(codonList, tables[1])
	# When calculating the average frequency for a degenerate codon 
	# some codons will be coded for multiple times. The average frequency 
	# is not weighted here, so duplicate nucleotide values are discarded. 
	seenCodons = []
	# The probability of each codon is appended to this value. After
	# all codons have been processed the value is averaged. 
	frequency = 0
	for codon in codonList:
	    #if codon in seenCodons: continue
	    seenCodons.append(codon)
	    frequency += freqTable[codon]
	frequency = float(frequency)/float(len(seenCodons))
	# Append all the amino acids encoded for by the single degenerate
	# nucleotide to allAA. 
	allAA = []
	for key, value in AAList.iteritems():
	    if key not in allAA: 
		allAA.append(key)
        allAA = sorted(allAA)
	imbScore = imbalance(AAList)
	bestAAList = "".join(allAA)
	bestDegenCodon = "".join(degenCodon)
	bestImbScore = imbScore
	bestFrequencey = frequency
	# Append the histogram row to resultTable. 
	resultTable.append((bestAAList, bestImbScore, bestDegenCodon, bestFrequencey))

    # Sort the result table by the amino acid list, this clusters 
    # all identical amino acid lists together. For each amino acid list the best
    # degenerate nucleotide is determined by frequency and imbalance scores. 
    # The best degenerate nucleotide in a given block will have the highest frequency. 
    sortedTable = sorted(resultTable, key = operator.itemgetter(0))
    # In order to differentiate blocks two internal structures are used. The first
    # is the boolean firstBlock, which indicates the start of the list. The second
    # is the currentBlock string.  CurrentBlock is set to the amino acid list using
    # the firstBlock boolean. For all iterations after this, if the currentBlock is
    # not the same as the amino acid set, the program prints the best degenerate
    # nucleotide histogram row from the current amino acid set. 
    firstBlock = True
    # bestFreq is used to determine the best frequency. For each block
    # bestFreq is reset to 0 then overwritten by the highest frequency. 
    # If a frequency is higher than highScore all the values that correspond
    # to that frequency are stored so that they can be printed when the block ends.  
    bestFreq = 0
    # The minimal imbalance for each block, it is overwritten by lower 
    # imbalances then returned and reset when the block ends.      
    minImbalance = 1
    # finalResult holds the data in four columns, AAlist is the amino acid
    # set, imbScore is the imbalance score, codon is the 'best' degenerate codon for the amino acid set, 
    # and frequency is the average frequency associated with the 'best' degenerate codon. 
    finalResult = []
    for AAlist, imbScore, codon, frequency in sortedTable:
	if firstBlock:
	    firstBlock = False
	    currentBlock = AAlist
	if AAlist != currentBlock:	
	    # This marks the end of the current block, the best row 
	    # is appended to the finalResult histogram and currentBlock,
	    # bestFreq and minImbalance are reset. 
	    currentBlock = AAlist
	    finalResult.append((bestAAList, bestImbScore, bestDegenCodon, bestFreq))
	    bestFreq = 0 
	    minImbalance = 1
	if imbScore <= minImbalance:
	 	# Only consider the rows with low imbalance
		minImablance = imbScore
		if (frequency > bestFreq):
		    # Store the row with the best frequency 
		    bestFreq = frequency 	
		    bestDegenCodon = codon
	            bestAAList = AAlist
	            bestImbScore = imbScore
    finalResult.append((bestAAList, bestImbScore, bestDegenCodon, bestFreq))
    # Sort the finalResult histogram by the amino acid list, this causes the final
    # output to be in alphabet order. 
    sortedFinalResult = sorted(finalResult, key = operator.itemgetter(3), reverse = True)
    for AAlist, imbScore, codon, frequency in sortedFinalResult:
	print (AAlist, imbScore, codon, frequency)

def readCodonTable(URL):
    """
    INPUT:
	inputFile - A file like object. 
    OUTPUT:
    The codon table must be in the format seen at http://www.kazusa.or.jp/codon
    /cgi-bin/showcodon.cgi?species=273057.  The format is codon, amino acid,
    fraction, frequency per thousand, number. The number element is assumed 
    to be in the format '( #number)'.
    """
    page = urllib.urlopen(URL)
    table = False
    # A trigger for identifying the table in the HTML code. The line
    # before the table starts will contain a <PRE> tag, this boolean
    # will keep track of whether the PRE tag has been seen. 
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
	    # After splitting the elements are grouped into six elements. 
	    # The first is codon, the second is the amino acid the third
	    # is the fraction, the fourth is the frequency per thousand 
	    # the fifth is '(' and the sixth is the raw number. Of these
	    # values only the codon, amino acid and raw number are needed.  
	    for item in splitLine:
	        itemNumber += 1
    		complement_table = string.maketrans("U", "T")
	        # Convert U to T so that the program corresponds to DNA. 
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
