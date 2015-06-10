#!/usr/bin/env python2.7
# palindromeModule.py REDO 1
# Chris Eisenhart  ceisenha@ucsc.edu/ceisenhart@soe.ucsc.edu 
# 10.28.2014 - 11.29.2014 
"""
This program contains functions intended for use with the 
Python programs palindromes, . 
All functions for the two programs will be stored here. 
These functions are accessed in other python programs by 
prefixing this program's name, for example
palindromeModule.reverse_comp(counts, True) is a valid 
function call. 
fastFunctions was writted by Robert Calef

The following functions are contained in this module. 
reverse_comp(dna)
	Returns the reverse compliment of the string DNA
expected_value(palindrome,kmer_counts) 
	Returns the expected value of the palindrome
get_results_and_stats(inputList, min_k, max_k, verbose,
		 alphabet, allPalindromes, hypothesesCount)
	Returns the output table and some general stats. 
calculate_results(representedPalindromes, kmer_counts, 
	    	 totalCharCount, hypothesesCount)
	Creates the output table 
get_represented_palindromes(kmerCounts, allPalindromes)
	Returns a list of all kmers that are also palindromes
print_palindromes_and_stats(resultsAndStatsAndStuff, min_k,
		 max_k, hypothesesCount,cutoff)
        Print the table to stdout
get_all_palindromes_and_hypotheses_count(min_k,max_k)
	Returns all palindromes within the range 
get_palindromes(counts, alpha, k)
	Returns a hash of palindromes of length k.
count_kmers_and_totalChars(inputFile,  min_k, max_k, 
		alpha, verbose)
	Returns a hash of kmers and counts, and an int.
"""
import sys, operator, fileinput, string, os.path
import re, argparse, itertools, gzip
import collections, math, fastFunctions

def reverse_comp(dna):
    """
    INPUT: 
	dna - A string
    OUTPUT:	
	Returns a string that is the reverse-complement of  
    string "dna". 
    """
    complement_table = string.maketrans("ACGT", "TGCA")
    return dna[::-1].translate(complement_table)

def expected_value(palindrome,kmer_counts):
    """
    INPUT:
	palindrome - A string
	kmer_counts - A collections.Counter object
    OUTPUT:
	Returns a float that is the expected count of 
	palindrome. 

    The expected value is calculated by identifying
    the counts of substrings of palindrome. 
    These substrings are 2 to k, 1 to (k-1) and 2 to (k-1).
    The count of substring 2 to k is multiplied with the count
    of substring 1 to (k-1) then divided by substring  2 to
    (k-1).  This value is the expected value and is returned 
    as a float. 
    """
    firstCount = 0
    # Store the count of the kmer created by 
    # removing the last character of palindrome.
    secondCount = 0
    # Store the count of the kmer created by 
    # removing the first character of palindrome.
    if (len(palindrome) % 2) == 1:
       # Use [:-1] to remove the last character from a string.
       # Use [1:] to remove the first character from a string. 
       kmerList = list(palindrome)
       if kmerList[(len(palindrome)/2)] == "W":
           kmerList[(len(palindrome)/2)] = "A"
           firstCount+=float(kmer_counts["".join(kmerList)[:-1]])
           secondCount+=float(kmer_counts["".join(kmerList)[1:]])
           kmerList[(len(palindrome)/2)] = "T"
           firstCount+=float(kmer_counts["".join(kmerList)[:-1]])
           secondCount+=float(kmer_counts["".join(kmerList)[1:]])
       else:
           kmerList[(len(palindrome)/2)] = "C"
           firstCount+=float(kmer_counts["".join(kmerList)[:-1]])
           secondCount+=float(kmer_counts["".join(kmerList)[1:]])
           kmerList[(len(palindrome)/2)] = "G"
           firstCount+=float(kmer_counts["".join(kmerList)[:-1]])
           secondCount+=float(kmer_counts["".join(kmerList)[1:]])
    else:
        firstCount = float(kmer_counts[palindrome[:-1]])
        secondCount = float(kmer_counts[palindrome[1:]])
    thirdCount = float(kmer_counts[palindrome[1:-1]] + .0000000001)
    # Store the count of the kmer created by 
    # Since the thirdCount will be dividing the first two it is 
    # incremented by 1 e -10 to prevent division by 0. 
    # removing the first and last character of palindrome.
    return float(firstCount * secondCount / thirdCount)
    
def get_results_and_stats(inputList, min_k, max_k, verbose,
		 alphabet, allPalindromes, hypothesesCount):
    """
    INPUT:
	inputList - A string
	min_k - A string
	max_k - A string
	verbose - Boolean
	alphabet - A string
	allPalindromes - A collections.Counter object
	hypothesesCount - An integer
    OUTPUT:
	resultsAndStats - A list of tuples, each tuple
		holds a palindrome, its observed and
		expected counts, z-score, and E-value. 
	fileNames - A string 
	totalCharCount - An int

    This function iterates over all files, counting the number
     of times each palindrome is seen,  calculating the
     expected number of times a palindrome should be seen,
     the Z-score for each palindrome, and the E-value for 
     each palindrome.  These values are stored in a list of 
     tuples, where each tuple is storing five fields, 
     palindrome name, observed count, expected count, z-score
     and e-value.  In addition the file names, and total 
     character counts are returned.  Therefore the object
     returned for resultsAndStats is tuple with three
     elements. The first is the tupleList, the second is the
     file names, and the third is the total character count.  
    """
    resultsAndStats = []
    # The tuple list, will store palindrome and stats
    kmerCounts = collections.Counter()
    # The kmer's seen and their counts
    totalCharCount = 0
    # The total number of characters seen
    fileNames = ""
    # All file names
    representedPalindromes = collections.Counter()
    # All palindromes that are present in the files. 
    for possibleFile in inputList[1:]:
        if possibleFile.endswith('.fa') or possibleFile.endswith('.fa.gz'):
	    newFile = possibleFile
            if newFile.endswith(".gz"): 
                inputFile = gzip.GzipFile(newFile, "r")
                fileNames += newFile
            else: 
                inputFile = open(newFile, "r")
                fileNames += newFile
            fileNames += " "
            kmerCountsAndTotalChars=count_kmers_and_totalChars(
			inputFile, min_k - 2, max_k,
			 alphabet, verbose)
            kmerCounts += kmerCountsAndTotalChars[0]
            totalCharCount += kmerCountsAndTotalChars[1]
            representedPalindromes = get_represented_palindromes(
			kmerCounts, allPalindromes)
    resultsAndStats += calculate_results(
		representedPalindromes, kmerCounts,
		  totalCharCount, hypothesesCount)
    return resultsAndStats, fileNames, totalCharCount


def calculate_results(representedPalindromes, kmer_counts, 
	    	  totalCharCount, hypothesesCount):
    """
    INPUT:
	representedPalindromes - A collections.Counter() 
		object
	kmer_counts - A collections.Counter() object. 
	totalCharCount - An integer 
	hypothesesCount - An integer
    OUTPUT:
	tupleList - A list of tuples, each tuple
        contains the palindrome name, observed count,
        expected count, z-score and E-value. 

    Takes in a collections.Counter representedPalindromes,
    a collections.Counter kmer_counts,
    an int totalCharCount, and an int hypothesesCount.  For 
    each palindrome in representedPalindromes the expected
    value, z-Score, and e-value is calculated. These three
    values are coupled with the palindrome name and observed
    count, to create a tuple with five elements.  These 
    tuples are loaded into a list, which is returned after 
    all palindromes in representedPalindromes has been 
    iterated over. 
    """
    tupleList = []
    # The list of tuples. Each element in this list is a 
    # tuple, with five elements. The tuple elements
    # are palindrome name, observed count, expected count,
    # z-score and E-value. 
    for key, value in representedPalindromes.items():
        if (len(key)%2) == 1:
            keyList = list(key)
            if(keyList[len(key)/2]) in "ATCG": continue
	    # Throw out odd palindromes with middle characters
	    # ATCG, they are not desired.  
        expectedValue = expected_value(key,kmer_counts)
        standardDeviation = (math.sqrt(expectedValue *
		 (1 - (expectedValue/ totalCharCount))))
        zScore = ((kmer_counts[key] - expectedValue)
		/standardDeviation)
        if zScore <= 0:
           e_Value = ((math.erfc(-zScore/math.sqrt(2))/2)*
		2*hypothesesCount)
        else:
           e_Value = ((math.erfc(zScore/math.sqrt(2))/2)
		*2*hypothesesCount)
        tupleList.append((key, str(value), str(expectedValue),
		 zScore, e_Value))
    return tupleList


def get_represented_palindromes(kmerCounts, allPalindromes):
    """
    INPUT:
	kmerCounts- A collections.Counter() object 
	allPalindromes- A collections.Counter() object
    OUTPUT:
	represented-palindromes - A collections.Counter() 
    			object
    Takes in two collections.Counters.  kmerCounts holds 
    all kmers seen in the file. allPalindromes holds all 
    palindromes that are being queried. If the kmer is 
    seen in all palindromes, then the kmer is a palindrome.
    The palindrome and its observed value are stored in
    a collections.Counter representedPalindromes. After
    every kmer in kmerCounts has been considered, 
    representedPalindromes is returned. 
    """
    representedPalindromes = collections.Counter()
    # Stores all palindromes that are seen in the 
    # input files, as well as their observed counts. 
    for kmer, value in kmerCounts.items():
        if kmer in allPalindromes:
            representedPalindromes[kmer] = value
    return representedPalindromes 

def print_palindromes_and_stats(resultsAndStatsAndStuff, min_k,
		 max_k, hypothesesCount,cutoff, time):
    """
    INPUT:
	resultsAndStatsAndStuff- A tuple with three elements
	min_k - An integer
	max_k - An integer
	hypothesesCount - An integer
	cutoff - A float
	time - An integer
    OUTPUT:
	sys.stdout - Outputs to terminal.

	Prints a table to standard out. The table has five 
	header lines, total character number, number of hypotheses,
	min and max palindromes considered, the E-value cutoff
	and files used to generate the table.  The table 
	has five columns, kmer, observed, expected, Z-score
	and E-value. The tuples containing the data are printed
        to the table. The output is sorted by Z-score, and only
	tuples with E-values below the cutoff are printed.
        There is a tail line, the time the program took to run
	excluding printing.   
    """
    resultsAndStats = resultsAndStatsAndStuff[0]
    # The list of tuples
    inputFileName = resultsAndStatsAndStuff[1]
    # String containing the file names
    totalCharCount = resultsAndStatsAndStuff[2]
    # Total character count across all files
    resultsAndStats = sorted(resultsAndStats, key = 
		operator.itemgetter(3))
    print ("#Reading from " + inputFileName)
    print ("#Processing " + str(totalCharCount) + " characters")
    print ("#Considering " + str(hypothesesCount) + " words (" 
		+ str(2*hypothesesCount) + " hypotheses)")
    print ("#Printing palindromes with E-value below "
		 + str(cutoff))
    print ("#kmer\t observed\t expected \t Z_score \t E_value")
    for key,value,expectedValue,zScore,eValue in resultsAndStats:
        if eValue < cutoff:
            print(key +"\t\t"+ value + "\t"+ expectedValue
                    +"\t" + str(zScore) + "\t" + str(eValue))
    print ("#Program took " + str(time) + " seconds to run" )

def get_all_palindromes_and_hypotheses_count(min_k,max_k):
    """
    INPUT:
	min_k - An integer
	max_k - An integer
    OUTPUT:
	allPalindromes - A collections.Counter object
	hypothesesCount - An integer
     This function generates all palindromes within the range
     of min_k and max_k.  The total count of these palindromes
     becomes the hypotheses Count.  The function returns a 
     tuple with two elements. The first is a 
     collections.Counter() object holding the palindromes, 
     the second is the integer hypotheses count.  
    """
    allPalindromes = collections.Counter()
    # Holds all palindromes of size min_k to max_k
    hypothesesCount = 0
    # The total number of palindromes being considered
    for i in range (max_k - (min_k -1)):
        if ((min_k+i)%2) == 1:
            hypothesesCount += (4**((min_k+i)/2))*2
        else:
            hypothesesCount += (4**((min_k+i)/2))
        palindromes = get_palindromes( "ACTG",
		min_k +i)   
        allPalindromes += palindromes
    return allPalindromes, hypothesesCount
    

def get_palindromes(alpha, k):
    """
    INPUT:
	alpha- A string
	k- An integer
    OUTPUT:
	palindromes- A collections.Counter object
    Takes in an empty list, the alphabet being used, and the 
    size of the palindromes to be produced. For even 
    palindromes, every possible k/2 mer is produced
    then combined with its corresponding reverse compliment. 
    For odd kmers every possible (k-1)/2 mer is produced.
    The reverse compliment is determined. Two bases are 
    considered as middle bases, this prevents double counting
    odd palindromes (AATCC and AAACC are internally identical).
    The mer is prefixed and the reverse compliment is appended.
    As the palindromes are determined they are added to a list
    after all palindromes have been considered the list is 
    returned 
    """
    palindromes = collections.Counter()
    # This is the empty list of palindromes. The palindromes
    # will be stored as strings in this list. 
    halfSize = 0
    # This is an integer, it will be used to define the 
    # halfway length of a kmer. 
    
    if (k % 2) == 1:
        halfSize = (k-1)/2
    else: halfSize = k/2
 
    kmers = itertools.product(alpha, repeat = halfSize)
    # Kmers is a list of tuples, each tuple is a possible kmer.
    # Every possible kmer is represented in this list of 
    # tuples.  For example ('A', 'B', 'C') corresponds to the
    # string 'ABC'.  The tuples need to be processed into 
    # strings before the reverse compliment is determined.
    
    for item in kmers:
    # Iterate over every element in kmers, when this loop
    # finishes every palindrome of size k will be represented
    # in palindromes. 
        kmerTuple = str(item).translate(None,"(),'")
        # Turn the tuple into a string and remove punctuation
        kmer = ""
        for char in kmerTuple.split():
            kmer += char
        compKmer = reverse_comp(kmer)
        if (k % 2) == 1:
        # If the palindrome is odd then consider internal
        # characters.  This list contains all possible 
        # palindromes, therefore all possible palindromes
        # are represented. These are the palindromes where
        # the middle character is one of the chars, SWATCG
        # where S stands for either G or C and W stands for
        # A or T.  
            for char in "SWATCG":
                palindromes[kmer+ char +compKmer] = 1
        else:
        # Even palindromes are generated by merging kmer and 
        # compKmer. 
            palindromes[kmer+compKmer] = 1
    return palindromes


def count_kmers_and_totalChars(inputFile,  min_k, max_k,  
		alpha, verbose): 
    """ 
    INPUT:
	inputFile - A file like object
	min_k - An integer
	max_k - An integer
	alpha - A string
	verbose - Boolean
    OUTPUT:
	kmerCounts - A collections.Counter() object
	totalCharCounts - An integer
    Takes in an opened input file, a collections.Counter, and
     a integer k.  K corresponds to the order of the markov 
    chain, as well as the kmer size.  The program counts the 
    occurences of each kmer of the specified size. The 
    collections.Counter is returned after the fasta file has 
    been fully  processed. The values in the Counter can be 
    accessed as follows; given a kmer "TTA", result["TTA"] 
    returns the count of "TTA". 
    """
    totalCharCount = 0
    # The total number of characters seen in the file.
    kmerCounts = collections.Counter()
    if alpha is None:
        alpha = string.ascii_uppercase
    for fasta in fastFunctions.readFasta(inputFile,alpha,
			verbose,True):
        if len(fasta.sequence)== 0: continue
        totalCharCount += len(fasta.sequence)
        # Throw away empty sequences.
        for i in range(max_k - (min_k-1)):
            # Iterate over all kmer sizes from min to max
            for start in range(len(fasta.sequence)-(min_k+i)):
                # Iterate over all kmers, incrementing the 
	        # count each time  a kmer is seen. 
                kmer = fasta.sequence[start:start+(min_k+i)]
	        if (len(kmer)%2) == 1:
		   # Handle odd palindromes individually.
		   # Merge A and T center characters into W
		   #  and CG center characters into S. 
		   kmerList = list(kmer)
                   if ( kmer[:(len(kmer)/2)]== reverse_comp(
				kmer[((len(kmer)+1)/2):])):
                       kmerCounts["".join(kmerList)] += 1
                       if kmerList[(len(kmer)/2)] in "AT":
                           kmerList[(len(kmer)/2)] = "W"
                           kmerCounts["".join(kmerList)] += 1
                       if kmerList[(len(kmer)/2)] in "CG":
                           kmerList[(len(kmer)/2)] = "S"
                           kmerCounts["".join(kmerList)] += 1
                   else: kmerCounts["".join(kmerList)] += 1
                else: kmerCounts[fasta.sequence[start:start
				+(min_k + i)]] += 1 
                # fasta.sequence[start:start+k] is a single 
		# kmer string. These kmers are the key in
		# counts, so counts["kmer"] returns the count 
		# of "kmer". 
    return kmerCounts, totalCharCount
