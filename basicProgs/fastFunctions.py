#!/usr/bin/env python2.7
# fastFunctions.py REDO
# Chris Eisenhart 12.04.2014
# ceisenhart@soe.ucsc.edu/ceisenha@ucsc.edu
from __future__ import print_function
import sys, string, operator, fileinput, collections, os.path
import  re, argparse
from itertools import izip
"""
A module for fasta and fastq file I/O.
The functions contained in this module are;

warning (string) - Prints a warning message

error (string)- Prints a warning message and closes the program

printFastqSeq (opened file, fastq, int) - Prints a fastq 
	sequence to the specified file

printFastaSeq (opened file, fastq)- Prints a fasta sequence to 
	the specified file

printQualSeq (opened file, fastq)- Prints a quality sequence to
	 the specified file

readQuality (opened file, bool) - Yields a quality sequence at 
	a time, iterates over a full file.

readFasta (opened file, string, bool, bool) - Yields a fasta 
	sequence at a time, iterates over a full file.
 
readFastq (opened file, string, bool, bool, int) - Yields a 
	fastq sequence at a time, iterates over a full file. 

readFastaWithQuality (opened file, opened file, string, bool, 
				bool)-
	Yields a fastq sequence at a time, iterates over a
	 full file. 
"""

class fq:
    header = None
    # Stores the header line as a a string.
    extraData = "" 
    # Any extra data on the header line is stored as a string
    sequence = ""
    # The DNA sequence stored internally as a string.
    quality = []
    # The quality sequence will be stored as a list of integers. 
    

def warning(*objs):
    """
    INPUT:
	objs - A string. 
    OUTPUT:
	None
    Prints the input string to sys.stderr.
    """
    print("WARNING: ", *objs, file=sys.stderr)

def error(*objs):
    """
    INPUT:
	objs - A string. 
    OUTPUT:
	None
    Prints the input string to sys.stderr. Then exit the
    program. 
    """
    print("ERROR: ", *objs, file=sys.stderr)
    sys.exit()

def printFastqSeq(output, fastqSeq, phred):
    """
    INPUT:
	output - A file like object.
	fastqSeq - A fastq object.
	phred - An integer. 
    OUTPUT:
	None
    Prints a single fastq sequence to output.
    The quality scores are offset by the integer phred. 
    """
    print ("@" + fastqSeq.header +  fastqSeq.extraData
		, file = output )
    dnaLine = ""
    for word in fastqSeq.sequence:
        dnaLine += word
    print (dnaLine, file = output)
    print ("+", file = output)
    quality = ""
    for score in fastqSeq.quality:
        quality += str(chr(score + phred))
    print (quality, file = output)

def printFastaSeq(output, fastaSeq, limitLineLength):
    """
    INPUT:
	output - A file like object. 
	fastaSeq- A fastq object. 
    OUTPUT:
	None
    Prints a single fasta sequence to output. The fasta 
    sequence is stored in fastaSeq, with empty quality 
    variables. 
    """
    print (">" + fastaSeq.header + " " + fastaSeq.extraData,
		 file = output)
    dnaLine = ""
    for word in fastaSeq.sequence:
        dnaLine += word
	if limitLineLength:
	    if len(dnaLine) == 80:
		print (dnaLine, file = output)
		dnaLine = ""
    print (dnaLine, file = output)

def printQualSeq(output, qual):
    """
    INPUT:
	output - A file like object.
	qual - A fastq object.
    OUTPUT:
	None
    Prints a single quality sequence to output. The quality 
    sequence is stored in qual, with empty sequence variables.
    """
    print (">" + qual.header + " " +qual.extraData,
		 file = output)
    quality = ""
    for score in qual.quality:
        quality += str(score)
        quality += "   "
    print (quality, file = output)

def readQuality(qualFile, verbose):
    """
    INPUT:
	qualFile - A file like object.
	verbose - A boolean. 
    OUTPUT:
	qual - A fastq object. 
    Reads through the quality file, one sequence at a time.
    Each sequence is stored in a fastq object which are yielded
    as they are read. 
    """
    qual = fq()
    # Stores the header and quality scores, will be yielded for
    # each quality sequence.  
    for line in qualFile:
        if line.startswith('>'):
	# This is a header line, yield completed sequences then
	# reset the fields of the fastq class. 
            if qual.header is None:
                qual.header = line[1:]
                continue
            else: 
                yield qual
                qual.header = line[1:]
                qual.quality = []
                continue
        else:    
	# This is a quality line, read the scores in and 
	# append them to the quality score list. 
	    qual.quality+=[int(score) for score in line.split()]
	    continue
    if qual.header is not None: yield qual 

def readFasta(inputFile, alphabet, verbose, multiCases): 
    """
    INPUT:
	inputFile - A file like object. 
	alphabet - A string.
	verbose - A boolean.
	multiCases - A boolean. 
    OUTPUT:
	fasta - A fastq class.
    Takes in a file like object and yield fasta sequences
    from the file.  The sequences are yielded in a fasta
    class one sequence at a time. If characters are found
    that are not in alphabet then a warning is printed. 
    If multicases is true, then the program prints characters
    in both capital and lower case. 
    """
    fasta = fq()
    # This holds all the DNA sequences, it is yielded and reset
    # for each fasta sequence. 
    for line in inputFile:
	# This handles the first fasta sequence header
        if line.startswith('>') and fasta.header is None:
            splitLine = re.search("[\s,]", line) 
            fasta.header = line[1:splitLine.start(0)]
            fasta.extraData = line[splitLine.start(0):-1]
            continue
	# This handles all fasta sequence headers after the 
	# first sequence. Additionally the previous fasta is
	# yielded and reset.
	if line.startswith('>') and fasta.header is not None: 
	    yield fasta
            splitLine = re.search("[\s,]", line) 
            fasta.header = line[1:splitLine.start(0)]
            fasta.extraData = line[splitLine.start(0):-1]
	    fasta.sequence = "" 
	    continue
	# This handles the DNA sequence, if the character
	# is in the alphabet then the capitalized version is 
	# added to the DNA sequence.  If multicases is true
	# then the character is not capitalized. 
        for char in line:
	    if char in alphabet:
	        if multiCases: fasta.sequence += char
	        else: fasta.sequence += char.upper()
		         
    if fasta.header is not None: yield fasta        


def readFastq(inputFile, alphabet, verbose, multiCases, phred): 
    """
    INPUT:
	inputFile - A file like object. 
	alphabet - A string.
	verbose - A boolean.
	multiCases - A boolean. 
	phred - An integer.
    OUTPUT:
	fastq - A fastq class.
    Takes in a file like object and yield fasta sequences
    from the file.  The sequences are yielded in a fastq
    class one sequence at a time. If characters are found
    that are not in alphabet then an error is printed. 
    If multicases is true, then the program prints characters
    in both capital and lower case. Phred corresponds to an
    ofset applies to quality scores when they are printed.   
    """
    fastq = fq()
    afterDel = False
    # This boolean keeps track of whether the fastq
    # is before or after the deliminator. 
    completedSequence = False
    # This boolean keeps track of whether the quality 
    # sequence length equals the dna sequence length, 
    # which marks the end of the fastq read. 
    for line in inputFile:
	# Handle the first fastq reads header
	if fastq.header is None:
            splitLine = re.search("[\s,]", line) 
            fastq.header = line[1:splitLine.start(0)]
            fastq.extraData = line[splitLine.start(0):-1]
            continue

	# Yields the completed sequence and handles the 
	# fastq header.  
	if (line.startswith("@") and len(fastq.sequence) == 
		    len(fastq.quality) and completedSequence):
            yield fastq
	    fastq.sequence = ""
	    fastq.quality = [] 
	    afterDel = False
	    completedSeqeunce = False
            splitLine = re.search("[\s,]", line) 
            fastq.header = line[1:splitLine.start(0)]
            fastq.extraData = line[splitLine.start(0):-1]
            continue

	#This block handles the dna sequence and the 
	# deliminator. 
        if not afterDel and fastq.header is not None:
	    if line.startswith("+"):
		# This breaks the DNA sequence growth
		# and signals that the next lines
		# will be quality sequence
	        afterDel = True
	        continue
	    for char in line:
	        if char in alphabet:
	            if multiCases: fastq.sequence += char
	            else: fastq.sequence += char.upper()
	    continue

	# This block handles the quality sequence
	if afterDel:
	    for score in line[:-1]:
		fastq.quality.append(ord(score) - phred)
	    if len(fastq.sequence) == len(fastq.quality):
		# When the quality sequence length is the 
		# same as the DNA sequence length then 
		# the fastq object is completed.  
		completedSequence = True
		afterDel = False

    if fastq.header is not None: yield fastq

def readFastaWithQuality(fastaFile, qualFile, alphabet,
			 verbose, multiCases):
    """
    INPUT:
	fastaFile - A file like object. 
	qualFile - A file like object. 
	alphabet - A string.
	verbose - A boolean.
	multiCases - A boolean. 
    OUTPUT:
	fastq - A fastq class.
    Takes in two file like objects and yields fastq sequences
    created.  The sequences are yielded in a fastq
    class one sequence at a time. If characters are found
    that are not in alphabet then a warning is printed. 
    If multicases is true, then the program prints characters
    in both capital and lower case. 
    """
    fastq = fq()
    # This will merge the fasta and qual reads into a single
    # fastq object. 
    for fastaLine, qualLine in izip(readFasta(fastaFile,
			 alphabet, verbose, multiCases),
			readQuality(qualFile, verbose)):
        fastq.header = fastaLine.header
        fastq.extraData = fastaLine.extraData
	fastq.sequence = fastaLine.sequence
        fastq.quality = qualLine.quality
        if len(fastq.sequence) != len(fastq.quality):
	    error("Sequence length and quality length differ."
		" The input files may have been corrrupted. ")		
	yield fastq
