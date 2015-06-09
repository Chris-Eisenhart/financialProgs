#!/usr/bin/env python2.7SD
# align.py
# Chris Eisenhart  ceisenha@ucsc.edu/ceisenhart@soe.ucsc.edu
# 11.20.2014-12.10.2014 
"""
This module contains functions used for sequence alignment, 
the sequences can be either protein or nucleotide. The user
must provide a score table. The primary function is alignFile, 
 see program test_align for an example. 

This module contains the following classes;

local_aligner(collections.Counter, int, int, int):

global_aligner(collections.Counter, int, int, int):

This module contains the following functions;

getScoreMatrix (string):

makeMatrix(collections.Counter, string, string, int, int, int, bool):

getMScore(tuple, int, bool):

getCScore(tuple, int, int, int):

getRScore(tuple, int, int, int):

setGlobalInitialConditions(dict, string, string, int, int):

getLocalScore(int, int, dict, collections.Counter, string, string, int, int, int):

getGlobalScore(int, int, dict, collections.Counter, string, string, int, int, int):

makeMatrix(collections.Counter, string, string, int, int, int, bool):

findHighScoreLocal(dict, string, string):

globalTraceBack(dict, collections.Counter, string, string, int, int, int): 

localTraceBack(dict, collections.Counter, string, string, int, int, int, tuple): 

scoreA2MLocal(string, string, collections.Counter, int, int, int):

scoreA2MGlobal(string, string, collections.Counter, int, int, int):

alignFile(file, int, int, int, string, collections.Counter, bool, bool):
"""

from __future__ import print_function  
import sys, operator, fileinput, collections, urllib
import string, os.path, re, argparse, bisect
import gzip, math, random, itertools, fastaFastqParser


class local_aligner:
    """
    This class is used for local alignments. The class has 
    four functions;
    init (collections,counter, int, int, int) - Initializes the
        class object, each object has a substitution matrix 
	as well as open, extend and double gap values. 
    score_a2m(string, string) Calculates and prints the A2M 
        local score associated with the two strings. 
    score_align(string, string) Generates a local alignment
        matrix for the two strings. The highest score is 
	located and returned. 
    traceback() Starts at the highest scoring location and
        traces backwards through the matrix. The alignment
	generated is returned in A2M format. 
    """
    def __init__(self, subst, openCost=-12, extendCost=-1, doubleGapCost=-3):
        """
	INPUT: 
   	    subst - A collections.Counter object, it keeps
	        substitution matrix scores for amino acid pairs. 
	        The amino acid pairs are represented as a string, 
	        the score for each pair is stored as an integer. 
	    openCost - An integer, the cost penalty for opening a
	        new gap.
	    extendCost - An integer, the cost penalty for extending
	        a gap.
	    doubleCost - An integer, the cost penalty for a double
	        gap.
	OUTPUT:
	    none
	Initializes values for this instance of the class. 
	The acceptable alphabet is identified from the 
	substitution matrix. 
	"""
	self.subst = subst
	self.alphabet = ""
	for key, value in subst.iteritems():
	    self.alphabet += key
	self.openCost = openCost
	self.extendCost = extendCost
	self.doubleCost = doubleGapCost
    
    def score_a2m(self,s1,s2):
	"""
	INPUT:
	    s1 - A string, this represents the master sequence.
	    s2 - A string, this represents the slave sequence. 
	OUTPUT:
	    none
	"""
	print (scoreA2MLocal (s1, s2, self.subst, self.doubleCost, self.extendCost, self.openCost))

    def score_align(self, s1, s2):
	"""
	INPUT:
	    s1 - A string, this represents the master sequence.
	    s2 - A string, this represents the slave sequence. 
	OUTPUT:
	    self.highScore[0] - An integer, the highest value 
		seen in the alignment matrix. 
	Takes in two strings and generates an alignment matrix
	for the two strings. The highest value in the matrix is
	identified and returned. 
	"""
	self.youngSlave = s2.replace(" ", "").upper()
        self.youngMaster = s1.replace(" ", "").upper()
	self.slave = ""
	self.master = ""
	# Remove characters that are not elements of a
	# substitution matrix index. 
	for char in self.youngSlave:
	    if char in self.alphabet:
		self.slave += char
	for char in self.youngMaster:
	    if char in self.alphabet:
		self.master += char
	self.matrix = makeMatrix(self.subst, self.master,
		 self.slave, self.openCost, self.extendCost,
		 self.doubleCost, True)
	self.highScore = findHighScoreLocal(self.matrix, self.master, self.slave)
	return (self.highScore[0])

    def traceback(self):	
	"""
	INPUT:
	    none
	OUTPUT:
	    self.result - A string, this string is the optimal
		local alignment in A2M format. 
	Takes the values currently stored in the alignment 
	object and performs a traceback. The best local 
	alignment is determined and returned in A2M format. 
	"""
	self.result = localTraceBack(self.matrix, self.subst, 
		self.master, self.slave, self.openCost,
		 self.extendCost, self.doubleCost, self.highScore) 
	return self.result
	
class global_aligner:
    """
    This class is used for global alignments. The class has 
    four functions;
    init (collections,counter, int, int, int) - Initializes the
        class object, each object has a substitution matrix 
	as well as open, extend and double gap values. 
    score_a2m(string, string) Calculates and prints the A2M 
        global score associated with the two strings. 
    score_align(string, string) Generates a global alignment
        matrix for the two strings. The highest score is 
	located and returned. 
    traceback() Starts at the highest scoring location and
        traces backwards through the matrix. The alignment
	generated is returned in A2M format. 
    """
    def __init__(self, subst, openCost=-12, extendCost=-1, doubleGapCost=-3):
        """
	INPUT: 
   	    subst - A collections.Counter object, it keeps
	        substitution matrix scores for amino acid pairs. 
	        The amino acid pairs are represented as a string, 
	        the score for each pair is stored as an integer. 
	    openCost - An integer, the cost penalty for opening a
	        new gap.
	    extendCost - An integer, the cost penalty for extending
	        a gap.
	    doubleCost - An integer, the cost penalty for a double
	        gap.
	OUTPUT:
	    none
	Initializes values for this instance of the class. 
	The acceptable alphabet is identified from the 
	substitution matrix. 
	"""
	self.subst = subst
	self.alphabet = ""
	for key, value in subst.iteritems():
	    self.alphabet += key
	self.openCost = openCost
	self.extendCost = extendCost
	self.doubleCost = doubleGapCost
    
    def score_a2m(self,s1,s2):
	"""
	INPUT:
	    s1 - A string, this represents the master sequence.
	    s2 - A string, this represents the slave sequence. 
	OUTPUT:
	    none
	"""
	print (scoreA2MGlobal (s1, s2, self.subst, self.doubleCost, self.extendCost, self.openCost))

    def score_align(self, s1, s2):
	"""
	INPUT:
	    s1 - A string, this represents the master sequence.
	    s2 - A string, this represents the slave sequence. 
	OUTPUT:
	    self.highScore - An integer, the highest last value 
		seen in the alignment matrix. 
	Takes in two strings and generates an alignment matrix
	for the two strings. The highest value in the matrix is
	identified and returned. 
	"""
	self.youngSlave = s2.replace(" ", "").upper()
        self.youngMaster = s1.replace(" ", "").upper()
	self.slave = ""
	self.master = ""
	# Remove characters that are not elements of a
	# substitution matrix index. 
	for char in self.youngSlave:
	    if char in self.alphabet:
		self.slave += char
	for char in self.youngMaster:
	    if char in self.alphabet:
		self.master += char
	self.matrix = makeMatrix(self.subst, self.master,
		 self.slave, self.openCost, self.extendCost,
		 self.doubleCost, False)
	self.highScore = self.matrix[(len(self.slave)-1, len(self.master)-1)]
	return max(self.highScore)

    def traceback(self):	
	"""
	INPUT:
	    none
	OUTPUT:
	    self.result - A string, this string is the optimal
		global alignment in A2M format. 
	Takes the values currently stored in the alignment 
	object and performs a traceback. The best global
	alignment is determined and returned in A2M format. 
	"""
	self.result = globalTraceBack(self.matrix, self.subst, 
		self.master, self.slave, self.openCost,
		 self.extendCost, self.doubleCost) 
	return self.result

def argmax(l):
    """
    INPUT:
	l - A list or iterable object. For this program
	    l will be a tuple with either three or four
	    elements. 
    OUTPUT:
	l.index(max(l)) - An integer, this is the location of  
	    the maximum element.
	max(l) - An integer, this is the maximum elements. 
    """
    return l.index(max(l)), max(l)

def getScoreMatrix(URL):
    """
    INPUT:
	URL- A string. 
    OUTPUT:	
	scores - A collections.Counter object. 
    Reads in a score matrix from the input URL. Use "None" as 
    the input URL to default to BLOSUM62. The score matrix is 
    returned as a collections.Counter object, where the key is
    two amino acids concattenated, and the value is the score 
    associated with the amino acids. 
    """
    alphabet = "ARNDCQEGHILKMFPSTWYVBZX*"
    page = urllib.urlopen(URL)
    # scoreMatrix will store all the substitution scores, two 
    # amino acids will have their score stored in scoreMatrix
    # as follows; scoreMatrix[AA] = 4
    scoreMatrix = collections.Counter()
    for line in page:
	# Skip comment lines. 
        if line.startswith("#"):continue
	# The first non comment line is a row of amino acids.
	# It is split and stored in splitLine. 
        splitLine = line.split()
	if splitLine[1] in alphabet: 
	    aminoAcids = splitLine
	    continue
	# This evaluates all rows after the first 
	firstAA = splitLine[0]
	charSpot = 0
	for char in splitLine[:-1]:
	    if char in aminoAcids: 
		continue
	    key = aminoAcids[charSpot] + firstAA
	    scoreMatrix[key] = int(char)
	    charSpot += 1
    return scoreMatrix

def getMScore(tupleCell, score, isLocal):
    """
    INPUT:
	tupleCell- The elements are
	    (rowScore, colScore, diagonalScore).
	    RowScore is the score obtained from the cell
	    in the x-1 position, colScore is the score obtained from 
	    the y-1 position, and diagonalScore is the score
	    obtained from the cell in the x-1 y-1 position. 
	score - An integer, the substitution matrix score
	    associated with the current position. 
	isLocal - A boolean, it is used to differentiate
	    between local and global alignments
    OUTPUT:
	result - An integer, this integer is combined with
	    two others to create a new tupleCell. This 
	    specific integer corresponds to the diagonalScore. 
    """
    tupleList = list(tupleCell)
    if isLocal:tupleList.append(0)
    result = max(tupleList) + score
    return result


def getCScore(tupleCell, doubleCost, extendCost, openCost):
    """
    INPUT:
	tupleCell- The elements are
	    (rowScore, colScore, diagonalScore).
	    RowScore is the score obtained from the cell
	    in the x-1 position, colScore is the score obtained from 
	    the y-1 position, and diagonalScore is the score
	    obtained from the cell in the x-1 y-1 position. 
	doubleCost - An integer, the cost penalty for a double
	    gap.
	extendCost - An integer, the cost penalty for extending
	    a gap.
	openCost - An integer, the cost penalty for opening a
	    new gap.
    OUTPUT:
	result - An integer, this integer is combined with
	    two others to create a new tupleCell. This 
	    specific integer corresponds to the colScore. 
    """
    tupleList = list(tupleCell)
    tupleList[0] += doubleCost
    tupleList[1] += extendCost
    tupleList[2] += openCost
    return max(tupleList)

def getRScore(tupleCell, doubleCost, extendCost, openCost):
    """
    INPUT:
	tupleCell- The elements are
	    (rowScore, colScore, diagonalScore).
	    RowScore is the score obtained from the cell
	    in the x-1 position, colScore is the score obtained from 
	    the y-1 position, and diagonalScore is the score
	    obtained from the cell in the x-1 y-1 position. 
	doubleCost - An integer, the cost penalty for a double
	    gap.
	extendCost - An integer, the cost penalty for extending
	    a gap.
	openCost - An integer, the cost penalty for opening a
	    new gap.
    OUTPUT:
	result - An integer, this integer is combined with
	    two others to create a new tupleCell. This 
	    specific integer corresponds to the rowScore. 
    """
    tupleList = list(tupleCell)
    tupleList[0] += extendCost
    tupleList[1] += doubleCost
    tupleList[2] += openCost
    return max(tupleList)

def setGlobalInitialConditions(matrix, master, slave, openCost, extendCost):
    """
    INPUT:
	matrix - A dict object, this stores the alignment matrix
	   for two sequences. 
	master - A string, this stores the master sequence.
	slave - A string, this stores the slave sequence. 
	openCost - An integer, the cost penalty for opening a
	    new gap.
	extendCost - An integer, the cost penalty for extending
	    a gap.
    OUTPUT:
	none
    """
    masterSpot = 0
    slaveSpot = 0
    matrix [(-1,-1)] = (-1000, -1000, 0)
    # Set the initial conditions for the first row
    for char in master:
        matrix[(-1,masterSpot)]=(openCost + masterSpot*(extendCost), -1000, -1000)
	masterSpot += 1
    # Set the initial conditions for the first col
    for char in slave:
        matrix[(slaveSpot,-1)]=(-1000, openCost +slaveSpot*(extendCost), -1000)
	slaveSpot += 1

def getLocalScore(row, col, matrix, subst, master, slave, openCost, extendCost, doubleCost):
    """
    INPUT:
	row - An integer, keeps track of the slave location. 
	col - An integer, keeps track of the master location. 
	matrix - A dict object, this stores the alignment matrix
	   for two sequences. 
	subst - A collections.Counter object, it keeps
	    substitution matrix scores for amino acid pairs. 
	    The amino acid pairs are represented as a string, 
	    the score for each pair is stored as an integer. 
	master - A string, this stores the master sequence.
	slave - A string, this stores the slave sequence. 
	openCost - An integer, the cost penalty for opening a
	    new gap.
	extendCost - An integer, the cost penalty for extending
	    a gap.
	doubleCost - An integer, the cost penalty for a double
	    gap.
    OUTPUT:
	result - An integer
    """
    master = list(master)
    slave = list(slave)
    score = subst[master[col]+slave[row]]
    if col is 0 and row is 0:
	return (-10000,-10000,score)
    # Handle the boundary conditions for the first row
    if row is 0 and col > 0:
	return (-10000, -10000, score)
    # Handle the boundary conditions for the first row
    if col is 0 and row > 0:
	return (-10000, -10000, score)
    # Calculate the row, column and diagonal scores. 
    rScore = getRScore(matrix[(row,col-1)], doubleCost, extendCost, openCost)
    cScore = getCScore(matrix[(row-1,col)], doubleCost, extendCost, openCost)
    mScore = getMScore(matrix[(row-1,col-1)], score, True)
    return (rScore, cScore, mScore)

def getGlobalScore(row, col, matrix, subst, master, slave, openCost, extendCost, doubleCost):
    """
    INPUT:
	row - An integer, keeps track of the slave location. 
	col - An integer, keeps track of the master location. 
	matrix - A dict object, this stores the alignment matrix
	   for two sequences. 
	subst - A collections.Counter object, it keeps
	    substitution matrix scores for amino acid pairs. 
	    The amino acid pairs are represented as a string, 
	    the score for each pair is stored as an integer. 
	master - A string, this stores the master sequence.
	slave - A string, this stores the slave sequence. 
	openCost - An integer, the cost penalty for opening a
	    new gap.
	extendCost - An integer, the cost penalty for extending
	    a gap.
	doubleCost - An integer, the cost penalty for a double
	    gap.
    OUTPUT:
	result - An integer
    """
    master = list(master)
    slave = list(slave)
    score = subst[master[col]+slave[row]]
    # Calculate the row, column and diagonal scores. 
    rScore = getRScore(matrix[(row,col-1)], doubleCost, extendCost, openCost)
    cScore = getCScore(matrix[(row-1,col)], doubleCost, extendCost, openCost)
    mScore = getMScore(matrix[(row-1,col-1)], score, False)
    return (rScore, cScore, mScore)

def makeMatrix(subst, master, slave, openCost, extendCost, doubleGapCost, isLocal):
    """
    INPUT:
	subst - A collections.Counter object, it keeps
	    substitution matrix scores for amino acid pairs. 
	    The amino acid pairs are represented as a string, 
	    the score for each pair is stored as an integer. 
	master - A string, this stores the master sequence.
	slave - A string, this stores the slave sequence. 
	openCost - An integer, the cost penalty for opening a
	    new gap.
	extendCost - An integer, the cost penalty for extending
	    a gap.
	doubleCost - An integer, the cost penalty for a double
	    gap.
	isLocal - A boolean, it is used to differentiate
	    between local and global alignments
    OUTPUT:
	matrix - A list of lists. 
    Takes in the score matrix and the strings to be scored.
    Before processing the strings a prefix and suffix character
    is added. This allows for the first column and row to be
    defaulted to -10000, which accurately models the biological
    alignment problem of starting gaps.  
    Each row is iterated over, for each cell the score is
    determined.
    """
    # matrix will store the alignment score matrix for the
    # two input strings. The matrix is indexed using tuples
    # with two elements, for example matrix[(-1,-1)] 
    # corresponds to the cell at x = -1 and y = -1. Each cell
    # is stored as a tuple with three elements, the first
    # is the row score, the second is the col score, the third
    # is the diagonal score. Therefore the matrix element
    # at matrix[(-1,-1)] could be (-1000, -1000, 0) (for a 
    # global alignment). 
    matrix = dict()
    if not isLocal:
	setGlobalInitialConditions(matrix, master, slave, openCost, extendCost)
    for row in range(len(slave)):
	for col in range(len(master)):
	    if isLocal: 
		matrix[(row,col)] = getLocalScore(row, col, matrix, subst, master, slave, openCost, extendCost, doubleGapCost)
	    else:
		matrix[(row,col)] = getGlobalScore(row, col, matrix, subst, master, slave, openCost, extendCost, doubleGapCost)
    return matrix
    
def findHighScoreLocal(matrix, master, slave):
    """
    INPUT:
	matrix - A dict object, this stores the alignment matrix
	   for two sequences. 
	master - A string, this stores the master sequence.
	slave - A string, this stores the slave sequence. 
    OUTPUT:
	maxScore - An integer, the highest score seen in the 
	    alignment matrix
	rowLocation - An integer, the location in the slave
	    that the highest score was seen at 
	colLocation - An integer, the location in the master
	    that the highest score was seen at 
    """
    maxScore = 0
    rowLocation = -1
    colLocation = -1 
    for key, value in matrix.iteritems():
	if max(value) > maxScore:
	    maxScore = max(value)
	    rowLocation = key[0]
	    colLocation = key[1]
    return maxScore, rowLocation, colLocation


def globalTraceBack(matrix, subst, master, slave, openCost, extendCost, doubleCost):
    """
    INPUT:
	matrix - A dict object, this stores the alignment matrix
	   for two sequences. 
	subst - A collections.Counter object, it keeps
	    substitution matrix scores for amino acid pairs. 
	    The amino acid pairs are represented as a string, 
	    the score for each pair is stored as an integer. 
	master - A string, this stores the master sequence.
	slave - A string, this stores the slave sequence. 
	openCost - An integer, the cost penalty for opening a
	    new gap.
	extendCost - An integer, the cost penalty for extending
	    a gap.
	doubleCost - An integer, the cost penalty for a double
	    gap.
    OUTPUT:
	result- A string, this string is the alignment generated
	    by tracing back through the matrix. 
    """ 
    # Turn the string into a list so that it can be directly
    # indexed
    master = list(master)
    slave = list(slave)
    # The result will store the alignment, characters will be 
    # appended as they are encountered.
    result = ""
    # For global alignments the traceback starts at the last
    # cell. currentSlaveLoc corresponds to the y axis, and
    # currentMasterLoc corresponds to the x axis. 
    currentSlaveLoc = len(slave) - 1
    currentMasterLoc = len(master) - 1
    maxes = argmax(matrix[(currentSlaveLoc,currentMasterLoc)])
    # nextLoc keeps track of the next cell in the traceback. 
    # Each state is given an integer value, the x-1 state 
    # is assigned 0, the y-1 state is 1, and the x-1, y-1 state
    # is assigned 2. 
    nextLoc = maxes[0]
    # Traceback untill the alignment runs off the matrix. 
    while currentMasterLoc >= 0 and currentSlaveLoc >= 0:
	    if nextLoc == 0:
		# This is a left jump
		values = list(matrix[(currentSlaveLoc,currentMasterLoc - 1)])
		values[0] += extendCost
		values[1] += doubleCost
		values[2] += openCost
		nextMax = argmax(values)
		nextLoc = nextMax[0]
		currentMasterLoc -= 1
		result += "-"
		continue
	    if nextLoc == 1:
		# This is an upward jump
		values = list(matrix[(currentSlaveLoc-1,currentMasterLoc)])
		values[0] += doubleCost
		values[1] += extendCost
		values[2] += openCost
		nextMax = argmax(values)
		nextLoc = nextMax[0]
		result += slave[currentSlaveLoc].lower()
		currentSlaveLoc -= 1
		continue
	    if nextLoc == 2:
		# This is a diagonal jump
		values = list(matrix[(currentSlaveLoc - 1, currentMasterLoc - 1)])
		score = subst[master[currentMasterLoc]+slave[currentSlaveLoc]]
		nextMax = argmax(values)
		values[2] += score
		nextLoc = nextMax[0]
		result += slave[currentSlaveLoc].upper()
		currentMasterLoc -= 1
		currentSlaveLoc -= 1
		continue
    # This is the case where the traceback runs off the
    # top edge
    if currentMasterLoc >= 0:
    	result += (currentMasterLoc+1)*"-"
    # This is the case where the traceback runs off the
    # left edge
    if currentSlaveLoc >= 0: 
	result += "".join(slave[0:1+currentSlaveLoc]).lower() 
    result = result[::-1]
    return result


def localTraceBack(matrix, subst, master, slave, openCost, extendCost, doubleCost, highScoreAndLocation): 
    """
    INPUT:
	matrix - A dict object, this stores the alignment matrix
	   for two sequences. 
	subst - A collections.Counter object, it keeps
	    substitution matrix scores for amino acid pairs. 
	    The amino acid pairs are represented as a string, 
	    the score for each pair is stored as an integer. 
	master - A string, this stores the master sequence.
	slave - A string, this stores the slave sequence. 
	openCost - An integer, the cost penalty for opening a
	    new gap.
	extendCost - An integer, the cost penalty for extending
	    a gap.
	doubleCost - An integer, the cost penalty for a double
	    gap.
	highScoreAndLocation - A tuple with three integer 
	    elements, the first is the max score. The second
	    and third identify where in the matrix the max
	    score was found. The second element is the 
	    location in the slave string, the third element
	    is the location in the master string. 
    OUTPUT:
	result- A string, this string is the alignment generated
	    by tracing back through the matrix. 
    """ 
    # Turn the string into a list so that it can be directly
    # indexed
    master = list(master)
    slave = list(slave)
    # The result will store the alignment, characters will be 
    # appended as they are encountered.
    result = ""
    # unpack the highScoreAndLocation tuple
    # highScore is an integer which represents the highest
    # score seen in the matrix
    highScore = highScoreAndLocation[0]
    # currentSlaveLoc in an integer which corresponds
    # to the y location of the high score in the alignment 
    # matrix. 
    currentSlaveLoc = highScoreAndLocation[1] 
    # currentMasterLoc in an integer which corresponds
    # to the x location of the high score in the alignment 
    # matrix. 
    currentMasterLoc = highScoreAndLocation[2]
    # The local alignment may not start at the end of both
    # sequences, suffix handles the characters after
    # the alignment in the query/slave sequence. endGap 
    # handles the characters after the alignment in the
    # master sequence. 
    suffix = "".join(slave[currentSlaveLoc + 1:]).lower()
    endGap = (len(master) - currentMasterLoc  - 1) * "-"
    # nextLoc keeps track of the next cell in the traceback. 
    # Each state is given an integer value, the x-1 state 
    # is assigned 0, the y-1 state is 1, and the x-1, y-1 state
    # is assigned 2. For local alignments the alignment always
    # ends on a match, so it begins at 2.  
    nextLoc = 2
    # The local alignment ends when the alignment runs off the matrix or 
    # the highScore drops to 0. 
    while currentMasterLoc > 0 and currentSlaveLoc > 0 and highScore >= 0:
	    if nextLoc == 0:
		# This is a left jump
		values = list(matrix[(currentSlaveLoc,currentMasterLoc - 1)])
		values[0] += extendCost
		values[1] += doubleCost
		values[2] += openCost
		nextMax = argmax(values)
		nextLoc = nextMax[0]
		if nextLoc == 0:
		    highScore -= extendCost 
		if nextLoc == 1:
		    highScore -= doubleCost
		if nextLoc == 2:
		    highScore -= openCost
		currentMasterLoc -= 1
		result += "-"
		continue
	    if nextLoc == 1:
		# This is an upward jump
		values = list(matrix[(currentSlaveLoc-1,currentMasterLoc)])
		values[0] += doubleCost
		values[1] += extendCost
		values[2] += openCost
		nextMax = argmax(values)
		nextLoc = nextMax[0]
		if nextLoc == 0:
		    highScore -= doubleCost 
		if nextLoc == 1:
		    highScore -= extendCost
		if nextLoc == 2:
		    highScore -= openCost
		result += slave[currentSlaveLoc].lower()
		currentSlaveLoc -= 1
		continue
	    if nextLoc == 2:
		# This is a diagonal jump
		values = list(matrix[(currentSlaveLoc - 1, currentMasterLoc - 1)])
		values.append(0)
		score = subst[master[currentMasterLoc]+slave[currentSlaveLoc]]
		nextMax = argmax(values)
		nextLoc = nextMax[0]
		highScore -= score
		result += slave[currentSlaveLoc].upper()
		currentMasterLoc -= 1
		currentSlaveLoc -= 1
		if nextMax[1] == 0: 
		    break
		continue
    # The local alignment may end one iteration early
    # handle the last iteration. 
    if currentSlaveLoc == 0:
        if nextLoc == 2:
	    result += slave[currentSlaveLoc].upper()
	    currentSlaveLoc -= 1
	    currentMasterLoc -= 1
    result = result[::-1] 
    # The local alignment may not have aligned to the
    # start of both sequences, in these cases the remaining 
    # slave characters are appended as insertions, and 
    # the master characters are appended as gaps. 
    if currentMasterLoc >= 0:
    	result = (currentMasterLoc+1)*"-" + result
    if currentSlaveLoc >= 0: 
	result = "".join(slave[0:1+currentSlaveLoc]).lower() + result
    result = result + endGap + suffix
    return result

def scoreA2MLocal(master, slave, subst, doubleCost, extendCost, openCost):
    """
    INPUT:
	master - A string, this will act as the reference
	slave - A string, this will act as a the query
	subst - A collections.Counter object, it keeps
	    substitution matrix scores for amino acid pairs. 
	    The amino acid pairs are represented as a string, 
	    the score for each pair is stored as an integer. 
	doubleCost - An integer, the cost penalty for a double
	    gap.
	extendCost - An integer, the cost penalty for extending
	    a gap.
	openCost - An integer, the cost penalty for opening a
	    new gap.
    OUTPUT:
	result - An integer, the alignment score associated  
	    with the master slave strings. 
    """
    # The result is initialized to 0, scores and penalties
    # will modify this value. 
    result = 0
    # Turn the string into a list so that it can be directly
    # indexed
    master = list(master.replace(" " , ""))
    slave = list(slave.replace(" " , ""))
    # This boolean will keep track of whether the sequence
    # is in a master gap
    masterGap = False
    # This boolean will keep track of whether the sequence
    # is in a slave insertion
    slaveInsertion = False
    masterCount = 0
    slaveCount = 0
    # For local alignments leading gaps/insertions are ignored
    # this increments the matrix indexes to correspond to the 
    # first aligned character.
    for char in slave:
        if char in string.ascii_lowercase:
	    slaveCount += 1
	if char == "-":
	    slaveCount += 1
	    masterCount += 1
	if char in string.ascii_uppercase: break
    alignEnd = 0
    # For local alignments trailing gaps/insertions are also ignored
    # this is handles by keeping a value to indicate the end of the
    # alignment.
    for char in slave[::-1]:
	alignEnd +=1
	if char in string.ascii_uppercase: break
    # endPoint marks the end of the alignment
    endPoint = len(slave) - alignEnd + 1
    # A boolean used to skip over leading gaps and insertions 
    beforeSeq = True
    # The current location in the slave sequence
    currentPlace = 0 
    # iterate over the slave sequence comparing each
    # character in the alignment to the master sequence.  
    for char in slave:
	# Break if the end of the alignment is reached
	if currentPlace == endPoint: break
	currentPlace +=1
	if char in string.ascii_uppercase:
	    # This is the situation where the alignments matched
	    # it also handles the start of the alignment, by switching
	    # beforeSeq to false. 
	    result += subst[master[masterCount] + slave[slaveCount]]
	    beforeSeq = False
	    slaveCount += 1
	    masterCount += 1
	    masterGap = False
	    slaveInsertion = False
	    continue
	if char == "-":
	    # This is the instance where the slave is in a gap
	    # If it is the first gap then an open cost is applied
	    # and the boolean masterGap is turned to true. 
	    # If the slave is in an insertion then a double gap
	    # is applied, and the boolean slaveInsertion is turned
	    # to false. 
	    if beforeSeq: continue
	    if masterGap:
		result += extendCost 
	    if slaveInsertion:
		result += doubleCost
	    if not masterGap and not slaveInsertion:
		result += openCost
	    masterGap = True
	    slaveInsertion = False
	    slaveCount += 1
	    masterCount += 1
	    continue
	if char in string.ascii_lowercase:
	    # This is the instance where the slave is in an insertion
	    # If it is the first insertion then an open cost is applied
	    # and the boolean slaveInsertion is turned to true. 
	    # If the slave is in a gap then a double gap
	    # is applied, and the boolean masterGap is turned
	    # to false. 
	    if beforeSeq: continue
	    if masterGap:
		result += doubleCost
	    if slaveInsertion:
		result += extendCost
	    if not masterGap and not slaveInsertion:
		result += openCost
	    masterGap = False
	    slaveInsertion = True
	    slaveCount += 1
	    continue
    return result

def scoreA2MGlobal(master, slave, subst, doubleCost, extendCost, openCost):
    """
    INPUT:
	master - A string, this will act as the reference
	slave - A string, this will act as a the query
	subst - A collections.Counter object, it keeps
	    substitution matrix scores for amino acid pairs. 
	    The amino acid pairs are represented as a string, 
	    the score for each pair is stored as an integer. 
	doubleCost - An integer, the cost penalty for a double
	    gap.
	extendCost - An integer, the cost penalty for extending
	    a gap.
	openCost - An integer, the cost penalty for opening a
	    new gap.
    OUTPUT:
	result - An integer, the alignment score associated  
	    with the master slave strings. 
    """
    # The result is initialized to 0, scores and penalties
    # will modify this value. 
    result = 0
    # Turn the string into a list so that it can be directly
    # indexed
    master = list(master.replace(" " , ""))
    slave = list(slave.replace(" " , ""))
    # This boolean will keep track of whether the sequence
    # is in a master gap
    masterGap = False
    # This boolean will keep track of whether the sequence
    # is in a slave insertion
    slaveInsertion = False
    # The masterCount and slaveCount integers will keep track
    # of the location in the master and slave sequences. 
    masterCount = 0
    slaveCount = 0
    for char in slave:
	if char in string.ascii_uppercase:
	    # This is the situation where the alignments matched
	    result += subst[master[masterCount] + slave[slaveCount]]
	    slaveCount += 1
	    masterCount += 1
	    masterGap = False
	    slaveInsertion = False
	    continue
	if char == "-":
	    # This is the instance where the slave is in a gap
	    # If it is the first gap then an open cost is applied
	    # and the boolean masterGap is turned to true. 
	    # If the slave is in an insertion then a double gap
	    # is applied, and the boolean slaveInsertion is turned
	    # to false. 
	    if masterGap:
		result += extendCost 
	    if slaveInsertion:
		result += doubleCost
	    if not masterGap and not slaveInsertion:
		result += openCost
	    masterGap = True
	    slaveInsertion = False
	    slaveCount += 1
	    masterCount += 1
	    continue
	if char in string.ascii_lowercase:
	    # This is the instance where the slave is in an insertion
	    # If it is the first insertion then an open cost is applied
	    # and the boolean slaveInsertion is turned to true. 
	    # If the slave is in a gap then a double gap
	    # is applied, and the boolean masterGap is turned
	    # to false. 
	    if masterGap:
		result += doubleCost
	    if slaveInsertion:
		result += extendCost
	    if not masterGap and not slaveInsertion:
		result += openCost
	    masterGap = False
	    slaveInsertion = True
	    slaveCount += 1
	    continue
    return result

def alignFile(inputFile, extendCost, doubleCost, openCost, alphabet, subst, isGlobal, score_file):
    """
    INPUT:
	inputFile - A file like object, the file should be
	    opened and in fasta format. 	
	extendCost - An integer, the cost penalty for extending
	    a gap.
	doubleCost - An integer, the cost penalty for a double
	    gap.
	openCost - An integer, the cost penalty for opening a
	    new gap.
	alphabet - The allowed alphabet, all characters that
	    have values in the substitution matrix. 
	subst - A collections.Counter object, it keeps
	    substitution matrix scores for amino acid pairs. 
	    The amino acid pairs are represented as a string, 
	    the score for each pair is stored as an integer. 
	isGlobal - A boolean, it is used to differentiate
	    between local and global alignments.
	score_file - A boolean, it is used primarily for 
	    debugging. The scores for each generated alignment
	    are printed below the alignment on a new line. 
    OUTPUT: 
	none
    This function prints to sys.stdout, it provides no internal 
    output.  The user specifies a global or local alignment
    and provides an input fasta file. The first sequence in the
    fasta file is the master sequence, all other sequences 
    are considered slave sequences.  The slave sequences are
    aligned to the master sequence. The new alignment is 
    printed to sys.stdout in fasta A2M format. If the value
    score_file is true then the score for each generated 
    alignment is printed on a new line beneath each fasta
    sequence. Therefore if the value for score_files is true
    the output will not be in fasta format. 
    """
    # This boolean is used to identify the first sequence in a 
    # file. The first sequence becomes the master sequence. 
    masterSeq = True
    # master will hold the master sequence
    master = ""
    # initiate a locAlignment class and a globAlignment class
    locAlignment = local_aligner(subst, openCost, extendCost ,doubleCost)
    globAlignment = global_aligner(subst, openCost, extendCost, doubleCost)
    # iterate over every sequence in the input fasta file.
    # Store the first sequence as the master sequence, all 
    # subsequent sequences are treates as slaves. Each slave
    # sequence is aligned to the master sequence, the 
    # optimal alignment is printed in A2M format. 
    for fasta in fastaFastqParser.read_fasta(inputFile, None,
		False, False): 
	if masterSeq:
	    # Handle the master sequence, check that each
	    # character is allowed before appending it to 
	    # master. After master is generated print the
	    # header and master. 
	    masterSeq = False
	    potentialMaster = fasta.sequence
	    for char in fasta.sequence:
	        if char in alphabet:
		     master += char
	    print (">"+fasta.identifier + fasta.comment)
	    print (master)
	else:
	    # Handle the slave sequences. First check that each
	    # character is allowed before appending it to slave
	    slave = ""
	    for char in fasta.sequence:
	        if char in alphabet:
		     slave += char
	    if isGlobal:
		# For global alignments, generate the matrix and
		# traceback from the last cell. Print the header
		# and the optimal alignment in A2M format. 
		print (">"+fasta.identifier + fasta.comment)
		globAlignment.score_align(master, slave)
		alignedSeq = globAlignment.traceback()
		print (alignedSeq)
		if score_file:
		    # Print the A2M score for the generated alignment
		    # and the master sequence. Primarily for debugging,
		    # this corrupts the output's fasta format. 
		    print (scoreA2MGlobal (master, alignedSeq, subst, doubleCost, extendCost, openCost))
 	    else:
		# For local alignments, generate the matrix and
		# traceback from the highScore. Print the header
		# and the optimal alignment in A2M format. 
		print (">"+fasta.identifier + fasta.comment)
		locAlignment.score_align(master, slave)
		alignedSeq = locAlignment.traceback()
		print (alignedSeq)
		if score_file:
		    # Print the A2M score for the generated alignment
		    # and the master sequence. Primarily for debugging,
		    # this corrupts the output's fasta format. 
		    print (scoreA2MLocal (master, alignedSeq, subst, doubleCost, extendCost, openCost))
