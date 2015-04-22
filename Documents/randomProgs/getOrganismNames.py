#!/usr/bin/env python2.7
# getOrganismNames
# Nelson Seilhan, Chris Eisenhart 02/16/2015
"""
This program is designed to parse the NCBI Genbank summary format. The user provides
an input file in NCBI Genbank summary format and the program identifies the total
number of base pairs present. Additionally a list of unique genera/species combinations
is printed. 
"""
from __future__ import print_function  
import  sys, operator, fileinput, collections, string, os.path
import  re, argparse, random

# USAGE:
# Check that the program, input file are in the current directory using the '-ls' command.
# Verify Python is installed using 'python' command. To run the program use the command;
#	python getOrganismNames < inputFile 
# To print the output to a new file; 
# 	python getOrganismNames < inputFile > outputFile
# To use the -dense option and print only genera; 
# 	python getOrganismNames < inputFile > outputFile --dense
#			OR
# 	python getOrganismNames < inputFile > outputFile -d

def parseArgs(args): 
    """
    Set the specifications for user provided options. The following 
    options are supported, 
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument ("--dense", "-d",
    help = " Print the output in dense format, showing only unique genera.",
    action = "store_true")
    parser.add_argument ("--verbose", "-v",
    help = " Print individual simulation statistics and all default conditions",
    action = "store_true")
    parser.set_defaults (verbose = False)
    parser.set_defaults (dense = False)
    options = parser.parse_args()
    return options

def main(args):
    """
    Reads in the user commands and options, opening
    any files if necessary otherwise reading from sys.stdin. 
    """
    options = parseArgs(args)
    sumTotal = 0 
    scores = collections.Counter()
    resultList = []
    denseResultList = []
    for line in sys.stdin:
	splitLine = line.split()
    	splitLineList = list(splitLine)
	location = 0
	for item in splitLine:
	    if splitLineList[location] == "ORGANISM":
		if (splitLineList[location+1],splitLineList[location+2]) in resultList: continue
		else: resultList.append((splitLineList[location+1],splitLineList[location+2]))
		if (splitLineList[location+1]) in denseResultList: continue
		else: denseResultList.append(splitLineList[location+1])
	    if splitLineList[location] == "bp" and splitLineList[location+1] == "DNA": 
		sumTotal += int(splitLineList[location-1])
	    location += 1
    if options.dense:
	for genera in sorted(denseResultList):
	    print(genera)
    else:
    	for genera, species in sorted(resultList):
	    print (genera, species)
    print ("The total number of bases is %s" %(sumTotal))

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
