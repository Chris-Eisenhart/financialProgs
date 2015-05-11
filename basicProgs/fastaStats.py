s#!/usr/bin/env python2.7
# fastaStats
# Chris Eisenhart 05/11/2015 
# ceisenha@ucsc.edu/ceisenhart@soe.ucsc.edu 
"""
This program runs on FASTA files and spits out some statistics. The largest contig, 
average contig size, and total bases in the Fasta file will be reported. 
USAGE: 
  python fastaStats.py < inputfile > outputfile 
"""

from __future__ import print_function  
import  sys, operator, fileinput, collections, string, os.path
import  re, argparse, fastFunctions

alphabet="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'"
# The accepted characters in a word. A word is defined as any
# collection of characters from this alphabet. 
def parseArgs(args): 
    """
    INPUT:
	args: The command line information at compile time. 
    OUTPUT: 
	args: A tuple with five elements. 
    Parses the user supplied commands. The result is a tuple
    with five elements, one for each option. The options
    are descend, ascend, alphabet, inputFile and outputFile.  
    """
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument ("--inputFile",
    help = " Specifies that the text should be read in from an"
		" input file",
    action = "store")
    parser.add_argument ("--outputFile",
    help = " Specifies that the word and count pairs should"
		" be printed to an output file",
    action = "store")
    options = parser.parse_args()
    return options


def main(args):
    """
    """
    options = parseArgs(args)
    # User options are stored here
    starts = ["ATG"]
    stops = ["TAG", "TAA", "TGA"]
    inputFile = sys.stdin
    # Read from sys.stdin as a default.
    if options.inputFile is not None:
        inputFile = open(options.inputFile, 'r')
    outputFile = sys.stdout
    totalBases = 0
    totalSeqs = 0
    longestContig = 0
    for fasta in fastFunctions.readFasta (inputFile, alphabet, False, True):
        totalBases += len(fasta.sequence)
        totalSeqs += 1
        if len(fasta.sequence) > longestContig: longestContig = len(fasta.sequence)
    print ("There are %i sequences with an average of %i bases. The longest contig is %i bases"
                " There are %i bases total"% (totalSeqs, totalBases/totalSeqs, longestContig, totalBases))

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
