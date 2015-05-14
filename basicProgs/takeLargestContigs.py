#!/usr/bin/env python2.7
# takeLargestContigs
# Chris Eisenhart 05/14/2015 
# ceisenha@ucsc.edu/ceisenhart@soe.ucsc.edu 
"""
This program runs on a fasta file. The largest 30 sequences are printed to the output file provided
USAGE: 
  python takeLargestContigs.py < inputfile > outputfile 
"""

from __future__ import print_function  
import  sys, operator, fileinput, collections, string, os.path
import  re, fastFunctions

alphabet="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'"
# The accepted characters in a word. A word is defined as any
# collection of characters from this alphabet. 

def biggerThan(length,bigList):
    """
    Sees if the length belongs on the list
    """
    for item in bigList:
        if length > item: return True 
    return False

def main(args):
    """
    Go over each sequence in the fasta file. Keep a list of the 30 largest sequences, for each 
    sequence seen append it to the list, then sort the list.  Take the 30 largest sequences
    and go onto the next sequence.  After the file is processed iterate over the list and print
    out each sequence in .fasta format. 
    """
    inputFile = sys.stdin
    outputList = []
    for fasta in fastFunctions.readFasta (inputFile, alphabet, False, True):
        outputList.append((len(fasta.sequence),fasta.sequence, fasta.header,fasta.extraData))
        outputList = sorted(outputList, key = operator.itemgetter(0), reverse = True)
        outputList = outputList[:30]
    i = 1
    for length, sequence, header, extraData in outputList:
        print(">contig%i %s%s length %i"%(i,header,extraData, length))
        print(sequence)
        i +=1
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
