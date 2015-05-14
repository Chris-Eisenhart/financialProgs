#!/usr/bin/env python2.7
# grabContigs
# Chris Eisenhart 05/11/2015 
# ceisenha@ucsc.edu/ceisenhart@soe.ucsc.edu 
"""
This program runs on a fasta file. Contigs of length 5000 or greater are coppied to the output file
USAGE: 
  python fastaStats.py < inputfile > outputfile 
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
    """
    inputFile = sys.stdin
    outputList = []
    for fasta in fastFunctions.readFasta (inputFile, alphabet, False, True):
        outputList.append((len(fasta.sequence),fasta.sequence, fasta.header,fasta.extraData))
        outputList = sorted(outputList, key = operator.itemgetter(0), reverse = True)
        outputList = outputList[:30]
    for length, sequence, header, extraData in outputList:
        print(">%s%s length %i"%(header,extraData, length))
        print(sequence)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
