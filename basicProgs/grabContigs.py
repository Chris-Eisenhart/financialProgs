#!/usr/bin/env python2.7
# grabContigs
# Chris Eisenhart 05/11/2015 
# ceisenha@ucsc.edu/ceisenhart@soe.ucsc.edu 
"""
This program runs on a fasta file. Contigs of length 10000 or greater are coppied to the output file
USAGE: 
  python fastaStats.py < inputfile > outputfile 
"""

from __future__ import print_function  
import  sys, operator, fileinput, collections, string, os.path
import  re, fastFunctions

alphabet="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'"
# The accepted characters in a word. A word is defined as any
# collection of characters from this alphabet. 

def main(args):
    """
    """
    inputFile = sys.stdin
    i = 1 
    for fasta in fastFunctions.readFasta (inputFile, alphabet, False, True):
        if len(fasta.sequence) > 10000: # Changing this number updates the contig size cutoff 
            print (">contig%i" %(i))
            print (fasta.sequence)
            i += 1
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
