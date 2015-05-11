#!/usr/bin/env python2.7
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
import  re, fastFunctions

alphabet="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'"
# The accepted characters in a word. A word is defined as any
# collection of characters from this alphabet. 

def main(args):
    """
    """
    inputFile = sys.stdin
    outputFile = sys.stdout
    totalBases = 0
    totalSeqs = 0
    longestContig = 0 
    contigs1k = 0
    contigs3k = 0
    contigs5k = 0
    contigs7k = 0
    contigs9k = 0
    for fasta in fastFunctions.readFasta (inputFile, alphabet, False, True):
        totalBases += len(fasta.sequence)
        totalSeqs += 1
        if len(fasta.sequence) > 1000: contigs1k +=1 
        if len(fasta.sequence) > 3000: contigs3k +=1 
        if len(fasta.sequence) > 5000: contigs5k +=1 
        if len(fasta.sequence) > 7000: contigs7k +=1 
        if len(fasta.sequence) > 9000: contigs9k +=1 
        if len(fasta.sequence) > longestContig:
            longestContig = len(fasta.sequence)
    print ("There are %i sequences with an average of %i bases. The longest contig is %i bases"
                " There are %i bases total"% (totalSeqs, totalBases/totalSeqs, longestContig, totalBases))
    print ("There are %i contigs of length 1000 or greater" % (contigs1k))
    print ("There are %i contigs of length 3000 or greater" % (contigs3k))
    print ("There are %i contigs of length 5000 or greater" % (contigs5k))
    print ("There are %i contigs of length 7000 or greater" % (contigs7k))
    print ("There are %i contigs of length 9000 or greater" % (contigs9k))
if __name__ == "__main__" :
    sys.exit(main(sys.argv))
