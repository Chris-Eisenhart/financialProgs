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
    totalBases = 0.0
    totalSeqs = 0
    longestContig = 0 
    contigs1k = 0
    contigs3k = 0
    contigs5k = 0
    contigs10k = 0
    contigs20k = 0
    contigs1kTotal = 0
    contigs3kTotal = 0
    contigs5kTotal = 0
    contigs10kTotal = 0
    contigs20kTotal = 0 
    for fasta in fastFunctions.readFasta (inputFile, alphabet, False, True):
        totalBases += len(fasta.sequence)
        totalSeqs += 1
        if len(fasta.sequence) > 1000:
            contigs1k +=1 
            contigs1kTotal += len (fasta.sequence)
        if len(fasta.sequence) > 3000: 
            contigs3k +=1 
            contigs3kTotal += len (fasta.sequence)
        if len(fasta.sequence) > 5000:
            contigs5k +=1 
            contigs5kTotal += len (fasta.sequence)
        if len(fasta.sequence) > 10000: 
            contigs10k +=1 
            contigs10kTotal += len (fasta.sequence)
        if len(fasta.sequence) > 20000:
            contigs20k +=1 
            contigs20kTotal += len (fasta.sequence)
        if len(fasta.sequence) > longestContig:
            longestContig = len(fasta.sequence)
    print ("There are %i sequences with an average of %i bases. The longest contig is %i bases"
                " There are %i bases total"% (totalSeqs, totalBases/totalSeqs, longestContig, totalBases))
    print ("There are %i sequences of length 1,000 or greater, the bases in these sequences account for %f0.1 of the genome." % (contigs1k, float(contigs1kTotal/totalBases )))
    print ("There are %i sequences of length 3,000 or greater, the bases in these sequences account for %f0.1 of the genome." % (contigs3k, float(contigs3kTotal/totalBases )))
    print ("There are %i sequences of length 5,000 or greater, the bases in these sequences account for %f0.1 of the genome." % (contigs5k, float(contigs5kTotal/totalBases )))
    print ("There are %i sequences of length 10,000 or greater, the bases in these sequences account for %f0.1 of the genome." % (contigs10k, float(contigs10kTotal/totalBases )))
    print ("There are %i sequences of length 20,000 or greater, the bases in these sequences account for %f0.1 of the genome." % (contigs20k, float(contigs20kTotal/totalBases )))

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
