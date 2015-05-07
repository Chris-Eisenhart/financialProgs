#!/usr/bin/env python2.7
# findORFs
# Chris Eisenhart 11/10/2014 
# ceisenha@ucsc.edu/ceisenhart@soe.ucsc.edu  Third Version
"""
"""

from __future__ import print_function  
import  sys, operator, fileinput, collections, string, os.path
import  re, argparse, fastFunctions

class frame:
    """
    This class holds all the start/stop codon information for a given codon frame in a sequence of DNA. 
    """
    def __init__(self):
        self.starts = [] # All start positions for the frame, a list of integers
        self.ends = [] # All end positions for the frame, a list of integers
        self.rStarts = [] # All start positions for the reverse complement frame, a list of integers
        self.rEnds = [] # All stop positions for the reverse complement frame, a list of integers


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

def printFrame(frameNum,frame):
    """
    """
    print("For frame", frameNum,"starts",frame.starts,"ends", frame.ends,"rStarts", frame.rStarts,"rEnds", frame.rEnds)

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
    for fasta in fastFunctions.readFasta (inputFile, alphabet, False, True):
      i = 0; 
      allFrames = {0:frame(), 1:frame(), 2:frame()}
      for char in fasta.sequence: 
          codon = fasta.sequence[i:i+3]
          complement_table = string.maketrans("ACGT", "TGCA")
          rCodon = codon[::-1].translate(complement_table)
          if codon in stops:
              allFrames[i%3].ends.append(i)
          if codon in starts:
              allFrames[i%3].starts.append(i)
          if rCodon in stops:
              allFrames[i%3].rEnds.append(i)
          if rCodon in starts:
              allFrames[i%3].rStarts.append(i)
          i += 1
      for key,value in allFrames.iteritems(): 
          printFrame(key,value)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
