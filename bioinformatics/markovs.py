#!/usr/bin/env python2.7
# markovs
# Chris Eisenhart  ceisenha@ucsc.edu/ceisenhart@soe.ucsc.edu 10.28.2014 
"""
This program contains functions intended for use with the Python programs 
count-kmers and coding-cost.  All functions for the two programs will be stored
here. These functions are accessed in other python programs by prefixing this 
program's name, for example markovs.sort_counts(counts, True) is a valid 
function call.  
"""
import sys, operator, fileinput, string, os.path
import re, argparse, fastFunctions, itertools
import collections, math



def sort_counts(counts, descend):
    """
    Sorts a collection.Counter object by the key, the output will be in 
    alphabetical order.  Returns a sorted list of the collection.Counter keys,
    the returned list is not a collections.Counter object. If the descence
    option is True the output will be printed in descending order.  
    """
    sortedCounts = sorted(counts.items(), key = operator.itemgetter(descend))
    return sortedCounts


def print_counter(counts):
    """
    Takes in a collections.Counter and iterates through it. Each
    key value pair is printed to a stdout as a single line with a space separating 
    the two values. This is used in count-kmers to print the kmer counts. 
    """
    for key, value in counts:
        print(str(key) + " " + str(value))


def print_scores(scores, kmerFile, fastaFile, k, alpha):
    """
    Takes in a scores tuple. Prints the tuple to sys.stdout with 
    corresponding name tags. 
    """
    print ("Train: " + kmerFile + " Test: " + fastaFile) 
    print ("Order: " + k + " Alphabet: " + alpha)
    print ("Total encoding cost for file: " + str(scores[0]))
    print ("Average cost per sequence: " + str(scores[1]))
    print ("Average cost per character: " + str(scores[2]))

def get_alpha(counts):
    """
    Takes in a collections.Counter and iterates through it.  Each
    key is iterated over, any new character is added to the result.
    The result contains all characters seen in the keys of counts,
    these characters are returned as the alphabet. 
    """
    alpha = "" # alpha is a string, each character will be added onto it
    for entry in counts.items():
        for char in entry[0]:
            if char not in alpha:
                alpha += char
    return alpha 

def read_counter_and_kmerSize(inputFile):
    """
    Takes in an opened file like object.  The file must have the format that
    the kmer and count are listed on a single line, separated by a space. 
    Returns a collections.Counter object with all key value pairs from the 
    input file and the size of the kmers. Additionally the kmer size is returned,
    the return value is a tuple, counts can be accessed with result[0] and 
    the length with result[1]
    """
    counts = collections.Counter()
    # counts["temp"] returns the number associated with temp
    for line in inputFile:
        splitLine = line.split()
        # splitLine[0] refers to the first half of the line (key),
        # splitLine[1] refers to the second half of the line (value). 
        counts[splitLine[0]] = int(splitLine[1])
    return counts, len(splitLine[0])

def count_kmers(inputFile, counts, k, alpha, verbose): 
    """ 
    Takes in an opened input file, a collections.Counter, and a integer k.  
    K corresponds to the order of the markov chain, as well as the kmer size. 
    The program counts the occurences of each kmer of the specified size.
    The collections.Counter is returned after the fasta file has been fully
    processed. The values in the Counter can be accessed as follows; given a 
    kmer "TTA", result["TTA"] returns the count of "TTA". 
    """
    if alpha is None:
    # If no alphabet is provided use all ascii upper case letters 
        alpha = string.ascii_uppercase
    for fasta in fastFunctions.readFasta(inputFile,alpha,verbose,True):
    # fasta is a fasta object, as defined in fastFunctions
        if len(fasta.sequence)== 0: continue
        # Throw away empty sequences.
	fasta = add_prefix_suffix(fasta, k)
        # Add prefix and suffix characters
        for start in range(len(fasta.sequence) - (k-1)):
        # Iterate over all kmers, incrementing the count each time
        # a kmer is seen. For example with 3-mers the last kmer is 
        # 'AA$' where A is any char in the alphabet and $ is the end char. 
            counts[fasta.sequence[start:start+k]] += 1
            # fasta.sequence[start:start+k] is a single kmer string.
            # These kmers are the key in counts, so counts["kmer"]
            # returns the count of "kmer". 
    return counts

def add_prefix_suffix(fasta, k):
    """
    Takes in a fasta file and the kmer size.  The fasta sequence will have
    k-1 $ characters prefixing it and k-1 * characters appended.  
    """
    prefix = ""
    suffix = ""
    # The prefix and suffix characters, stored internally as a strings
    if k == 1:
    # Handle the 0 order case
        prefix = "$"
        suffix = "*"
    for i in range(k - 1):
    # Handle all orders 1 and greater
        prefix += "$"
        suffix += "*"
    fasta.sequence += suffix
    # Append the suffix to the fasta sequence, the sequence now ends
    # TTAS** for a 3 mer. 
    prefix += fasta.sequence
    # Append the fasta sequence to the prefix string. The prefix string
    # now holds the full fasta sequence with suffix and prefix. 
    fasta.sequence = prefix
    # Set the fasta.sequence to point to the modified sequence. 
    return fasta 

  
def score_fasta_file(inputFile, scoreMatrix, alpha, k, verbose):
    """
    Takes in a Fasta input file and a score matrix.  Iterates over every
    fasta sequence in the fasta file and calculated the encoding cost.  The
    encoding cost for a character x using model M is -log2P_M(x), the encoding
    cost for a sequence is the summation of these individual scores.  Likewise
    the encoding cost for a fasta file is the summation of the encoding costs
    of each sequence. 
    """
    fileCost = 0
    # A float that will keep track of the total file cost
    sequenceCount = 0
    # An int that keeps count of the number of sequences seen
    charCount = 0
    # An int that keeps count of the total number of characters seen
    for fasta in fastFunctions.readFasta(inputFile,alpha,verbose,True):
        charCount += len(fasta.sequence) 
        sequenceCount += 1
        sequenceCost = 0
        if len(fasta.sequence) == 0: continue
        # Throw away the empty sequence 
        if k >> 1:
        # Append prefix and suffix characters, but handle the 0 order case
        # specifically. The 0 order model only has end characters added.
            fasta = add_prefix_suffix(fasta,k)
        else: 
            fasta.sequence += "*"
        for start in range(len(fasta.sequence)-(k - 1)): 
        # Iterate over all kmers, calculating the sequence cost and adding it to 
        # the running total.
            sequenceCost += scoreMatrix[fasta.sequence[start:start+k]]
        fileCost += sequenceCost
    return fileCost, fileCost/sequenceCount, fileCost/(charCount)

def make_score_matrix(counts, k):
    """
    Takes in a sorted collection.Counter() object (use sort_counts(counts,False)).
    The counter will be used to make a score matrix.
    Elements in the counter are keyed by kmer, stored as a string. The value is 
    the count of each kmer, in floats. It is assumed that all kmers are represented
    and have non 0 counts. K is the length of the kmers. 
    A collections.Counter() object is returned, where result["kmer"] returns the
     - log_2_ (probability) of the kmer appearing. 
    """
    kmerBlock = ""
    # A string that stores the name of a kmer block. A kmer block is a 
    # collection of kmers that begin with the same k-1 characters, where k
    # is the kmer length. 
    
    totalCount = 0
    # This is the total count for a kmer block, it will be uses to normalize the
    # kmer probabilities within kmer blocks. 
    
    allKmerBlocks = collections.Counter()
    # An intermediate hash table, this table stores all kmer blocks and the 
    # associated total count for a kmer block. For example given kmers TGB, TGC,
    #  and TGD, with counts 1, 2, and 3. The kmer block is  "TG" and the count
    # is 6. These values are stored in allKmerBlocks so that allKmerBlocks["TG"] 
    # will return 6.  
    for element in counts:
        if element[0] is "$": continue
        # This addresses the 0 order case, where the start char does not apply. 
        if kmerBlock != element[0][0:k-1]:
            allKmerBlocks[kmerBlock] = totalCount
            totalCount = 0
        kmerBlock = element[0][0:k-1]
        totalCount += element[1]
    allKmerBlocks[kmerBlock] = totalCount   
    resultMatrix = collections.Counter()    
    # The result hash table, stored internally as a collections.Counter()
    # resultMatrix["TTA"] returns the a probability value between 0 and 1.  
    # The values of keys that start with "TT" will sum to 1. This defines
    # the probability of any character given "TT".  Every kmer is assigned
    # a value.
    for element in counts: 
       if element[0][0:k-1] in allKmerBlocks:
           if element[0] == "$": continue
           # This addresses the 0 order case, where the start char does not apply. 
           resultMatrix[element[0]] = -math.log(float(element[1])/float(allKmerBlocks[element[0][0:k-1]]), 2)
    return resultMatrix

 
 
def pseudo_count(counts, pseudo, alpha, k):
    """
    Takes in a collections.Counter object, a float pseudo, the alphabet, and k
    mer size.  Every possible kmer will be determined from the alphabet, and given 
    a value of pseudo. These kmers and scores are stored in a collecition.Counter 
    object, which is merged with counts to create the result. 
    """
    alpha = alpha.translate(None,string.punctuation)
    # Remove the start and end characters, they will be handled specifically
    pseudoCounts = collections.Counter()
    # A collections.Counter object. All possible kmers will be inserted into 
    # pseudoCounts as keys, the value will be psuedo. This collections.Counter
    # object will be merged with the one provided to create the result.
    kmers = itertools.product(alpha, repeat = k)
    # Kmers is a list of tuples, each tuple is a possible kmer.  Every possible
    # kmer is represented in this list of tuples.  For example
    # ('A', 'B', 'C') corresponds to the string 'ABC'.  The tuples need to be
    # processed into strings before they can be entered as keys in psuedoCounts.
    for item in kmers:
        kmerTuple = str(item).translate(None,"(),'")
        # Turn the tuple into a string and remove tuple punctuation
        kmer = ""
        # Store the kmer as a string, and reset it each iteration. 
        for char in kmerTuple.split():
        # Remove the while spaces and add the characters to the kmer
            kmer += char
        pseudoCounts[kmer] = pseudo
   
    if k == 1: 
    # Handle the 0 order case. 
        pseudoCounts["$"] = pseudo
        pseudoCounts["*"] = pseudo
    
    # This block of code pseudo counts the prefix edge kmers. For example
    # consider a 2nd order markov chain  all prefix edge kmers are  of the form
    #  $$A and $AA where  A is any character in the alphabet.
    prefix = 1
    for i in range(k):
        if k-prefix > 0:
            for item in itertools.product(alpha,repeat = k-prefix):
                # This iterates over all possible kmers of size k - prefix
                kmerTuple = str(item).translate(None,"(),'")
                # Turn the tuple into a string and remove tuple punctuation
                # Reset the kmer
                kmer = ""
                for char in kmerTuple.split():
                # Remove the while spaces and add the characters to the kmer
                    kmer += char
                prefixChars = ""
                for int in range(prefix):
                # Add the correct number of prefix character 
                    prefixChars += "$"
                prefixChars += kmer
                # Append the kmer to the prefix characters, now prefixChars
                # is the completed kmer.  
                pseudoCounts[prefixChars] = pseudo
                # Insert the new kmer into the pseudo count collection, with a 
                # value of pseudo (defaults to 1). 
            prefix += 1
            # At this point all strings with k-1 alphabet characters have been 
            # processed, and we move onto k-2.  For example, a 3-mer, after one
            # pass all '$AA' kmers (where A is any character) have been added. Prefix
            # is incremented, and the next iteration will process all '$$A' kmers. 

    # This block of code pseudo counts the suffix edge kmers. For example
    # consider a 2nd order markov chain all suffix edge kmers are of the form
    #  A** and AA* where  A is any character in the alphabet.
    suffix = 1
    # This keeps track of how many suffix characters are needed
    for i in range(k):
        if k-suffix > 0:
            # for item in itertools.permutations(alpha,k-suffix):
            for item in itertools.product(alpha, repeat = k-suffix):
                # This iterates over all possible kmers of size k - prefix
                kmerTuple = str(item).translate(None,"(),'")
                # Turn the tuple into a string and remove tuple punctuation
                kmer = ""
                # Reset the kmer, the kmer is a string. 
                for char in kmerTuple.split():
                # Remove the white spaces and add the characters to the kmer
                    kmer += char
                suffixChars = ""
                # Stores the total suffix characters for the current edge level
                for int in range(suffix):
                    suffixChars += "*"
                kmer += suffixChars
                # Append the correct number of suffix characters
                pseudoCounts[kmer] = pseudo
                # Insert the new kmer into the pseudo count collection, with a 
                # value of pseudo (defaults to 1). 
            suffix += 1
            # At this point all strings with k-1 alphabet characters have been 
            # processed, and we move onto k-2.  For example, a 3-mer, after one
            # pass all 'AA*' kmers (where A is any character) have been added. Suffix
            # is incremented, and the next iteration will process all 'A**' kmers. 
    return counts + pseudoCounts #Merges the two collection.Counter objects
