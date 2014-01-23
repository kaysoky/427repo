import sys
import os
import argparse
import random

import numpy
from numpy import zeros, array, argmax, size, unravel_index

"""
A local file containing the score matrix 
for alignment of various proteins.
"""
BLOSUM = 'BLOSUM62.txt'

"""
A map from protein letter to index into the BLOSUM score matrix
"""
PROT_INDEX = {}

"""
Letter within the PROT_INDEX used to match all other letters
"""
WILDCARD = "*"

"""
The cost of aligning a gap with any letter
"""
GAP_COST = -4

"""
The maximum length of the protein label used when printing alignments
"""
LABEL_LENGTH = 6

"""
The number of letters per line of sequence printing
"""
SEQUENCE_FRAGMENT = 60

"""
File extension that triggers some additional processing
The first line is ignored and whitespace is stripped out
"""
FASTA = '.fasta'
    
# Read in the score matrix
with open(BLOSUM, 'r') as f:
    content = f.readlines()
content = [line.strip() for line in content]
content = filter(lambda x: x[0] != "#", content)

# The first line of the score matrix holds labels
letters = content[0].split()
for iter in range(len(letters)):
    PROT_INDEX[letters[iter].upper()] = iter
    PROT_INDEX[letters[iter].lower()] = iter

# Initialize the score matrix
BLOSUM = zeros((len(letters), len(letters)))
for iter in range(len(content) - 1):
    line = content[iter + 1]
    line = line.split()[1:(1 + len(letters))]
    BLOSUM[iter:] = [int(token) for token in line]

def do_align(sequenceA, sequenceB):
    """
    Takes two strings and performs Smith-Waterman local alignment
    Returns the calculated score matrix
    """
    
    # Pad both sequences with a leading space
    # Doing so aligns the letters of the sequence with the indices of the matrix
    sequenceA = " " + sequenceA
    sequenceB = " " + sequenceB
    
    # Initialize the local alignment matrix
    alignments = zeros((len(sequenceA), len(sequenceB)))

    # Perform the alignment, column by column
    for a in range(1, len(sequenceA)):
        for b in range(1, len(sequenceB)):
            # Find the associated score in the BLOSUM matrix
            paths = _calculate_costs(alignments, (a, b), (sequenceA[a], sequenceB[b]))
            alignments[a, b] = max(paths)
    return alignments
            
def do_traceback(alignments, sequenceA, sequenceB, rowColumn=None):
    """
    Performs a traceback
    If rowColumn (tuple or array of 2 values) is provided, 
        then the trace will start from that index in the matrix
    
    :return: A tuple of 3 values (an aligned sequence for A, 
                                  a comparison sequence, 
                                  and an aligned sequence for B)
    """
    
    # Pad both sequences with a leading space
    # Doing so aligns the letters of the sequence with the indices of the matrix
    sequenceA = " " + sequenceA
    sequenceB = " " + sequenceB
    assert len(sequenceA) == size(alignments, 0)
    assert len(sequenceB) == size(alignments, 1)
    
    # Get the starting indices
    if rowColumn is None or len(rowColumn) != 2:
        traceRow, traceCol = unravel_index(argmax(alignments), 
                (size(alignments, 0), size(alignments, 1)))
    else:
        traceRow, traceCol = rowColumn
    
    tracedA = ''
    middle = ''
    tracedB = ''
    while alignments[traceRow, traceCol] != 0:
        paths = _calculate_costs(alignments, (traceRow, traceCol), (sequenceA[traceRow], sequenceB[traceCol]))
        path = argmax(paths)
        
        if path == 0:
            tracedA = sequenceA[traceRow] + tracedA
            tracedB = sequenceB[traceCol] + tracedB
            if sequenceA[traceRow] == sequenceB[traceCol]:
                middle = sequenceA[traceRow] + middle
            elif alignments[traceRow, traceCol] - alignments[traceRow - 1, traceCol - 1] > 0:
                middle = '+' + middle
            else:
                middle = ' ' + middle
                
            traceRow -= 1
            traceCol -= 1
            
        elif path == 1:
            tracedA = sequenceA[traceRow] + tracedA
            tracedB = '-' + tracedB
            middle = ' ' + middle
            
            traceRow -= 1
        
        elif path == 2:
            tracedA = '-' + tracedA
            tracedB = sequenceB[traceCol] + tracedB
            middle = ' ' + middle
            
            traceCol -= 1
        
        elif path == 3:
            print "Impossible case"
            sys.exit(0)
            
    return (tracedA, middle, tracedB)
    
def _calculate_costs(alignments, rowColumn, letters):
    """
    Helper for calculating the cost of a particular cell
    
    :param alignments: The score matrix to use
    :param rowColumn:  A tuple of the indices to calculate
    :param letters:    A tuple of the letters to compare
    :return:           An array of 4 values
                       Index 0 -> cost from above left
                       Index 1 -> cost from above
                       Index 2 -> cost from left
                       Index 3 -> zero
    """
    
    matchScore = GAP_COST
    if all([letter in PROT_INDEX for letter in letters]):
        matchScore = BLOSUM[PROT_INDEX[letters[0]], PROT_INDEX[letters[1]]]
    elif letters[0] == letters[1]:
        matchScore = BLOSUM[PROT_INDEX[WILDCARD], PROT_INDEX[WILDCARD]]
    
    # Calculate the scores
    a, b = rowColumn
    return array([
        alignments[a - 1, b - 1] + matchScore, 
        alignments[a - 1, b] + GAP_COST, 
        alignments[a, b - 1] + GAP_COST, 
        0])

def print_alignments(labels, alignments, originals):
    """
    Neatly prints out the alignments
    :param labels: a tuple of two values, only the first 5 letters will be used
    :param alignments: a tuple of 3 values (sequence A, middle, sequence B)
    :param originals: a tuple of 2 values containing the orignal sequences
    """
    
    labelA, labelB = [label[0:LABEL_LENGTH] if len(label) > LABEL_LENGTH else label for label in labels]
    traceA, middle, traceB = alignments
    sequenceA, sequenceB = originals
    
    assert(len(traceA) == len(traceB))
    
    # Find the starting index for both sequences
    # Add one to make this a 1-based index
    startA = sequenceA.find(traceA.replace('-', '')) + 1
    startB = sequenceB.find(traceB.replace('-', '')) + 1
    
    # Do the printing
    for i in range(0, len(traceA), SEQUENCE_FRAGMENT):
        subTraceA = traceA[i:(i + SEQUENCE_FRAGMENT)]
        subTraceB = traceB[i:(i + SEQUENCE_FRAGMENT)]
        print "%*s: %*s %s" % (LABEL_LENGTH, labelA, LABEL_LENGTH, str(startA).center(LABEL_LENGTH, ' '), subTraceA)
        print "%*s %s" % (LABEL_LENGTH * 2 + 2, '', middle[i:(i + SEQUENCE_FRAGMENT)])
        print "%*s: %*s %s" % (LABEL_LENGTH, labelB, LABEL_LENGTH, str(startB).center(LABEL_LENGTH, ' '), subTraceB)
        print ''
        startA += len(subTraceA.replace('-', ''))
        startB += len(subTraceB.replace('-', ''))
    
def shuffle_string(text):
    """
    Explicitly implements shuffle to match assignment specs
    Normally, this would be sufficient:
        ''.join(random.shuffle(list(text)))
    """
    
    text = list(text)
    for i in range(len(text)):
        randIndex = random.randint(i, len(text) - 1)
        temp = text[i]
        text[i] = text[randIndex]
        text[randIndex] = temp
    return ''.join(text)
        
def calculate_empirical_probability(sequenceA, sequenceB, optimal, num, isVerbose=False):
    """
    Aligns sequenceA to some number of random permutations of sequenceB
    """
    
    betterCount = 0
    for i in range(num):
        permutation = shuffle_string(sequenceB)
        scores = do_align(sequenceA, permutation)
        best = numpy.amax(scores)
        if isVerbose:
            print "Permutation %d score: %d" % (i, best)
        if best >= optimal:
            betterCount += 1
    
    if num > 0:
        return float(betterCount) / num
    return float(betterCount)

def process_fasta(text):
    """Removes the first line and removes newlines"""
    lines = text.split('\n')
    return ''.join(lines[1:])
    
def do_main(sequenceA, sequenceB, labelA, labelB, isVerbose=False, num=0):
    # Print the score matrix
    scores = do_align(sequenceA, sequenceB)
    if isVerbose:
        print "Score matrix:"
        print scores
        print ''
    
    # Print the alignments
    alignments = do_traceback(scores, sequenceA, sequenceB)
    print "Alignment:"
    print_alignments((labelA, labelB), alignments, (sequenceA, sequenceB))
    
    # Print the optimal score
    optimal = numpy.amax(scores)
    print "Optimal score: %d\n" % optimal
    
    # Calculate the empirical probability
    if num > 0:
        probability = calculate_empirical_probability(sequenceA, sequenceB, optimal, num, isVerbose)
        print "Empirical probability: %f\n" % probability

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Performs Smith-Waterman local alignment on sequences found in two files')
    parser.add_argument('sequenceA', type=str)
    parser.add_argument('sequenceB', type=str)
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('-n', type=int, default=0)
    args = parser.parse_args()

    # Read the two sequences in as strings
    with open(args.sequenceA) as f:
        sequenceA = f.read().strip()
    with open(args.sequenceB) as f:
        sequenceB = f.read().strip()
        
    # Handle FASTA files
    if os.path.splitext(args.sequenceA)[1] == FASTA:
        sequenceA = process_fasta(sequenceA)
    if os.path.splitext(args.sequenceB)[1] == FASTA:
        sequenceB = process_fasta(sequenceB)
    
    do_main(sequenceA, sequenceB, args.sequenceA, args.sequenceB, args.verbose, args.n)
