import argparse

import numpy
from numpy import zeros, array, argmax

"""
A local file containing the score matrix 
for alignment of various proteins.
"""
BLOSUM = "BLOSUM62.txt"

"""
A map from protein letter to index into the BLOSUM score matrix
"""
PROT_INDEX = {}

"""
The cost of aligning a gap with any letter
"""
GAP_COST = -4
    
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
# Get rid of the wildcard "*"
del PROT_INDEX['*']

# Initialize the score matrix
BLOSUM = zeros((len(letters), len(letters)))
for iter in range(len(content) - 1):
    line = content[iter + 1]
    line = line.split()[1:(1 + len(letters))]
    BLOSUM[iter:] = [int(token) for token in line]

def do_align(sequenceA, sequenceB):
    """
    Takes two strings and performs Smith-Waterman local alignment
    Returns an optimal alignment
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
            matchScore = GAP_COST
            if sequenceA[a] in PROT_INDEX and sequenceB[b] in PROT_INDEX:
                matchScore = BLOSUM[PROT_INDEX[sequenceA[a]], PROT_INDEX[sequenceB[b]]]
            
            # Calculate the scores
            paths = array([
                alignments[a - 1, b - 1] + matchScore, 
                alignments[a - 1, b] + GAP_COST, 
                alignments[a, b - 1] + GAP_COST, 
                0])
                
            alignments[a, b] = max(paths)
    
    # Perform a trace-back
    print alignments
    ##TODO

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description='Performs Smith-Waterman local alignment on sequences found in two files')
    # parser.add_argument('scoreMatrix', type=file)
    parser.add_argument('sequenceA', type=file)
    parser.add_argument('sequenceB', type=file)
    args = parser.parse_args()

    # Read the two sequences in as strings
    sequenceA = args.sequenceA.read().strip()
    sequenceB = args.sequenceB.read().strip()
    do_align(sequenceA, sequenceB)
