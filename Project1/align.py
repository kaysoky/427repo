import argparse

import numpy
from numpy import zeros

"""
A local file containing the score matrix 
for alignment of various proteins.
"""
BLOSUM = "BLOSUM62.txt"

"""
A map from protein letter to index into the BLOSUM score matrix
"""
PROT_INDEX = {}

parser = argparse.ArgumentParser(
        description='Performs BLAST (Smith-Waterman) on sequences found in two files')
# parser.add_argument('scoreMatrix', type=file)
parser.add_argument('sequenceA', type=file)
parser.add_argument('sequenceB', type=file)
args = parser.parse_args()

# Read the two sequences in as strings
sequenceA = args.sequenceA.read().strip()
sequenceB = args.sequenceB.read().strip()

# Read in the score matrix
with open(BLOSUM, 'r') as f:
    content = f.readlines()
content = [line.strip() for line in content]
content = filter(lambda x: x[0] != "#", content)

# The first line of the score matrix holds labels
letters = content[0].split()
for iter in range(len(letters)):
    PROT_INDEX[letters[iter]] = iter
# Get rid of the wildcard "*"
del PROT_INDEX['*']

# Initialize the score matrix
BLOSUM = zeros((len(PROT_INDEX), len(PROT_INDEX)))
for iter in range(len(content) - 1):
    line = content[iter + 1]
    line = line.split()[1:(1 + len(PROT_INDEX))]
    BLOSUM[iter:] = [int(token) for token in line]

# Initialize the local alignment matrix
alignments = zeros((len(sequenceA), len(sequenceB)))

# Perform the alignment

