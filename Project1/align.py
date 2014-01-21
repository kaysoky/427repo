import argparse

import numpy
from numpy import zeros

parser = argparse.ArgumentParser(
        description='Performs BLAST (Smith-Waterman) on sequences found in two files')
# parser.add_argument('scoreMatrix', type=file)
parser.add_argument('sequenceA', type=file)
parser.add_argument('sequenceB', type=file)
args = parser.parse_args()

# Read the two sequences in as strings
sequenceA = args.sequenceA.read()
sequenceB = args.sequenceB.read()

# TODO
