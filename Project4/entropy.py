import sys
import os
import argparse
import json
import time
import numpy

# Import some helper functions and globals
from filter import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Takes a weight matrix model and a background model and computes the relative entropy')
    parser.add_argument('wmm', type=str, help='A weight matrix model')
    parser.add_argument('background', type=str, help='A filtered SAM file')

    args = parser.parse_args()

    # Keep track of how much time this filter requires
    startTime = 0
    if args.verbose:
        startTime = time.clock()

    # Open the WMM
    fileext = os.path.splitext(args.wmm)[1]
    model = None
    if fileext == JSON_FILE:
        with open(args.wmm, 'r') as f:
            model = numpy.array(json.load(f))
    else:
        print 'Unknown input file format'
        exit()

    # Handle SAM input
    fileext = os.path.splitext(args.background)[1]
    background = None
    if fileext == JSON_FILE:
        with open(args.background, 'r') as f:
            background = numpy.array(json.load(f))
    else:
        print 'Unknown input file format'
        exit()

    ## TODO
    ## Run MEME
