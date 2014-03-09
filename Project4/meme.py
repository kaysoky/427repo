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
            description='Takes a weight matrix model and a filtered SAM file (exported in JSON) and runs one iteration of the MEME algorithm on the WMM and data.  Outputs a new WMM')
    parser.add_argument('wmm', type=str, help='A weight matrix model')
    parser.add_argument('file', type=str, help='A filtered SAM file')
    parser.add_argument('output', type=str, help='Where to store the calculated weight matrix model')
    parser.add_argument('--verbose', action='store_true', help='Should progress and summary statistics be printed?')

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
    fileext = os.path.splitext(args.file)[1]
    iterator = None
    if fileext == JSON_FILE:
        iterator = json_generator(args.file)
    else:
        print 'Unknown input file format'
        exit()

    # Open the WMM output file
    output = open(args.output, 'w')

    # Use a generator to read in and process the file line by line
    for data in iterator:
        ## TODO
        ## Run MEME

    # Close the WMM output file
    output.close()

    if args.verbose:
        # Print how long the filter took
        elapsed = (time.clock() - startTime)
        print 'Analysis finished in %f seconds' % elapsed
