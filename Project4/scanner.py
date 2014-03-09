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
            description='Takes a weight matrix model and a filtered SAM file (exported in JSON) and applies the WMM on the data.  Outputs the number of model "hits", average distance from hit to cleave site, and a histogram of hit positions')
    parser.add_argument('wmm', type=str, help='A weight matrix model')
    parser.add_argument('file', type=str, help='A filtered SAM file')
    parser.add_argument('output', type=str, help='File to store the histogram table')
    parser.add_argument('--limit', type=int, required=False, help='How many lines of input should be read?')
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

    # Open the histogram file
    output = open(args.output, 'w')

    # Use a generator to read in and process the file line by line
    counter = 0
    for data in iterator:
        # Limit the number of input lines read (for debugging purposes)
        if args.limit is not None and counter >= args.limit:
            break
        counter += 1

        ## TODO
        ## Compute number of motif hits
        ##   average distance to cleave site
        ##   histogram of hit positions

    # Close the histogram output file
    output.close()

    if args.verbose:
        # Print how long the filter took
        elapsed = (time.clock() - startTime)
        print 'Analysis finished in %f seconds' % elapsed
