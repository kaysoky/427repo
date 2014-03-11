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
    parser.add_argument('background', type=str, help='A wight matrix model of the background')
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
    assert_is_json_file(args.wmm)
    model = None
    with open(args.wmm, 'r') as f:
        model = json.load(f)
    model = normalize_wmm(model)

    # Open the background model
    assert_is_json_file(args.background)
    background = None
    with open(args.background, 'r') as f:
        background = json.load(f)
    background = normalize_wmm(background)

    # Handle SAM input
    assert_is_json_file(args.file)
    iterator = json_generator(args.file)

    # Open the histogram file
    output = open(args.output, 'w')
    
    # Declare the statistics we're looking for
    motifHits = 0
    motifDistance = 0
    motifHistogram = numpy.zeros(101)

    # Use a generator to read in and process the file line by line
    counter = 0
    for data in iterator:
        # Limit the number of input lines read (for debugging purposes)
        if args.limit is not None and counter >= args.limit:
            break
        counter += 1

        # Find the location of the tail
        tail = POLY_A_TAIL_SEARCH_REGEX.search(data[SAM_SEQ])
        tail = len(data[SAM_SEQ]) if tail is None else tail.start()
        
        # Extract the UTR region of the sequence
        sequence = matrixify_sequence(data[SAM_SEQ][:tail])
        
        # Apply the WMM and background to the sequence
        scores = apply_wmm_to_sequence(model, sequence) / apply_wmm_to_sequence(background, sequence)
        
        # Check for a motif hit
        scores = numpy.log(scores)
        maxScore = numpy.max(scores)
        if maxScore <= 0:
            continue
            
        # There's a hit, so tabulate it
        lastHitIndex = numpy.where(scores == maxScore)[-1][-1]
        lastHitDistance = tail - lastHitIndex
        motifHits += 1
        motifDistance += lastHitDistance
        motifHistogram[lastHitDistance] += 1
        
    ##TODO
    # Output the number of hits and average distance
    # Output the histogram

    # Close the histogram output file
    output.close()

    if args.verbose:
        # Print how long the filter took
        elapsed = (time.clock() - startTime)
        print 'Analysis finished in %f seconds' % elapsed
