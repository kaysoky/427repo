import sys
import os
import argparse
import json
import time
import numpy

# Import some helper functions and globals
from shared import *

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
    assert_is_json_file(args.wmm)
    model = None
    with open(args.wmm, 'r') as f:
        model = json.load(f)
    model = normalize_wmm(model)

    # Handle SAM input
    assert_is_json_file(args.file)
    iterator = json_generator(args.file)

    # Open the WMM output file
    output = open(args.output, 'w')

    # Use a generator to read in and process the file line by line
    memeModel = numpy.zeros((4, WMM_LENGTH))
    for data in iterator:
        # Find the location of the tail
        tail = POLY_A_TAIL_SEARCH_REGEX.search(data[SAM_SEQ])
        tail = len(data[SAM_SEQ]) if tail is None else tail.start()

        # Extract the UTR region of the sequence
        sequence = matrixify_sequence(data[SAM_SEQ][:tail])
        
        # Apply the WMM to the sequence
        scores = apply_wmm_to_sequence(model, sequence)
        
        # Add up the weighted scores and add that to the new model
        aggregator = get_wmm_count_aggregator(probabilities=scores)
        modelDelta = numpy.dot(sequence, aggregator)
        modelDelta = normalize_wmm(modelDelta)
        memeModel += modelDelta

    # Close the WMM output file
    json.dump(memeModel.tolist(), output)
    output.close()

    if args.verbose:
        # Print how long the filter took
        elapsed = (time.clock() - startTime)
        print 'Iteration finished in %f seconds' % elapsed
