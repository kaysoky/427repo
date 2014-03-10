import sys
import os
import argparse
import json
import numpy

# Import some helper functions and globals
from shared import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Takes a weight matrix model and a background model and computes the relative entropy')
    parser.add_argument('wmm', type=str, help='A weight matrix model')
    parser.add_argument('background', type=str, help='A wight matrix model of the background')

    args = parser.parse_args()

    # Open the WMM
    fileext = os.path.splitext(args.wmm)[1]
    model = None
    if fileext == JSON_FILE:
        with open(args.wmm, 'r') as f:
            model = json.load(f)
    else:
        print 'Unknown input file format'
        exit()

    # Open the background model
    fileext = os.path.splitext(args.background)[1]
    background = None
    if fileext == JSON_FILE:
        with open(args.background, 'r') as f:
            background = json.load(f)
    else:
        print 'Unknown input file format'
        exit()

    # Normalize the models
    model = normalize_wmm(model)
    background = normalize_wmm(background)

    # Compute the relative entropy
    entropy = numpy.sum(model * numpy.log(model / background) / numpy.log(2))
    
    print "Relative entropy: %f" % entropy
