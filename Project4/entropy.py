import sys
import os
import argparse
import json
import numpy

# Import some helper functions and globals
from shared import *

"""
After taking the logarithm of the model / background, 
  replace any infinities or NAN's with this value

"""
PSUEDOCOUNT = -100

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Takes a weight matrix model and a background model and computes the relative entropy')
    parser.add_argument('wmm', type=str, help='A weight matrix model')
    parser.add_argument('background', type=str, help='A wight matrix model of the background')

    args = parser.parse_args()

    # Open the WMM
    assert_is_json_file(args.wmm)
    model = None
    with open(args.wmm, 'r') as f:
        model = json.load(f)

    # Open the background model
    assert_is_json_file(args.background)
    background = None
    with open(args.background, 'r') as f:
        background = json.load(f)

    # Normalize the models
    model = normalize_wmm(model)
    background = normalize_wmm(background)

    # Compute the relative entropy
    entropy = numpy.log(model / background)
    entropy = numpy.where(numpy.isinf(entropy), PSUEDOCOUNT, entropy)
    entropy = numpy.where(numpy.isnan(entropy), PSUEDOCOUNT, entropy)
    entropy = numpy.sum(model * entropy / numpy.log(2))
    
    print "Relative entropy: %f" % entropy
