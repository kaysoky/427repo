import sys
import os
import argparse
import re
import numpy

from math import log
from numpy import zeros

"""
File extension that triggers some additional processing
Lines starting with '>' are removed
    and whitespace is stripped out
"""
FASTA = '.fna'

"""
The four possible nucleotides
Also the four emissions of the Hidden Markov Model
"""
NUCLEOTIDES = ['A', 'C', 'G', 'T']

#######################################
## States of the Hidden Markov Model ##
#######################################
# Note: The states should be consecutive numbers, in order, starting from 0
STATE_LOW_GC = 0
STATE_HIGH_GC = 1

"""
Convenient collection of all states in the model
"""
STATES = [STATE_LOW_GC, STATE_HIGH_GC]

"""
Initial emission probability matrix
This will be updated by any iterative Viterbi training
"""
EMISSION_PROBABILITY = {
    STATE_LOW_GC : {
        'A' : 0.25,
        'C' : 0.25,
        'G' : 0.25,
        'T' : 0.25
    },
    STATE_HIGH_GC : {
        'A' : 0.20,
        'C' : 0.30,
        'G' : 0.30,
        'T' : 0.20
    }
}

"""
Initial belief of state
This remains constant
"""
INITIAL_STATE_PROBABILITY = {
    STATE_LOW_GC : 0.9999,
    STATE_HIGH_GC: 0.0001
}

"""
Initial transition probability matrix
This will be updated by any iterative Viterbi training
"""
TRANSITION_PROBABILITY = {
    STATE_LOW_GC : {
        STATE_LOW_GC : 0.9999,
        STATE_HIGH_GC: 0.0001
    },
    STATE_HIGH_GC : {
        STATE_LOW_GC : 0.01,
        STATE_HIGH_GC: 0.99
    }
}

def process_fasta(text):
    """
    Removes the comments and whitespace
    Also case-desensitizes the input and
        replaces all non-ACGT characters with 'T'
    """

    lines = text.split('\n')
    lines = filter(lambda x: not x.startswith('>'), lines)
    lines = [line.strip() for line in lines]
    text = ''.join(lines)
    text = text.upper()
    return re.sub(r'[^ACGT]', 'T', text)

def run_viterbi(sequence, verbose=False):
    """
    Runs the Viterbi algorithm on the given sequence
        using the global probability matrices
    Returns a dictionary containing one entry per state
        Each entry is a list of 2-tuples denoting indices (inclusive, exclusive)
            which are predicted to be in the given state
    Also returns the log probability of the Viterbi path
    """

    # Initialize space to hold the calculated Viterbi probabilities
    max_log_probabilities = zeros((len(STATES), len(sequence)))
    previous_state        = zeros((len(STATES), len(sequence)))
    #Note: The first column of 'previous_state' is meaningless

    # Fill in the start state probability
    max_log_probabilities[:, 0] = [
        log(INITIAL_STATE_PROBABILITY[STATE_LOW_GC] * EMISSION_PROBABILITY[STATE_LOW_GC][sequence[0]]),
        log(INITIAL_STATE_PROBABILITY[STATE_HIGH_GC] * EMISSION_PROBABILITY[STATE_HIGH_GC][sequence[0]])
    ]

    # Fill in the remainder of the probabilities
    probabilities = zeros((len(STATES), len(STATES)));
    for index in range(1, len(sequence)):
        # Calculate the len(STATES) ** 2 number of probabilities
        # Holds the probability[End state, Start state]
        probabilities[:,:] = [
            [
                max_log_probabilities[STATE_LOW_GC, index - 1]
                    + log(TRANSITION_PROBABILITY[STATE_LOW_GC][STATE_LOW_GC]
                        * EMISSION_PROBABILITY[STATE_LOW_GC][sequence[index]]),
                max_log_probabilities[STATE_HIGH_GC, index - 1]
                    + log(TRANSITION_PROBABILITY[STATE_HIGH_GC][STATE_LOW_GC]
                        * EMISSION_PROBABILITY[STATE_LOW_GC][sequence[index]]),
            ], [
                max_log_probabilities[STATE_LOW_GC, index - 1]
                    + log(TRANSITION_PROBABILITY[STATE_LOW_GC][STATE_HIGH_GC]
                        * EMISSION_PROBABILITY[STATE_HIGH_GC][sequence[index]]), 
                max_log_probabilities[STATE_HIGH_GC, index - 1]
                    + log(TRANSITION_PROBABILITY[STATE_HIGH_GC][STATE_HIGH_GC]
                        * EMISSION_PROBABILITY[STATE_HIGH_GC][sequence[index]])
            ]
        ]

        # Keep the maximums
        max_log_probabilities[:, index] = numpy.max(probabilities, 1)
        previous_state[:, index] = numpy.argmax(probabilities, 1)
        if verbose and index % 1000 == 0:
            print index, sequence[index], previous_state[:, index], max_log_probabilities[:, index]

    # Back-trace to find all the predicted states
    # First initialize the space required
    viterbi_path = zeros(len(sequence))

    # Fill in the last value
    viterbi_path[len(sequence) - 1] = numpy.argmax(max_log_probabilities[:, len(sequence) - 1])

    # Back-fill in the remainder of the path
    for index in range(len(sequence) - 2, -1, -1):
        viterbi_path[index] = previous_state[viterbi_path[index + 1], index + 1]

    # Yay, now I can use functional programming to do the remaining transformation
    # Namely, the transformation of the path to the return value
    # First, find all the indices where there is a state transition 
    #   This will be the inclusive end of all the contiguous state chunks
    thresholds = filter(lambda x: viterbi_path[x] != viterbi_path[x + 1], 
                        range(0, len(viterbi_path) - 1))
    # Pretend that there's also a state transition at the end
    thresholds.append(len(viterbi_path) - 1)
    
    # Now transform these indices into tuples of (inclusive start, exclusive end)
    thresholds = map(lambda x:(0 if x <= 0 else (thresholds[x - 1] + 1), thresholds[x] + 1), 
                     range(len(thresholds)))
     
    # And filter them into the proper buckets
    results = {}
    for state in STATES:
        results[state] = filter(lambda x: viterbi_path[x[0]] == state, thresholds)
        
    return results, numpy.max(max_log_probabilities[:, len(sequence) - 1])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='TODO')
    parser.add_argument('sequence', type=str)
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    # Read the sequence in as a string
    with open(args.sequence) as f:
        sequence = f.read().strip()

    # Handle sequence files
    if os.path.splitext(args.sequence)[1] == FASTA:
        sequence = process_fasta(sequence)
    else:
        print 'Unknown sequence file format'
        exit()

    results, max_prob = run_viterbi(sequence, args.verbose)
    print 'Viterbi Probability: %f' % max_prob
    for state in results.keys():
        print "-----"
        print state
        print "-----"
        for item in results[state]:
            print item
