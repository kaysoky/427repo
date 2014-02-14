import sys
import os
import argparse
import re
import numpy
import json
import time

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
Emission probability matrix
"""
EMISSION_PROBABILITY = {}

"""
Initial belief of state
"""
INITIAL_STATE_PROBABILITY = {
    STATE_LOW_GC : 0.9999,
    STATE_HIGH_GC: 0.0001
}

"""
Transition probability matrix
"""
TRANSITION_PROBABILITY = {}

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

def run_viterbi(sequence):
    """
    Runs the Viterbi algorithm on the given sequence
        using the global probability matrices
    Returns a dictionary containing one entry per state
        Each entry is a list of 2-tuples denoting indices (inclusive, exclusive)
            which are predicted to be in the given state
    Also returns the log probability of the Viterbi path
    """

    # Initialize space to hold the calculated Viterbi probabilities
    max_log_probabilities = zeros(len(STATES))
    previous_state        = zeros((len(STATES), len(sequence)))
    #Note: The first column of 'previous_state' is meaningless

    # Fill in the start state probability
    max_log_probabilities[:] = [
        log(INITIAL_STATE_PROBABILITY[STATE_LOW_GC] * EMISSION_PROBABILITY[STATE_LOW_GC][sequence[0]]),
        log(INITIAL_STATE_PROBABILITY[STATE_HIGH_GC] * EMISSION_PROBABILITY[STATE_HIGH_GC][sequence[0]])
    ]

    # Define how to step through the Viterbi algorithm and fill in the appropriate values
    probabilities = zeros((len(STATES), len(STATES)));
    def viterbi_step(previous, index):
        # Calculate the len(STATES) ** 2 number of probabilities
        # Holds the probability[End state, Start state]
        for endState in STATES:
            for startState in STATES:
                probabilities[endState, startState] = \
                    previous[startState] \
                    + log(TRANSITION_PROBABILITY[startState][endState] \
                        * EMISSION_PROBABILITY[endState][sequence[index]])

        # Keep the maximums
        previous_state[:, index] = numpy.argmax(probabilities, 1)
        return numpy.max(probabilities, 1)

    # Fill in the remainder of the probabilities
    max_log_probabilities[:] = reduce(viterbi_step, range(1, len(sequence)), max_log_probabilities)

    # Back-trace to find all the predicted states
    # First initialize the space required
    viterbi_path = zeros(len(sequence))

    # Fill in the last value
    viterbi_path[len(sequence) - 1] = numpy.argmax(max_log_probabilities)

    # Define how to step back through the states and determine the Viterbi path
    def viterbi_backstep(previous, index):
        viterbi_path[index] = previous_state[previous, index + 1]
        return viterbi_path[index]

    # Back-fill in the remainder of the path
    reduce(viterbi_backstep, range(len(sequence) - 2, -1, -1), viterbi_path[len(sequence) - 1])

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

    return results, numpy.max(max_log_probabilities)

def run_transition_training(viterbi):
    """
    Takes the first output result of running the Viterbi algorithm
    And calculates a new transition probability matrix
        Note: this replaces the global probability matrices
    """

    # First zero-out the matrix
    for start in STATES:
        TRANSITION_PROBABILITY[start] = {}
        for end in STATES:
            TRANSITION_PROBABILITY[start][end] = 0

    # Next, let's unwrap the results
    #   from: {state: (start, end), ...}
    #     to: [(start, end, state), ...]
    unwrapped = []
    for state in STATES:
        unwrapped.extend([(item[0], item[1], state) for item in viterbi[state]])

    # Sort the unwrapped results in ascending order
    unwrapped.sort(key=lambda x: x[0])

    # Add up transition counts
    for index in range(len(unwrapped)):
        start = unwrapped[index][2]

        # Add a single transition from one state to another
        if index < len(unwrapped) - 1:
            end = unwrapped[index + 1][2]
            TRANSITION_PROBABILITY[start][end] += 1

        # Add transitions from one state to itself
        TRANSITION_PROBABILITY[start][start] += unwrapped[index][1] - unwrapped[index][0]

    # Normalize the transition counts to probabilities
    for start in STATES:
        total = 0.0
        for end in STATES:
            total += TRANSITION_PROBABILITY[start][end]
        for end in STATES:
            TRANSITION_PROBABILITY[start][end] /= total

def run_emission_training(sequence, viterbi):
    """
    Takes the first output result of running the Viterbi algorithm
        And the sequence used to calculate the output
    And calculates a new emission probability matrix
        Note: this replaces the global probability matrices
    """

    # First zero-out the matrix
    for start in STATES:
        EMISSION_PROBABILITY[start] = {}
        for emit in NUCLEOTIDES:
            EMISSION_PROBABILITY[start][emit] = 0

    # Add up the emission counts in each state
    for state in viterbi:
        # Concatenate each subsequence together
        subsequence = []
        for result in viterbi[state]:
            subsequence.append(sequence[result[0]:result[1]])
        subsequence = ''.join(subsequence)

        # Count the number of each emission
        for emit in NUCLEOTIDES:
            EMISSION_PROBABILITY[state][emit] += reduce(
                    lambda x, y: x + (1 if subsequence[y] == emit else 0),
                    range(len(subsequence)),
                    0) # Initializer, necessary otherwise the first letter will be ignored

    # Normalize the transition counts to probabilities
    for start in STATES:
        total = 0.0
        for emit in NUCLEOTIDES:
            total += EMISSION_PROBABILITY[start][emit]
        for emit in NUCLEOTIDES:
            EMISSION_PROBABILITY[start][emit] /= total

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Runs the Viterbi algorithm on the given sequence, looking for areas of high GC content')
    parser.add_argument('sequence', type=str)
    parser.add_argument('--inE', type=str, help='JSON matrix used to initialize the emission probability matrix')
    parser.add_argument('--inT', type=str, help='JSON matrix used to initialize the transition probability matrix')
    parser.add_argument('--outE', type=str, required=False,
        help='File to save the trained emission probability matrix')
    parser.add_argument('--outT', type=str, required=False,
        help='File to save the trained transition probability matrix')
    parser.add_argument('--verbose', action='store_true', help='Print every high GC hit?')
    parser.add_argument('--time', action='store_true', help='Time the algorithm?')
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

    # Load the the emission probability matrix
    with open(args.inE) as f:
        emission = json.load(f)

        # JSON doesn't support integer keys
        # But the states are defined as integers
        # So convert it here
        for key in emission.keys():
            EMISSION_PROBABILITY[int(key)] = emission[key]

    # Load the the transition probability matrix
    with open(args.inT) as f:
        transition = json.load(f)

        # Again, JSON doesn't support integer keys
        for key in transition.keys():
            subTransition = {}
            for subkey in transition[key]:
                subTransition[int(subkey)] = transition[key][subkey]
            TRANSITION_PROBABILITY[int(key)] = subTransition

    # If specified, time the iteration
    if args.time:
        startTime = time.clock()

    # Run the Viterbi algorithm
    results, max_prob = run_viterbi(sequence)

    # Print the results
    print 'Viterbi log probability: %f' % max_prob
    print 'Number of high GC content hits: %d' % len(results[STATE_HIGH_GC])
    if args.verbose:
        print 'High GC content hits (one-based index):'
        print '     Start | End '
        for item in results[STATE_HIGH_GC]:
            print '%*d | %d' % (10, item[0] + 1, item[1] + 1)

    # Train the new emission matrix
    if args.outE is not None:
        run_emission_training(sequence, results)
        with open(args.outE, 'w') as f:
            json.dump(EMISSION_PROBABILITY, f, indent=2)

    # Train the new transition matrix
    if args.outT is not None:
        run_transition_training(results)
        with open(args.outT, 'w') as f:
            json.dump(TRANSITION_PROBABILITY, f, indent=2)

    # If specified, print the time elapsed
    if args.time:
        elapsed = (time.clock() - startTime)
        print 'Iteration time %f seconds' % elapsed

