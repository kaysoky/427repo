import sys
import os
import argparse
import re
import json
import time

"""
File extension indicating that an input file should be parsed as a SAM file
"""
SAM_FILE = '.sam'

"""
File extension indicating that an input file should be parsed as a JSON file
"""
JSON_FILE = '.json'

#####################
## SAM File Fields ##
#####################
# Note: Refer to section 1.4 in the SAM file specification
SAM_QNAME = 'QNAME'
SAM_FLAG  = 'FLAG'
SAM_RNAME = 'RNAME'
SAM_POS   = 'POS'
SAM_MAPQ  = 'MAPQ'
SAM_CIGAR = 'CIGAR'
SAM_RNEXT = 'RNEXT'
SAM_PNEXT = 'PNEXT'
SAM_TLEN  = 'TLEN'
SAM_SEQ   = 'SEQ'
SAM_QUAL  = 'QUAL'
SAM_A_SCR = 'AS:i:'
SAM_NUMMM = 'NM:i:'
SAM_MSMAT = 'MD:Z:'

"""
Bit mask for the SAM_FLAG parameter
Indicates if the sequence is unmapped to the reference
"""
SAM_UNMAPPED_FLAG_MASK = 0x4

"""
Bit mask for the SAM_FLAG parameter
Indicates if the sequence is a reverse complement
"""
SAM_REV_COMPLEMENT_FLAG_MASK = 0x10

"""
Usage: POLY_A_TAIL_SEARCH_REGEX.search(data[SAM_SEQ])
Finds the poly-A tail of the given sequence
"""
POLY_A_TAIL_SEARCH_REGEX = re.compile(r"(A+)$")
POLY_T_TAIL_SEARCH_REGEX = re.compile(r"(T+)$")

"""
Usage: MISMATCH_SEARCH_REGEX.findall(data[SAM_MSMAT])
Tokenizes the string representing mismatching positions of a sequence mapping
"""
MISMATCH_SEARCH_REGEX = re.compile(r"([A-Z]|\^[A-Z]+|[0-9]+)")

def parse_SAM_data(line):
    """
    Takes a tab-delimited line of processed mapping data in SAM format
    And parses it into a dictionary of appropriate values
    """

    result = {}
    tokens = line.strip().split('\t')

    # Capture the 11 mandatory fields in order
    result[SAM_QNAME] = tokens[ 0]
    result[SAM_FLAG ] = tokens[ 1]
    result[SAM_RNAME] = tokens[ 2]
    result[SAM_POS  ] = tokens[ 3]
    result[SAM_MAPQ ] = tokens[ 4]
    result[SAM_CIGAR] = tokens[ 5]
    result[SAM_RNEXT] = tokens[ 6]
    result[SAM_PNEXT] = tokens[ 7]
    result[SAM_TLEN ] = tokens[ 8]
    result[SAM_SEQ  ] = tokens[ 9]
    result[SAM_QUAL ] = tokens[10]

    # Parse the remaining interesting optional fields
    for token in tokens[11:]:
        if token.startswith(SAM_A_SCR):
            result[SAM_A_SCR] = token[len(SAM_A_SCR):]

        elif token.startswith(SAM_NUMMM):
            result[SAM_NUMMM] = token[len(SAM_NUMMM):]

        elif token.startswith(SAM_MSMAT):
            result[SAM_MSMAT] = token[len(SAM_MSMAT):]

    return result

def sam_generator(filename):
    """
    Opens the given SAM file and iterates over the file
    Each iteration parses and returns a single sequence with its mapping data
    """

    with open(filename, 'r') as file:
        for line in file:
            # Skip headers
            if line.startswith('@'):
                continue

            # Parse sequence mapping data
            else:
                yield parse_SAM_data(line)

def json_generator(filename):
    """
    Opens a processed SAM file that was exported in JSON
    And iterates over the array, returning each value
    """

    with open(filename, 'r') as file:
        data = json.load(file)
        for item in data:
            if item: # Not empty
                yield item

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Opens a SAM file or a processed SAM file (exported in JSON) and performs the specified set of filtering operations')
    parser.add_argument('file', type=str, help='A SAM or processed SAM file')
    parser.add_argument('output', type=str, help='Output JSON file')
    parser.add_argument('--limit', type=int, required=False, help='How many lines of input should be read?')
    parser.add_argument('--verbose', action='store_true', help='Should progress and summary statistics be printed?')

    # Filtering parameters
    parser.add_argument('--matches_only', action='store_true', help='Filters out data that does not match map to the reference genome')
    parser.add_argument('--min_mismatch', type=int, help='Filters out data with less than (exclusive) the given number of mismatches')
    parser.add_argument('--max_align_score', type=int, help='Filters out data with an alignment score greater than (exclusive) the given value; Note: scores are negative')
    parser.add_argument('--min_polyAlen', type=int, help='Filters out data with a trailing poly-A tail of less than (exclusive) the given length')
    parser.add_argument('--max_non_tail_mismatches', type=int, help='Filters out data which contains more than (exclusive) the given number of mismatches in the non-tail region')

    args = parser.parse_args()
    
    # Keep track of how much time this filter requires
    startTime = 0
    if args.verbose:
        startTime = time.clock()

    # Handle input parsing
    fileext = os.path.splitext(args.file)[1]
    iterator = None
    if fileext == SAM_FILE:
        iterator = sam_generator(args.file)
    elif fileext == JSON_FILE:
        iterator = json_generator(args.file)
    else:
        print 'Unknown input file format'
        exit()

    # Open up the output file for writing
    fileext = os.path.splitext(args.output)[1]
    if fileext != JSON_FILE:
        print 'Only JSON file format is acceptable for output'
        exit()
    output = open(args.output, 'w')
    output.write('[\n')

    # Since the input file might not fit in memory
    # Use a generator to read in and process the file line by line
    counter = 0
    outputLines = 0
    for data in iterator:
        # Limit the number of input lines read (for debugging purposes)
        if args.limit is not None and counter >= args.limit:
            break
        counter += 1

        # Remove all sequences that are unmapped
        if args.matches_only:
            if int(data[SAM_FLAG]) & SAM_UNMAPPED_FLAG_MASK:
                continue

        # Remove all sequences that strongly match the reference
        if args.min_mismatch:
            if SAM_NUMMM in data and int(data[SAM_NUMMM]) < args.min_mismatch:
                continue

        # Filter out sequences with high scores
        if args.max_align_score:
            if SAM_A_SCR in data and int(data[SAM_A_SCR]) > args.max_align_score:
                continue

        # Find the poly-A tail region if necessary
        tail = None
        if args.min_polyAlen or args.max_non_tail_mismatches:
            regex = POLY_A_TAIL_SEARCH_REGEX

            # Account for reverse complements
            if int(data[SAM_FLAG]) & SAM_REV_COMPLEMENT_FLAG_MASK:
                regex = POLY_T_TAIL_SEARCH_REGEX

            tail = regex.search(data[SAM_SEQ])

        # Remove all sequences without a significant poly-A tail
        if args.min_polyAlen:
            if tail is None or len(tail.groups()[0]) < args.min_polyAlen:
                continue

        # Remove all sequences with major mismatching in the 3' UTR
        if args.max_non_tail_mismatches:
            tailIndex = len(data[SAM_SEQ]) if tail is None else tail.start()

            # Count mismatches up until the tail index is hit
            misCount = 0
            misIndex = 0
            mismatches = MISMATCH_SEARCH_REGEX.findall(data[SAM_MSMAT])
            for token in mismatches:
                if misIndex >= tailIndex:
                    break
                    
                ##TODO

            if misCount > args.max_non_tail_mismatches:
                continue

        # Data passed the filter, so save it
        output.write(json.dumps(data))
        output.write(',\n')
        outputLines += 1

    # Close the JSON array with an empty item
    output.write('{}]')
    output.close()
    
    if args.verbose:
        # Print how long the filter took
        elapsed = (time.clock() - startTime)
        print 'Filter finished in %f seconds' % elapsed
        
        # Print summary statistics about the result of the filter
        print '%d lines of input -> %d lines of output' % (counter, outputLines)
