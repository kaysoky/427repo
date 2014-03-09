import sys
import os
import argparse
import re
import json
import time
import numpy

"""
File extension indicating that an input file should be parsed as a SAM file
"""
SAM_FILE = '.sam'

"""
File extension indicating that an input file should be parsed as a JSON file
"""
JSON_FILE = '.json'

"""
Expected length of sequences within the SAM file
"""
SEQUENCE_LENGTH = 75

"""
Motif length the Weight Matrix Model will represent
"""
WMM_LENGTH = 6

"""
IUPAC nucleotide base notation
"""
NUCLEOTIDES = {
    'A': 0, 
    'C': 1, 
    'G': 2, 
    'T': 3,
    'R': 4, # A or G
    'Y': 5, # C or T
    'S': 6, # G or C
    'W': 7, # A or T
    'K': 8, # G or T
    'M': 9, # A or C
    'B': 10, # C or G or T
    'D': 11, # A or G or T
    'H': 12, # A or C or T
    'V': 13, # A or C or G
    'N': 14} # any base
    
"""
Used to transform a 15x75 matrix of bases type x location in sequence
  To a 15x6 matrix of bases in a WMM
"""
WMM_BASE_COUNT_AGGREGATOR = numpy.ones((SEQUENCE_LENGTH, WMM_LENGTH))
WMM_BASE_COUNT_AGGREGATOR[0:WMM_LENGTH, 0:WMM_LENGTH] += numpy.tril(numpy.ones((WMM_LENGTH, WMM_LENGTH))) - 1
WMM_BASE_COUNT_AGGREGATOR[(SEQUENCE_LENGTH - WMM_LENGTH):, 0:WMM_LENGTH] += numpy.triu(numpy.ones((WMM_LENGTH, WMM_LENGTH))) - 1
    
"""
Used to reduce the number of bases from 15 to 4 in a 15x6 WMM
"""
WMM_STANDARD_COUNT = numpy.array([
    [1, 0, 0, 0, 0.5, 0  , 0  , 0.5, 0  , 0.5, 0   , 0.33, 0.33, 0.33, 0.25],
    [0, 1, 0, 0, 0  , 0.5, 0.5, 0  , 0  , 0.5, 0.33, 0   , 0.33, 0.33, 0.25],
    [0, 0, 1, 0, 0.5, 0  , 0.5, 0  , 0.5, 0  , 0.33, 0.33, 0   , 0.33, 0.25],
    [0, 0, 0, 1, 0  , 0.5, 0  , 0.5, 0.5, 0  , 0.33, 0.33, 0.33, 0   , 0.25]
])

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

# Optional fields
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
Custom bit mask for the SAM_FLAG parameter
Indicates if the sequence was a reverse complement 
  but has been transformed by this filter into a 'ordinary' sequence
"""
SAM_TRANSFORMED_REV_COMP_FLAG_MASK = 0x1000

"""
Usage: POLY_A_TAIL_SEARCH_REGEX.search(data[SAM_SEQ])
Finds the poly-A tail of the given sequence
Note: This is extremely lenient and counts uncertain A's as part of the tail
"""
POLY_A_TAIL_SEARCH_REGEX = re.compile(r"([ARWMDHVN]+)$")
POLY_T_TAIL_SEARCH_REGEX = re.compile(r"^([TYWKBDHN]+)")

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
    And iterates over the values and returns the decoded JSON
    """

    with open(filename, 'r') as file:
        for line in file:
            yield json.loads(line)
            
def complement_sequence(sequence):
    """
    Takes a nucleotide sequence and complements it
    """
    
    # Replace nucleotides with a lowercase nucleotide (to prevent re-replacement)
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('R', 'y')
    sequence = sequence.replace('Y', 'r')
    sequence = sequence.replace('K', 'm')
    sequence = sequence.replace('M', 'k')
    sequence = sequence.replace('B', 'v')
    sequence = sequence.replace('D', 'h')
    sequence = sequence.replace('H', 'd')
    sequence = sequence.replace('V', 'b')
    
    # Uppercase the replaced string
    return sequence.upper()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Opens a SAM file or a processed SAM file (exported in JSON) and performs the specified set of filtering operations')
    parser.add_argument('file', type=str, help='A SAM or previously filtered SAM file')
    parser.add_argument('output', type=str, help='Output JSON file.  Note: every line of the file will contain a JSON string')
    parser.add_argument('--limit', type=int, required=False, help='How many lines of input should be read?')
    parser.add_argument('--verbose', action='store_true', help='Should progress and summary statistics be printed?')
    parser.add_argument('--dereverse', action='store_true', help='Should any reverse complements be reversed and complemented into "ordinary" sequences?')

    # Filtering parameters
    parser.add_argument('--matches_only', action='store_true', help='Filters out data that does not match map to the reference genome')
    parser.add_argument('--min_mismatch', type=int, help='Filters out data with less than (exclusive) the given number of mismatches')
    parser.add_argument('--max_align_score', type=int, help='Filters out data with an alignment score greater than (exclusive) the given value; Note: scores are negative')
    parser.add_argument('--min_polyAlen', type=int, help='Filters out data with a trailing poly-A tail of less than (exclusive) the given length')
    parser.add_argument('--min_UTRlen', type=int, help='Filters out data with a leading 3\' UTR region of less than (exclusive) the given length')
    parser.add_argument('--max_non_tail_mismatches', type=int, help='Filters out data which contains more than (exclusive) the given number of mismatches in the non-tail region')
    parser.add_argument('--compute_background', type=str, help='For all data not passing the filter, adds the data to a 6mer weight matrix model of the background')

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

    # Since the input file might not fit in memory
    # Use a generator to read in and process the file line by line
    counter = 0
    outputLines = 0
    backgroundModel = numpy.zeros((4, WMM_LENGTH))
    for data in iterator:
        # Limit the number of input lines read (for debugging purposes)
        if args.limit is not None and counter >= args.limit:
            break
        counter += 1
        
        # Transform reverse complements into a non-reverse complement
        if args.dereverse and int(data[SAM_FLAG]) & SAM_REV_COMPLEMENT_FLAG_MASK:
            # Reverse and complement the relevant fields
            data[SAM_SEQ] = complement_sequence(data[SAM_SEQ][::-1])
            data[SAM_QUAL] = data[SAM_QUAL][::-1]
            if SAM_MSMAT in data:
                mismatches = ''.join(MISMATCH_SEARCH_REGEX.findall(data[SAM_MSMAT])[::-1])
                data[SAM_MSMAT] = complement_sequence(mismatches)
            
            # Flip the relevant flag bits
            data[SAM_FLAG] = int(data[SAM_FLAG]) ^ (SAM_REV_COMPLEMENT_FLAG_MASK | SAM_TRANSFORMED_REV_COMP_FLAG_MASK)
        
        backgroundDelta = None
        if args.compute_background:
            # Flatten the sequence into a 15x75 matrix
            dataSum = numpy.zeros((len(NUCLEOTIDES), SEQUENCE_LENGTH))
            for index in range(SEQUENCE_LENGTH):
                base = data[SAM_SEQ][index]
                dataSum[NUCLEOTIDES[base], index] = 1
                
            # Aggregate each base into the appropriate bucket in the Weight Matrix
            backgroundDelta = numpy.dot(dataSum, WMM_BASE_COUNT_AGGREGATOR)
            
            # Remove non ACGT bases
            backgroundDelta = numpy.dot(WMM_STANDARD_COUNT, backgroundDelta)
            
            # Save the delta
            backgroundModel += backgroundDelta

        # Remove all sequences that are unmapped
        if args.matches_only:
            if int(data[SAM_FLAG]) & SAM_UNMAPPED_FLAG_MASK:
                continue

        # Remove all sequences that strongly match the reference
        if args.min_mismatch and SAM_NUMMM in data:
            if int(data[SAM_NUMMM]) < args.min_mismatch:
                continue

        # Filter out sequences with high scores
        if args.max_align_score and SAM_A_SCR in data:
            if int(data[SAM_A_SCR]) > args.max_align_score:
                continue

        # Find the poly-A tail region if necessary
        tail = None
        if args.min_polyAlen or args.max_non_tail_mismatches:
            regex = POLY_A_TAIL_SEARCH_REGEX

            # Account for reverse complements
            if int(data[SAM_FLAG]) & SAM_REV_COMPLEMENT_FLAG_MASK:
                regex = POLY_T_TAIL_SEARCH_REGEX

            tail = regex.search(data[SAM_SEQ])
        if tail is not None:
            tail = tail.groups()[0]

        # Remove all sequences without a significant poly-A tail
        if args.min_polyAlen:
            if tail is None or len(tail) < args.min_polyAlen:
                continue
                
        # Remove all sequences with a short UTR region
        if args.min_UTRlen and tail is not None:
            if len(data[SAM_SEQ]) - len(tail) < args.min_UTRlen:
                continue

        # Remove all sequences with major mismatching in the 3' UTR
        # Note: to prevent too much code duplication, 
        #         this filter only runs on "ordinary" sequences
        if args.max_non_tail_mismatches and SAM_MSMAT in data \
                and int(data[SAM_FLAG]) & SAM_REV_COMPLEMENT_FLAG_MASK:
            tailIndex = len(data[SAM_SEQ]) if tail is None else tail.start()

            # Count mismatches up until the tail index is hit
            misCount = 0
            misIndex = 0
            mismatches = MISMATCH_SEARCH_REGEX.findall(data[SAM_MSMAT])
            for token in mismatches:
                if misIndex >= tailIndex:
                    break
                    
                # Token indicates a number of matches
                try:
                    misIndex += int(token)
                    continue
                except exceptions.ValueError:
                    pass
                    
                # Token indicates a deletion
                if len(token) > 1:
                    misCount += len(token) - 1
                    misIndex += len(token) - 1
                    
                # Token indicates a mismatch
                else:
                    misCount += 1
                    misIndex += 1

            if misCount > args.max_non_tail_mismatches:
                continue

        # Data passed the filter, so save it
        output.write(json.dumps(data))
        output.write('\n')
        outputLines += 1
        
        # Data passed the filter, so exclude it from the background
        if args.compute_background:
            backgroundModel -= backgroundDelta

    # Close the JSON output file
    output.close()
    
    # Save the background model
    if args.compute_background:
        with open(args.compute_background, 'w') as f:
            json.dump(backgroundModel.tolist(), f)
    
    if args.verbose:
        # Print how long the filter took
        elapsed = (time.clock() - startTime)
        print 'Filter finished in %f seconds' % elapsed
        
        # Print summary statistics about the result of the filter
        print '%d lines of input -> %d lines of output' % (counter, outputLines)
