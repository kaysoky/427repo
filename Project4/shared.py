import os 
import re
import json
import numpy

"""
File extension indicating that an input file should be parsed as a JSON file
"""
JSON_FILE = '.json'

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

def assert_is_json_file(filename):
    fileext = os.path.splitext(filename)[1]
    assert fileext == JSON_FILE, 'File "%s" requires %s extension' % (filename, JSON_FILE)

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

def matrixify_sequence(sequence):
    """
    Converts a sequence into a 4 by N matrix of nucleotides to sequence position
    """

    # Flatten the sequence into a 15 by N matrix
    flattened = numpy.zeros((len(NUCLEOTIDES), len(sequence)))
    for index in range(len(sequence)):
        base = sequence[index]
        flattened[NUCLEOTIDES[base], index] = 1

    # Convert the 15 IUPAC nucleotides into the 4 canonical ones
    return numpy.dot(WMM_STANDARD_COUNT, flattened)

COUNT_AGGREGATOR_CACHE = {}
def get_wmm_count_aggregator(seqLen):
    """
    Returns a seqLen by WMM_LENGTH matrix used to aggregate nucleotide counts into a WMM
    """

    if seqLen in COUNT_AGGREGATOR_CACHE:
        return COUNT_AGGREGATOR_CACHE[seqLen]

    '''
    numpy.triu() and numpy.tril() return the upper and lower triangular portions of a matrix
    I use this to construct the middle matrix, which is subtracted from the left matrix:
        | 1 1 1 1 1 1 |   | 0 1 1 1 1 1 |   | 1 0 0 0 0 0 |
        | 1 1 1 1 1 1 |   | 0 0 1 1 1 1 |   | 1 1 0 0 0 0 |
        | 1 1 1 1 1 1 |   | 0 0 0 1 1 1 |   | 1 1 1 0 0 0 |
        | 1 1 1 1 1 1 |   | 0 0 0 0 1 1 |   | 1 1 1 1 0 0 |
        | 1 1 1 1 1 1 |   | 0 0 0 0 0 1 |   | 1 1 1 1 1 0 |
        | 1 1 1 1 1 1 |   | 0 0 0 0 0 0 |   | 1 1 1 1 1 1 |
        |     ...     | - |     ...     | = |     ...     |
        | 1 1 1 1 1 1 |   | 0 0 0 0 0 0 |   | 1 1 1 1 1 1 |
        | 1 1 1 1 1 1 |   | 1 0 0 0 0 0 |   | 0 1 1 1 1 1 |
        | 1 1 1 1 1 1 |   | 1 1 0 0 0 0 |   | 0 0 1 1 1 1 |
        | 1 1 1 1 1 1 |   | 1 1 1 0 0 0 |   | 0 0 0 1 1 1 |
        | 1 1 1 1 1 1 |   | 1 1 1 1 0 0 |   | 0 0 0 0 1 1 |
        | 1 1 1 1 1 1 |   | 1 1 1 1 1 0 |   | 0 0 0 0 0 1 |
    '''
    square = numpy.ones((WMM_LENGTH, WMM_LENGTH))
    aggregator = numpy.ones((seqLen, WMM_LENGTH))
    aggregator[0:WMM_LENGTH, 0:WMM_LENGTH] += numpy.tril(square) - 1
    aggregator[(seqLen - WMM_LENGTH):, 0:WMM_LENGTH] += numpy.triu(square) - 1
    
    COUNT_AGGREGATOR_CACHE[seqLen] = aggregator
    return aggregator
    
def normalize_wmm(wmm):
    """
    Normalizes and returns the provided Weight Matrix Model
    Each column of the returned model will add up to one
    Note: NAN's are not expected and therefore not handled
    """
    
    wmm = numpy.array(wmm)
    assert wmm.shape[1] == WMM_LENGTH, "Weight matrix model must have %d columns" % WMM_LENGTH
    
    return wmm / numpy.sum(wmm, 0)
    
def apply_wmm_to_sequence(wmm, sequence):
    """
    Applies the given WMM to the given matrixified sequence 
        (See: matrixify_sequence(...))
    Returns a Numpy array of length sequence.shape()[1] - WMM_LENGTH + 1
        where each probabiltity corresponds to the probability 
        of the WMM matching the sequence at that index
    """
    
    '''
    The result of [sequence]^T * [wmm] will be a matrix like:
        | 1 0 0 0 0 0 |
        | 2 1 0 0 0 0 |
        | 3 2 1 0 0 0 |
        | 4 3 2 1 0 0 |
        | 5 4 3 2 1 0 |
        | 6 5 4 3 2 1 |
        | 7 6 5 4 3 2 |
        | 0 7 6 5 4 3 |
        | 0 0 7 6 5 4 |
        | 0 0 0 7 6 5 |
        | 0 0 0 0 7 6 |
        | 0 0 0 0 0 7 |
    Where 0 indicates a value that makes no sense as part of the WMM score
        and 1-7 indicate partial scores of a particular position in the sequence
    The result must be shifted up and summed before returning
    '''
    scores = numpy.dot(numpy.transpose(sequence), wmm)
    
    resultLength = scores.shape[0] - WMM_LENGTH + 1
    assert resultLength > 0, "Not enough data to apply the WMM against"
    
    # Shift columns upwards
    for index in range(1, WMM_LENGTH):
        scores[0:resultLength, index] = scores[index:(index + resultLength), index]
    
    # Add up the relevant rows and normalize the probabilities
    scores = numpy.sum(scores[0:resultLength, :], 1)
    return scores / numpy.sum(scores)
