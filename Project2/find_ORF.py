import sys
import os
import argparse
import re

from math import log

"""
File extension that triggers some additional processing
Lines starting with '>' are removed
    and whitespace is stripped out
"""
FASTA = '.fna'

"""
File extension that triggers some additional processing
All lines that do not start with 'CDS' are removed
"""
GENEBANK = '.gbk'

"""
Threshold beyond which an ORF is considered to be a gene
"""
GENE_THRESHOLD = 1400

"""
Threshold below which an ORF is considered to not be a gene
"""
NOT_GENE_THRESHOLD = 50

"""
Default Markov chain degree
"""
MARKOV_CHAIN_DEGREE = 3

"""
Key into a Markov chain matrix that denotes
the probability of starting with the given key
"""
START = 'START'

"""
The four possible nucleotides
"""
NUCLEOTIDES = ['A', 'C', 'G', 'T']

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

def process_genebank(text):
    """Extracts all coding sequences from the genebank file"""

    lines = text.split('\n')
    lines = [line.strip() for line in lines]
    lines = filter(lambda x: x.startswith('CDS'), lines)

    # Ignore complementary strands as per the assignment
    lines = filter(lambda x: "complement" not in x, lines)

    # Extract the coding sequence indices
    regex = re.compile('(\d*)\.\.(\d*)')
    ORFs = [regex.search(line).groups() for line in lines]

    # Transform the sequences to integers and 0-based indices
    return [(int(ORF[0]) - 1, int(ORF[1]) - 1) for ORF in ORFs]

def _find_stops(sequence):
    """
    Helper for find_ORFs
    Returns the indices of all stop codons TAA, TAG, or TGA
    The returned list is in ascending order
    """

    matches = []
    regex = re.compile('(TAA|TAG|TGA)')
    match = regex.search(sequence)
    while (match is not None):
        index = match.start()
        matches.append(index)
        match = regex.search(sequence, index + 1)

    # This should already be in ascending order
    matches.sort()
    return matches

def _find_ORFs_offset(stops, offset):
    """
    Helper for find_ORFs
    Takes a list of stop codon positions
    Returns a list of the open reading frames
        starting at some offset
    See find_ORFs for the return format
    """

    # Filter out other offsets and adjacent stops
    stops = filter(lambda x: (x % 3) == offset, stops)
    stops = [stops[0]] + map(lambda x: stops[x],
                    filter(lambda x: (stops[x] - 3) != stops[x - 1],
                            range(1, len(stops))))
    return map(lambda x: (offset, stops[x]) if x == 0 else (stops[x - 1] + 3, stops[x]),
                range(len(stops)))

def find_ORFs(sequence):
    """
    Returns a sorted list of the open reading frames
        Tuple format: (start index, end index)
        The end excludes the stop codon
    """

    stops = _find_stops(sequence)
    ORFs =        _find_ORFs_offset(stops, 0)
    ORFs = ORFs + _find_ORFs_offset(stops, 1)
    ORFs = ORFs + _find_ORFs_offset(stops, 2)
    ORFs.sort(key=lambda x: x[0])

    return ORFs

def _compute_markov_chain(sequence, ORFs, degree):
    """
    Helper for compare_ORFs
    Looks up the ORFs in the sequence
        and calculates the posterior transition probabilities of each state
    The degree determines the total number of possible states
    Returned values are all log-probabilities
    """

    # Setup the "matrix" of counts
    # Since we're working with strings,
    #   the matrix is really a hash table of hash tables
    #   of bounded size
    counts = {}
    for level in range(len(NUCLEOTIDES) ** degree):
        key = []
        for base in range(1, degree + 1):
            key.append(NUCLEOTIDES[
                (level % (len(NUCLEOTIDES) ** base))
                / len(NUCLEOTIDES) ** (base - 1)])
        key = ''.join(key)
        counts[key] = {}
        for base in NUCLEOTIDES:
            counts[key][base] = 0

    # Insert all the counts into the matrix
    for ORF in ORFs:
        ORF = sequence[ORF[0]:ORF[1]]
        for index in range(len(ORF) - degree - 1):
            key = ORF[index:(index + degree)]
            next = ORF[index + degree + 1]
            counts[key][next] += 1

    # Calculate the probability of starting with a particular sequence
    total_counts = 0
    for start in counts:
        count = 0
        for base in NUCLEOTIDES:
            count += counts[start][base]
        counts[start][START] = count
        total_counts += count
    total_counts = float(total_counts)
    for start in counts:
        counts[start][START] = log(counts[start][START] / total_counts)

    # Transform counts into probabilities
    for start in counts:
        total_counts = 0
        for base in NUCLEOTIDES:
            total_counts += counts[start][base]
        total_counts = float(total_counts)
        for base in NUCLEOTIDES:
            if counts[start][base] > 0:
                counts[start][base] = log(counts[start][base] / total_counts)

    return counts
    
def _calculate_log_ratio(sequence, gene_probs, not_gene_probs):
    """
    Helper for compare_ORFs
    Takes a sequence and two Markov chain probability matrices
        and calculates the log ratio of probabilities
    """
    
    # Determine the degree of the chain
    # We assume that the inputs are correct
    #   i.e. the keys are the same length, 
    #        all keys of the given length exist, 
    #        all probabilities are logged, 
    #        and both matrices are similar
    degree = len(gene_probs.keys()[0])
    
    # Now sum up all the associated log probabilities
    key = sequence[0:degree]
    ratio = gene_probs[key][START] - not_gene_probs[key][START]
    for index in range(len(sequence) - degree - 1):
        key = sequence[index:(index + degree)]
        next = sequence[index + degree + 1]
        ratio += gene_probs[key][next] - not_gene_probs[key][next]
    
    return ratio

def compare_ORFs(sequence, ORFs, annotations, output_LaTeX):
    """
    Declares an ORF to be a "gene"
        iff the stop index matches a stop index of an annotation
    Also declares an ORF to be a "gene" based on Markov chains
    Prints out how many ORFs of a given length are and are not "genes"
        and the average log ratio of Markov chain probabilities
        and the number of ORFs with positive log ratios
    """
    
    # Calculate the Markov chain probabilities
    gene_probs = _compute_markov_chain(sequence,
            filter(lambda x: (x[1] - x[0]) > GENE_THRESHOLD, ORFs),
            MARKOV_CHAIN_DEGREE)
    not_gene_probs = _compute_markov_chain(sequence,
            filter(lambda x: (x[1] - x[0]) < NOT_GENE_THRESHOLD, ORFs),
            MARKOV_CHAIN_DEGREE)

    # Extract all the stop codons within the annotated genes
    annot_stops = set()
    for annot in annotations:
        stops = _find_stops(sequence[annot[0]:(annot[1] + 3)])
        stops = [annot[0] + stop for stop in stops]
        annot_stops = annot_stops.union(stops)

    # Store the comparison in a map from ORF length to a hash
    comparison = {}
    SIMPLE_GENE = 'SIMPLE_GENE'
    NOT_SIMPLE_GENE = 'NOT_SIMPLE_GENE'
    AVERAGE_LOG_RATIO = 'AVERAGE_LOG_RATIO'
    POSITIVE_LOG_RATIO = 'POSITIVE_LOG_RATIO'
    POSITIVE_HIT = 'POSITIVE_HIT'

    for ORF in ORFs:
        length = ORF[1] - ORF[0]
        if length not in comparison:
            comparison[length] = {SIMPLE_GENE:0, 
                                  NOT_SIMPLE_GENE:0, 
                                  AVERAGE_LOG_RATIO:0, 
                                  POSITIVE_LOG_RATIO:0, 
                                  POSITIVE_HIT:0}

        # Increment the counts for the simple heuristic
        if ORF[1] in annot_stops:
            comparison[length][SIMPLE_GENE] += 1
        else:
            comparison[length][NOT_SIMPLE_GENE] += 1
            
        # Update the values for the Markov heuristic
        ratio = _calculate_log_ratio(sequence[ORF[0]:ORF[1]], gene_probs, not_gene_probs)
        comparison[length][AVERAGE_LOG_RATIO] += ratio
        if ratio > 0:
            comparison[length][POSITIVE_LOG_RATIO] += 1
            if ORF[1] in annot_stops:
                comparison[length][POSITIVE_HIT] += 1
    
    # Average out the AVERAGE_LOG_RATIO field
    for length in comparison.keys():
        comparison[length][AVERAGE_LOG_RATIO] /= comparison[length][SIMPLE_GENE] + comparison[length][NOT_SIMPLE_GENE]

    # Output the "histogram" in LaTeX (PGFPlots package)
    # This will not include the results of the Markov heuristic
    if output_LaTeX:
        with open('Histogram.tex', 'w') as file:
            file.write('\\begin{tikzpicture}\n')
            file.write('\\begin{axis}[stack plots=x, '
                    'area style, '
                    'enlarge x limits=false, '
                    'enlarge y limits=false, '
                    'xmode=log, '
                    'ymode=log, '
                    'xlabel=Number of ORFs, '
                    'ylabel=Length of ORF, '
                    'width=\\textwidth, '
                    'height=\\textheight'
                    ']\n')
            file.write('\\addplot coordinates\n{')
            for length in sorted(comparison.keys()):
                file.write('(%d, %d)' % (comparison[length][SIMPLE_GENE], length))
            file.write('}\n\\closedcycle;\n')
            file.write('\\addlegendentry{Match}\n')
            file.write('\\addplot coordinates {\n')
            for length in sorted(comparison.keys()):
                file.write('(%d, %d)' % (comparison[length][NOT_SIMPLE_GENE], length))
            file.write('}\n\\closedcycle;\n')
            file.write('\\addlegendentry{No match}\n')
            file.write('\\end{axis}\n')
            file.write('\\end{tikzpicture}\n')

    # Also output the plain old text
    print 'ORF length: Match | No Match | Average Log Ratio | Positive Log Ratio | Positive Hits'
    for length in sorted(comparison.keys()):
        print '%*d: %*d | %s | %*f | %*d | %d' % (10, length, 
                                  5, comparison[length][SIMPLE_GENE], 
                                  str(comparison[length][NOT_SIMPLE_GENE]).ljust(8), 
                                  17, comparison[length][AVERAGE_LOG_RATIO], 
                                  18, comparison[length][POSITIVE_LOG_RATIO], 
                                  comparison[length][POSITIVE_HIT])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Finds open reading frames in the given sequence')
    parser.add_argument('sequence', type=str)
    parser.add_argument('annotations', type=str)
    parser.add_argument('--LaTeX', action='store_true')
    args = parser.parse_args()

    # Read the sequence in as a string
    with open(args.sequence) as f:
        sequence = f.read().strip()

    # Read the annotation information in as a string
    with open(args.annotations) as f:
        annotations = f.read().strip()

    # Handle sequence files
    if os.path.splitext(args.sequence)[1] == FASTA:
        sequence = process_fasta(sequence)
    else:
        print 'Unknown sequence file format'
        exit()

    # Handle GeneBank files
    if os.path.splitext(args.annotations)[1] == GENEBANK:
        annotations = process_genebank(annotations)
    else:
        print 'Unknown genebank file format'
        exit()

    ORFs = find_ORFs(sequence)
    compare_ORFs(sequence, ORFs, annotations, args.LaTeX)
