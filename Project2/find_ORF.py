import sys
import os
import argparse
import re

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

def compare_ORFs_simple(sequence, ORFs, annotations):
    """
    Declares an ORF to be a "gene"
        iff the stop index matches a stop index of an annotation
    Prints out how many ORFs of a given length are and are not "genes"
    """

    # Extract all the stop codons within the annotated genes
    annot_stops = set()
    for annot in annotations:
        stops = _find_stops(sequence[annot[0]:(annot[1] + 3)])
        stops = [annot[0] + stop for stop in stops]
        annot_stops = annot_stops.union(stops)

    # Store the comparison in a map from ORF length to a tuple
    #   of (# of genes, # of not genes)
    comparison = {}

    for ORF in ORFs:
        length = ORF[1] - ORF[0]
        if length not in comparison:
            comparison[length] = (0, 0)

        if ORF[1] in annot_stops:
            comparison[length] = (comparison[length][0] + 1, comparison[length][1])
        else:
            comparison[length] = (comparison[length][0], comparison[length][1] + 1)

    # Output the "histogram" in LaTeX (PGFPlots package)
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
            file.write('(%d, %d)' % (comparison[length][0], length))
        file.write('}\n\\closedcycle;\n')
        file.write('\\addlegendentry{Match}\n')
        file.write('\\addplot coordinates {\n')
        for length in sorted(comparison.keys()):
            file.write('(%d, %d)' % (comparison[length][1], length))
        file.write('}\n\\closedcycle;\n')
        file.write('\\addlegendentry{No match}\n')
        file.write('\\end{axis}\n')
        file.write('\\end{tikzpicture}\n')
    
    # Also output the plain old text
    print 'ORF length: Match - No Match'
    for length in sorted(comparison.keys()):
        print "%*d: %*d - %d" % (10, length, 5, comparison[length][0], comparison[length][1])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Finds open reading frames in the given sequence')
    parser.add_argument('sequence', type=str)
    parser.add_argument('annotations', type=str)
    parser.add_argument('--simple', action='store_true')
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

    if args.simple:
        compare_ORFs_simple(sequence, ORFs, annotations)
