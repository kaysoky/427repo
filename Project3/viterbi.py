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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='TODO')
    parser.add_argument('sequence', type=str)
    parser.add_argument('--LaTeX', action='store_true')
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
    
    #TODO: do something
