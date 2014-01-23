import align

FILES = ['P15172.fasta',
         'P17542.fasta',
         'P10085.fasta',
         'P16075.fasta',
         'P13904.fasta',
         'Q90477.fasta',
         'Q8IU24.fasta',
         'P22816.fasta',
         'Q10574.fasta',
         'O95363.fasta']

if __name__ == "__main__":
    # Load the files
    sequences = []
    for file in FILES:
        with open(file) as f:
            text = f.read()
        sequences.append(align.process_fasta(text))
    
    for i in range(len(FILES)):
        for j in range(i, len(FILES)):
            print '-----%s ~ %s-----' % (FILES[i][0:-6], FILES[j][0:-6])
            align.do_main(sequences[i], sequences[j], FILES[i], FILES[j])
