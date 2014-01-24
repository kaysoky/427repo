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
    for i in range(len(FILES)):
        with open(FILES[i]) as f:
            text = f.read()
        sequences.append(align.process_fasta(text))
        
        # Since the extension isn't important (and gets in the way)
        # Delete the extension
        FILES[i] = FILES[i][0:-6]
    
    # Write the scores to a LaTeX table
    table = open("Proteins.tex", "w")
    table.write("\\scalebox{0.7}{\n")
    table.write("\\begin{tabular}{r|*{10}{c|}}\n")
    
    # Setup the first row of the table
    for i in range(len(FILES)):
        table.write("& %s " % FILES[i])
    table.write("\\\\ \\hline\n")
    
    for i in range(len(FILES)):
        # Write the first column of the table
        table.write("%s " % FILES[i])
        
        # Fill in blank spaces (for an upper triangular matrix)
        for j in range(i):
            table.write("& ")
            
        for j in range(i, len(FILES)):
            print '-----%s ~ %s-----' % (FILES[i], FILES[j])
            _, _, score, _ = align.do_main(sequences[i], sequences[j], FILES[i], FILES[j])
            
            # Write the score
            table.write("& %d" % score)
        
        # Start the next row
        table.write("\\\\ \\hline\n")
    
    table.write("\\end{tabular}}")
    table.close()
