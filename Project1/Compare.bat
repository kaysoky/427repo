:: Test case
python align.py TestA.txt TestB.txt -n 100 > AB.out

:: Protein cases
python compare_proteins.py > Proteins.out

:: Empirical p-value cases
python align.py P15172.fasta Q10574.fasta -n 2000 --verbose > Empirical_P15172_Q10574.out
python align.py P15172.fasta O95363.fasta -n 2000 --verbose > Empirical_P15172_O95363.out
