:: This file contains all the necessary commands to generate 
::   all the derivative files referenced or used in the report
python filter.py -h | tee Usage.out
python filter.py --min_mismatch 3 --min_polyAlen 3 --compute_background 00-01-background.json --verbose all.sam 00-01-MMM3-MPA3.json | tee 00-01-MMM3-MPA3.out