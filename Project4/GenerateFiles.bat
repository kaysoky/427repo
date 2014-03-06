:: This file contains all the necessary commands to generate 
::   all the derivative files referenced or used in the report
python filter.py -h | tee Usage.out
python filter.py --min\_mismatch 3 --min\_polyAlen 3 --verbose all.sam 01-MMM3-MPA3.json | tee 01-MMM3-MPA3.out