:: This file contains all the necessary commands to generate 
::   all the derivative files referenced or used in the report
python filter.py -h | tee Usage.out
python filter.py --min_mismatch 3 --min_polyAlen 3 --compute_background 00-01-background.json --verbose all.sam 00-01-MMM3-MPA3.json | tee 00-01-MMM3-MPA3.out
python filter.py --dereverse --verbose 00-01-MMM3-MPA3.json 01-02-DRV.json
python filter.py --min_polyAlen 10 --verbose 01-02-DRV.json 02-03-MPA10.json | tee 02-03-MPA10.out
python filter.py --min_polyAlen 10 --verbose 00-01-MMM3-MPA3.json 01-02-TEST.json
python filter.py --max_non_tail_mismatches 6 --min_UTRlen 18 --verbose 02-03-MPA10.json 03-04-NTM10-MUL18.json | tee 03-04-NTM10-MUL18.out
