:: This file contains all the necessary commands to generate 
::   all the derivative files referenced or used in the report

:: Help for each tool
python filter.py -h | tee Usage-Filter.out
python scanner.py -h | tee Usage-Scanner.out
python meme.py -h | tee Usage-Meme.out
python entropy.py -h | tee Usage-Entropy.out

:: Filter
python filter.py --min_mismatch 3 --min_polyAlen 3 --compute_background 00-01-background.json --verbose all.sam 00-01-MMM3-MPA3.json | tee 00-01-MMM3-MPA3.out
python filter.py --dereverse --verbose 00-01-MMM3-MPA3.json 01-02-DRV.json
python filter.py --min_polyAlen 10 --verbose 01-02-DRV.json 02-03-MPA10.json | tee 02-03-MPA10.out
python filter.py --min_polyAlen 10 --verbose 00-01-MMM3-MPA3.json 01-02-TEST.json
python filter.py --max_non_tail_mismatches 4 --min_UTRlen 18 --verbose 02-03-MPA10.json 03-04-NTM4-MUL18.json | tee 03-04-NTM4-MUL18.out

:: WMM0
python scanner.py --verbose 00-WMM0.json 00-01-background.json 03-04-NTM4-MUL18.json 03-04-WMM0-background.tsv | tee 03-04-WMM0-background.out
python entropy.py 00-WMM0.json 00-01-background.json | tee 00-WMM0-background-entropy.out
python scanner.py --verbose 00-WMM0.json 00-UniformBackground.json 03-04-NTM4-MUL18.json 03-04-WMM0-uniform.tsv | tee 03-04-WMM0-uniform.out
python entropy.py 00-WMM0.json 00-UniformBackground.json | tee 00-WMM0-uniform-entropy.out

:: WMM1
python scanner.py --verbose 00-WMM1.json 00-01-background.json 03-04-NTM4-MUL18.json 03-04-WMM1-background.tsv | tee 03-04-WMM1-background.out
python entropy.py 00-WMM1.json 00-01-background.json | tee 00-WMM1-background-entropy.out
python scanner.py --verbose 00-WMM1.json 00-UniformBackground.json 03-04-NTM4-MUL18.json 03-04-WMM1-uniform.tsv | tee 03-04-WMM1-uniform.out
python entropy.py 00-WMM1.json 00-UniformBackground.json | tee 00-WMM1-uniform-entropy.out

:: WMM2
python meme.py --verbose 00-WMM1.json 03-04-NTM4-MUL18.json 00-WMM2.json | tee 00-WMM2.out
python scanner.py --verbose 00-WMM2.json 00-01-background.json 03-04-NTM4-MUL18.json 03-04-WMM2-background.tsv | tee 03-04-WMM2-background.out
python entropy.py 00-WMM2.json 00-01-background.json | tee 00-WMM2-background-entropy.out
python scanner.py --verbose 00-WMM2.json 00-UniformBackground.json 03-04-NTM4-MUL18.json 03-04-WMM2-uniform.tsv | tee 03-04-WMM2-uniform.out
python entropy.py 00-WMM2.json 00-UniformBackground.json | tee 00-WMM2-uniform-entropy.out

:: Conclusions
python scanner.py --verbose 00-01-background.json 00-UniformBackground.json 03-04-NTM4-MUL18.json 03-04-background-uniform.tsv
