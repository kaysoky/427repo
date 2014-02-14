:: 10 iterations
python viterbi.py NC_000909.fna --inE Emission00.json --inT Transition00.json --outE Emission01.json --outT Transition01.json --time --verbose | tee Results00.out
python viterbi.py NC_000909.fna --inE Emission01.json --inT Transition01.json --outE Emission02.json --outT Transition02.json --time | tee Results01.out
python viterbi.py NC_000909.fna --inE Emission02.json --inT Transition02.json --outE Emission03.json --outT Transition03.json --time | tee Results02.out
python viterbi.py NC_000909.fna --inE Emission03.json --inT Transition03.json --outE Emission04.json --outT Transition04.json --time | tee Results03.out
python viterbi.py NC_000909.fna --inE Emission04.json --inT Transition04.json --outE Emission05.json --outT Transition05.json --time | tee Results04.out
python viterbi.py NC_000909.fna --inE Emission05.json --inT Transition05.json --outE Emission06.json --outT Transition06.json --time | tee Results05.out
python viterbi.py NC_000909.fna --inE Emission06.json --inT Transition06.json --outE Emission07.json --outT Transition07.json --time | tee Results06.out
python viterbi.py NC_000909.fna --inE Emission07.json --inT Transition07.json --outE Emission08.json --outT Transition08.json --time | tee Results07.out
python viterbi.py NC_000909.fna --inE Emission08.json --inT Transition08.json --outE Emission09.json --outT Transition09.json --time | tee Results08.out
python viterbi.py NC_000909.fna --inE Emission09.json --inT Transition09.json --verbose | tee Results09.out
