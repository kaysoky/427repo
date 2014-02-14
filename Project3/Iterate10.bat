:: 10 iterations
python viterbi.py TestViterbi.fna --inE Emission00.json --inT Transition00.json --outE Emission01.json --outT Transition01.json --time > Results00.out
python viterbi.py TestViterbi.fna --inE Emission01.json --inT Transition01.json --outE Emission02.json --outT Transition02.json --time > Results01.out
python viterbi.py TestViterbi.fna --inE Emission02.json --inT Transition02.json --outE Emission03.json --outT Transition03.json --time > Results02.out
python viterbi.py TestViterbi.fna --inE Emission03.json --inT Transition03.json --outE Emission04.json --outT Transition04.json --time > Results03.out
python viterbi.py TestViterbi.fna --inE Emission04.json --inT Transition04.json --outE Emission05.json --outT Transition05.json --time > Results04.out
python viterbi.py TestViterbi.fna --inE Emission05.json --inT Transition05.json --outE Emission06.json --outT Transition06.json --time > Results05.out
python viterbi.py TestViterbi.fna --inE Emission06.json --inT Transition06.json --outE Emission07.json --outT Transition07.json --time > Results06.out
python viterbi.py TestViterbi.fna --inE Emission07.json --inT Transition07.json --outE Emission08.json --outT Transition08.json --time > Results07.out
python viterbi.py TestViterbi.fna --inE Emission08.json --inT Transition08.json --outE Emission09.json --outT Transition09.json --time > Results08.out
python viterbi.py TestViterbi.fna --inE Emission09.json --inT Transition09.json --verbose > Results09.out
