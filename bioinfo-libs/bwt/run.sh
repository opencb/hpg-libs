#!/bin/bash

make test1-cpu 2> toto.kk ; grep error toto.kk 
echo "HOMO SAPIENS RUN (chromosoma:20 position:327752-527886/1)"
./test1-cpu /home/hmartinez/human_popcnt/ /home/hmartinez/human_popcnt/ CCCGCCTTTGCTCGGCGGAGACAGCAGGCAGAGAGGTGAGCTTAGCCCTGCCCCACGCGCGGCCAGGCCCCAGCCCCAGCCCCTGGAGAACCCCCGCGCTCTGCCCGCATCCTCAGCCCG

echo "FRUIT FLY RUN 1 (chromosoma:2L position:11177498_11177561/1)"
./test1-cpu /home/hmartinez/fly/index/ /home/hmartinez/fly/index/ TATTTATTCGCAAATGCATAACTTGAAAGGCTATTTTGTAGATATTAAATGCAACTAATTTGCC 

echo "FRUIT FLY RUN 2 (chromosoma:3R position:15260115_15260182/1)"
./test1-cpu /home/hmartinez/fly/index/ /home/hmartinez/fly/index/ CACTTGGCTACGCAAGTTCAAATTACGCAAAAATCAAATTCAGGACCTGAGATTGCGTTC 
