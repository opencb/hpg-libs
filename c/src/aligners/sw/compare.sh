#!/bin/bash

#for simd in SSE AVX2 AVX512; do
for simd in SSE AVX2; do
    make clean ; make TIMING=1 SIMD=$simd COMPILER=$1
    for b in 128 1024 4096; do
	n=1
	echo "running $simd batch-size:$b num-threads:$n"
	./bin/hpg-sw -q $2 -r $3 -s dnafull -b $b -p 10 -e 0.5 -n $n > results/${simd}_${b}_${n}.times
#	for n in {2..272..2}; do
	for n in {2..16..2}; do
	    echo "running $simd batch-size:$b num-threads:$n"
	    ./bin/hpg-sw -q $2 -r $3 -s dnafull -b $b -p 10 -e 0.5 -n $n > results/${simd}_${b}_${n}.times
	done 
    done 
done 

