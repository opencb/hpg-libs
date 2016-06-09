#!/bin/bash


for simd in SSE AVX2; do
    make clean ; make METHOD=OMP TIMING=1 SIMD=$simd
    for b in 128 1024 4096; do
	for n in 1 2 4 6 8 10 12 14 16; do
	    echo "running $simd batch-size:$b num-threads:$n"
	    ./bin/hpg-sw -q /ramdisks/queries.10M.txt -r /ramdisks/refs.10M.txt -d /ramdisks/ -o out.results -s dnafull -b $b -p 10 -e 0.5 -n $n > results/${simd}_${b}_${n}.times
	    exit
	done 
    done 
done 

# use predefined variables to access passed arguments
#echo arguments to the shell
#echo $1 $2 $3 ' -> echo $1 $2 $3'

# We can also store arguments from bash command line in special array
args=("$@")
#echo arguments to the shell
#echo ${args[0]} ${args[1]} ${args[2]} ' -> args=("$@"); echo ${args[0]} ${args[1]} ${args[2]}'

#use $@ to print out all arguments at once
#echo $@ ' -> echo $@'

# use $# variable to print out
# number of arguments passed to the bash script
#echo Number of arguments passed: $# ' -> echo Number of arguments passed: $#' 
