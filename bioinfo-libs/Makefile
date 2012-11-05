SIMD =

@if [ `grep avx /proc/cpuinfo` ];\
then\
        SIMD = -mavx\
fi

CC = gcc
CC = /opt/intel/bin/icc

CFLAGS = -std=c99 -O3 -fopenmp $(SIMD)

SW_DIR = ./sw

# global targets
all: build_sw

clean: clean_sw


# Smith-Waterman targets
build_sw:
	@echo "Invoking Smith-Waterman makefile to build"
	(cd $(SW_DIR); $(MAKE))

clean_sw:
	@echo "Invoking Smith-Waterman makefile to clean"
	(cd $(SW_DIR); $(MAKE) clean)

