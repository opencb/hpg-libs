SIMD =

@if [ `grep avx /proc/cpuinfo` ];\
then\
        SIMD = -mavx\
fi

CC = gcc
CFLAGS = -std=c99 -O3 -fopenmp $(SIMD)

INCLUDES = -I. 

all: smith_waterman.o smith_waterman_simd.o

smith_waterman.o:
	$(CC) $(CFLAGS) -c smith_waterman.c $(INCLUDES)

smith_waterman_simd.o:
	$(CC) $(CFLAGS) -c smith_waterman_simd.c $(INCLUDES)

clean:
	-rm smith_waterman.o smith_waterman_simd.o

