#CC = icc
#CFLAGS = -std=c99 -O3 -openmp -Wall -Wno-missing-declarations -fast -mavx2 -vec-report1

#CC = gcc
ifeq ($(COMPILER), intel)
	CC = icc 
else
	CC = gcc
endif

#CC = /opt/intel/bin/icc

SW_VER = SW_VERSION_2
ifeq ($(SW_VERSION), 1)
    SW_VER = SW_VERSION_1
endif

ifeq ($(SIMD), AVX2)
    CFLAGS = -std=c99 -O3 -fopenmp -g -mavx2 -DTIMING -DSW_$(SIMD) -D$(SW_VER)
    OBJS = smith_waterman.o sse.o avx2.o macros.o
else
    CFLAGS = -std=c99 -O3 -fopenmp -g -msse4.2 -DTIMING -DSW_$(SIMD) -D$(SW_VER)
    OBJS = smith_waterman.o sse.o macros.o
endif

INCLUDES = -I .
LIBS = -lxml2 -lm
TEST = test

all: hpg-sw hpg-sw-bench

hpg-sw: $(OBJS) $(TEST)/hpg_sw.h
	$(CC) $(CFLAGS) $(OBJS) $(TEST)/hpg_sw.c -o bin/hpg-sw $(INCLUDES) $(LIBS)

hpg-sw-bench: $(OBJS) emboss.o $(TEST)/hpg_sw_bench.h
	$(CC) $(CFLAGS) $(OBJS) emboss.o $(TEST)/hpg_sw_bench.c -o bin/hpg-sw-bench $(INCLUDES) $(LIBS)

macros.o: macros.c macros.h
	$(CC) $(CFLAGS) -c macros.c $(INCLUDES)

emboss.o: emboss.c emboss.h
	$(CC) $(CFLAGS) -c emboss.c $(INCLUDES)

sse.o: sse.c sse.h
	$(CC) $(CFLAGS) -c sse.c $(INCLUDES)

avx2.o: avx2.c avx2.h
	$(CC) $(CFLAGS) -c avx2.c $(INCLUDES)

smith_waterman.o: smith_waterman.c smith_waterman.h
	$(CC) $(CFLAGS) -c smith_waterman.c $(INCLUDES)

clean:
	-rm -rf main *.o bin/hpg-sw bin/hpg-sw-bench




