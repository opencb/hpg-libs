
SIMD = 

@if [ `grep avx /proc/cpuinfo` ];\
then\
	SIMD = -mavx\
fi

#SIMD = -mavx

#CC = icc
#CFLAGS = -std=c99 -O3 -openmp -Wall -Wno-missing-declarations -fast -xavx -vec-report1

CC = gcc
CFLAGS = -std=c99 -O3 -fopenmp $(SIMD)

INCLUDES = -I.
LIBS = -lxml2 -lm

all: build-all

build-all: smith_waterman.o sse.o emboss.o macros.o main.h
	$(CC) $(CFLAGS) smith_waterman.o macros.o emboss.o sse.o main.c -o main $(INCLUDES) $(LIBS)

macros.o: macros.c macros.h
	$(CC) $(CFLAGS) -c macros.c $(INCLUDES)

emboss.o: emboss.c emboss.h
	$(CC) $(CFLAGS) -c emboss.c $(INCLUDES)

sse.o: sse.c sse.h
	$(CC) $(CFLAGS) -c sse.c $(INCLUDES)

avx.o: avx.c avx.h
	$(CC) $(CFLAGS) -c avx.c $(INCLUDES)

smith_waterman.o: smith_waterman.c smith_waterman.h
	$(CC) $(CFLAGS) -c smith_waterman.c $(INCLUDES)

clean:
	-rm -rf main *.o




