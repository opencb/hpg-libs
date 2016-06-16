OBJS = smith_waterman.o sse.o avx2.o avx512.o macros.o sw_omp.o sw_commons.o

# SW version parameter
SW_VER = SW_VERSION_2
ifeq ($(SW_VERSION), 1)
    SW_VER = SW_VERSION_1
endif

# TIMING parameter
ifeq ($(TIMING), 1)
    FTIMING = -DTIMING
endif

# compiler parameter
ifeq ($(COMPILER), intel)
	CC = icc
	# SIMD parameter
	ifeq ($(SIMD), AVX512)
	    CFLAGS = -std=c99 -O3 -fopenmp -g -axCORE-AVX512 $(FTIMING) -D__AVX512__ -DSW_$(SIMD) -D$(SW_VER)
	else
		ifeq ($(SIMD), AVX2)
		    CFLAGS = -std=c99 -O3 -fopenmp -g -march=core-avx2 -ipo -no-prec-div $(FTIMING) -DSW_$(SIMD) -D$(SW_VER)
		else
		    CFLAGS = -std=c99 -O3 -fopenmp -g -msse4.2 $(FTIMING) -DSW_$(SIMD) -D$(SW_VER)
		endif
	endif
else
	CC = gcc
	# SIMD parameter
	ifeq ($(SIMD), AVX2)
	    CFLAGS = -std=c99 -O3 -fopenmp -g -mavx2 $(FTIMING) -DSW_$(SIMD) -D$(SW_VER)
	else
	    CFLAGS = -std=c99 -O3 -fopenmp -g -msse4.2 $(FTIMING) -DSW_$(SIMD) -D$(SW_VER)
	endif
endif

INCLUDES = -I .
LIBS = -lm
TEST = test

all: hpg-sw hpg-sw-bench

hpg-sw: $(OBJS) $(TEST)/hpg_sw.h
	$(CC) $(CFLAGS) $(OBJS) $(TEST)/hpg_sw.c -o bin/hpg-sw $(INCLUDES) $(LIBS)

hpg-sw-bench: $(OBJS) emboss.o $(TEST)/hpg_sw_bench.h
	$(CC) $(CFLAGS) $(OBJS) emboss.o $(TEST)/hpg_sw_bench.c -o bin/hpg-sw-bench $(INCLUDES) $(LIBS)

macros.o: macros.c macros.h
	$(CC) $(CFLAGS) -c macros.c $(INCLUDES)

sw_commons.o: sw_commons.c sw_commons.h
	$(CC) $(CFLAGS) -c sw_commons.c $(INCLUDES)

emboss.o: emboss.c emboss.h
	$(CC) $(CFLAGS) -c emboss.c $(INCLUDES)

sse.o: sse.c sse.h
	$(CC) $(CFLAGS) -c sse.c $(INCLUDES)

avx2.o: avx2.c avx2.h
	$(CC) $(CFLAGS) -c avx2.c $(INCLUDES)

avx512.o: avx512.c avx512.h
	$(CC) $(CFLAGS) -c avx512.c $(INCLUDES)

smith_waterman.o: smith_waterman.c smith_waterman.h
	$(CC) $(CFLAGS) -c smith_waterman.c $(INCLUDES)

sw_omp.o: sw_omp.c sw_omp.h
	$(CC) $(CFLAGS) -c sw_omp.c $(INCLUDES)

clean:
	-rm -rf main *.o bin/hpg-sw bin/hpg-sw-bench




