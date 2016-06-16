#ifndef SW_AVX512_H
#define SW_AVX512_H

#include <immintrin.h>

#include "macros.h"

//------------------------------------------------------------------------
// AVX512 functions
//------------------------------------------------------------------------

void avx512_matrix2(int num_seqs,
		    char **q, int *q_len, int max_q_len,
		    char **r, int *r_len, int max_r_len,
		    float match, float mismatch,
		    float gap_open, float gap_extend,
		    float *H, float *F, int *C, float *max_score);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // SW_AVX512_H
