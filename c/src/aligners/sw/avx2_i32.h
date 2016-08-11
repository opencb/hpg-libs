#ifndef SW_AVX2_I32_H
#define SW_AVX2_I32_H

#include <immintrin.h>

#include "macros.h"

//------------------------------------------------------------------------
// AVX2 functions
//------------------------------------------------------------------------

void avx2_matrix2_i32(int num_seqs,
		      char **q, int *q_len, int max_q_len,
		      char **r, int *r_len, int max_r_len,
		      int match, int mismatch,
		      int gap_open, int gap_extend,
		      int *H, int *F, int *C, int *max_score);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // SW_AVX2_I32_H
