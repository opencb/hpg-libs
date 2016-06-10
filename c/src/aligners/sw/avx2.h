#ifndef SW_AVX2_H
#define SW_AVX2_H

#include <immintrin.h>

#include "macros.h"

//------------------------------------------------------------------------
// AVX2 functions
//------------------------------------------------------------------------

void avx2_matrix(int num_seqs,
				char **q, int *q_len, int max_q_len,
				char **r, int *r_len, int max_r_len,
				float profile[128][128], float gap_open, float gap_extend,
				float *H, float *F, int *C, float *max_score);

//------------------------------------------------------------------------

void avx2_matrix2(int num_seqs,
				  char **q, int *q_len, int max_q_len,
				  char **r, int *r_len, int max_r_len,
				  float match, float mismatch,
				  float gap_open, float gap_extend,
				  float *H, float *F, int *C, float *max_score);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // SW_AVX2_H
