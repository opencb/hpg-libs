#ifndef SW_AVX2_I16_H
#define SW_AVX2_I16_H

#include <immintrin.h>

#include "macros.h"

//------------------------------------------------------------------------
// AVX2 i16 functions
//------------------------------------------------------------------------

void avx2_matrix_i16(int num_seqs,
		     short *qq, int max_q_len,
		     short *rr, int max_r_len,
		     short match, short mismatch, short gap_open, short gap_extend,
		     short *H, short *F, int *C, short *max_score);


//------------------------------------------------------------------------

void simd_traceback_i16(int depth, int num_seqs,
			char **q, int *q_len, int max_q_len,
			char **r, int *r_len, int max_r_len,
			short gap_open, short gap_extend,
			short *H, short *C,
			short *max_score,
			char **q_alig, int *q_start,
			char **r_alig, int *r_start,
			int *len_alig,
			char *q_aux, char *r_aux);

//------------------------------------------------------------------------

void simd_find_position_i16(int depth, int index, char *q, int q_len, char *r, int r_len,
			    short *H, int cols, int rows, short score,
			    int *q_pos, int *r_pos);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // SW_AVX2_I16_H
