#ifndef SW_I32_H
#define SW_I32_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <omp.h>

#include "sw_commons.h"

#include "sse.h"
#ifdef __AVX2__
#include "sw_commons_i32.h"
#endif
#ifdef __AVX512__
#include "avx512_i32.h"
#endif

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void sw_mqmr_i32(char **query_p, char **ref_p, unsigned int num_queries,
		 sw_optarg_t *optarg_p, sw_multi_output_t *output_p,
		 sw_mem_i32_t *mem);

//------------------------------------------------------------------------------------

void matrix_avx2_i32(int num_seqs,
		     int *qq, int *q_len, int max_q_len,
		     int *rr, int *r_len, int max_r_len,
		     int match, int mismatch, int gap_open, int gap_extend,
		     int *H, int *F, int *C, int *max_score);

//------------------------------------------------------------------------------------

void backtracking_i32(int depth, int num_seqs,
		      char **q, int *q_len, int max_q_len,
		      char **r, int *r_len, int max_r_len,
		      int gap_open, int gap_extend,
		      int *H, int *C,
		      int *max_score,
		      char **q_alig, int *q_start,
		      char **r_alig, int *r_start,
		      int *len_alig,
		      char *q_aux, char *r_aux);

//------------------------------------------------------------------------------------

void find_position_i32(int depth, int index, char *q, int q_len, char *r, int r_len,
		       int *H, int cols, int rows, int score,
		       int *q_pos, int *r_pos);

//------------------------------------------------------------------------------------

#endif // SW_I32_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
