#ifndef SW_F32_H
#define SW_F32_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <omp.h>

#include "sw_commons.h"

#include "sse.h"
#ifdef __AVX2__
#include "avx2_i16.h"
#include "sw_commons_f32.h"
#endif
#ifdef __AVX512__
#include "avx512_i16.h"
#endif

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void sw_mqmr_f32(char **query_p, char **ref_p, unsigned int num_queries,
		 sw_optarg_t *optarg_p, sw_multi_output_t *output_p,
		 sw_mem_f32_t *mem);

//------------------------------------------------------------------------------------

void matrix_avx2_f32(int num_seqs,
		     int *qq, int *q_len, int max_q_len,
		     int *rr, int *r_len, int max_r_len,
		     float match, float mismatch, float gap_open, float gap_extend,
		     float *H, float *F, int *C, float *max_score);

//------------------------------------------------------------------------------------

void backtracking_f32(int depth, int num_seqs,
		      char **q, int *q_len, int max_q_len,
		      char **r, int *r_len, int max_r_len,
		      float gap_open, float gap_extend,
		      float *H, int *C,
		      float *max_score,
		      char **q_alig, int *q_start,
		      char **r_alig, int *r_start,
		      int *len_alig,
		      char *q_aux, char *r_aux);

//------------------------------------------------------------------------------------

void find_position_f32(int depth, int index, char *q, int q_len, char *r, int r_len,
		       float *H, int cols, int rows, float score,
		       int *q_pos, int *r_pos);

//------------------------------------------------------------------------------------

#endif // SW_F32_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
