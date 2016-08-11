#ifndef SW_I16_H
#define SW_I16_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <omp.h>

#include "sw_commons.h"

#include "sse.h"
#ifdef __AVX2__
#include "avx2_i16.h"
#include "sw_commons_i16.h"
#endif
#ifdef __AVX512__
#include "avx512_i16.h"
#endif

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void sw_mqmr_i16(char **query_p, char **ref_p, unsigned int num_queries,
		 sw_optarg_t *optarg_p, sw_multi_output_t *output_p,
		 sw_mem_i16_t *mem);

//------------------------------------------------------------------------------------

void matrix_avx2_i16(int num_seqs,
		     short *qq, int *q_len, int max_q_len,
		     short *rr, int *r_len, int max_r_len,
		     short match, short mismatch, short gap_open, short gap_extend,
		     short *H, short *F, short *C, short *max_score);

//------------------------------------------------------------------------------------

void backtracking_i16(int depth, int num_seqs,
		      char **q, int *q_len, int max_q_len,
		      char **r, int *r_len, int max_r_len,
		      short gap_open, short gap_extend,
		      short *H, short *C,
		      short *max_score,
		      char **q_alig, int *q_start,
		      char **r_alig, int *r_start,
		      int *len_alig,
		      char *q_aux, char *r_aux);

//------------------------------------------------------------------------------------

void find_position_i16(int depth, int index, char *q, int q_len, char *r, int r_len,
		       short *H, int cols, int rows, short score,
		       int *q_pos, int *r_pos);

//------------------------------------------------------------------------------------

#endif // SW_I16_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
