#ifndef SW_I8_H
#define SW_I8_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <omp.h>

#include "sw_commons.h"

#include "sse.h"
#ifdef __AVX2__
#include "sw_commons_i8.h"
#endif
#ifdef __AVX512__
#include "avx512_i16.h"
#endif

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void sw_mqmr_i8(char **query_p, char **ref_p, unsigned int num_queries,
		sw_optarg_t *optarg_p, sw_multi_output_t *output_p,
		sw_mem_i8_t *mem);

//------------------------------------------------------------------------------------

void matrix_avx2_i8(int num_seqs,
		    char *qq, int *q_len, int max_q_len,
		    char *rr, int *r_len, int max_r_len,
		    char match, char mismatch, char gap_open, char gap_extend,
		    char *H, char *F, char *C, char *max_score);

//------------------------------------------------------------------------------------

void backtracking_i8(int depth, int num_seqs,
		     char **q, int *q_len, int max_q_len,
		     char **r, int *r_len, int max_r_len,
		     char gap_open, char gap_extend,
		     char *H, char *C,
		     char *max_score,
		     char **q_alig, int *q_start,
		     char **r_alig, int *r_start,
		     int *len_alig,
		     char *q_aux, char *r_aux);

//------------------------------------------------------------------------------------

void find_position_i8(int depth, int index, char *q, int q_len, char *r, int r_len,
		      char *H, int cols, int rows, char score,
		      int *q_pos, int *r_pos);

//------------------------------------------------------------------------------------

#endif // SW_I8_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
