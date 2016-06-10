#ifndef SMITH_WATERMAN_H
#define SMITH_WATERMAN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <omp.h>

#include "sw_commons.h"
#include "sse.h"
#ifdef __AVX2__
#include "avx2.h"
#endif

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void smith_waterman_mqmr(char **query_p, char **ref_p, unsigned int num_queries, 
			 sw_optarg_t *optarg_p, unsigned int num_threads, 
			 sw_multi_output_t *output_p);

//------------------------------------------------------------------------------------

void smith_waterman_mqsr(char **query_p, char *ref_p, unsigned int num_queries, 
			 sw_optarg_t *optarg_p, unsigned int num_threads, 
			 sw_multi_output_t *output_p);

//------------------------------------------------------------------------------------

void reallocate_memory(int max_q_len, int max_r_len, int simd_depth, 
		       int *H_size, float **H, int **C, int *F_size, float **F, 
		       int *aux_size, char **q_aux, char **r_aux);

//-------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------

#endif // SMITH_WATERMAN_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
