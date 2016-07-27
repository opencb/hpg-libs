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
#endif
#ifdef __AVX512__
#include "avx512_i16.h"
#endif

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void sw_mqmr_i16(char **query_p, char **ref_p, unsigned int num_queries,
		 sw_optarg_t *optarg_p, sw_multi_output_t *output_p);

//------------------------------------------------------------------------------------

void reallocate_memory_i16(int max_q_len, int max_r_len, int simd_depth, 
			   short *H_size, short **H, short **C, short *F_size, short **F, 
			   int *aux_size, char **q_aux, char **r_aux);

//------------------------------------------------------------------------------------

#endif // SW_I16_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
