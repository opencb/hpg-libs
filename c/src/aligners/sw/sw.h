#ifndef SW_H
#define SW_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <omp.h>

#include "sw_commons.h"
#include "sw_i8.h"
#include "sw_i16.h"

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

typedef struct sw_mem {
  int max_q_len;
  int max_r_len;
  int q_size;
  int r_size;
  int q_lens[32];
  int r_lens[32];
  int H_size;
  int F_size;
  int aux_size;
  short *qq;
  short *rr;
  short *H;
  short *F;
  short *C;
  char *q_aux;
  char *r_aux;
  char *q[32];
  char *r[32];
} sw_mem_t;

sw_mem_t *sw_mem_new();
void sw_mem_free(sw_mem_t *mem);

//------------------------------------------------------------------------------------

void sw_mqmr(char **query_p, char **ref_p, unsigned int num_queries,
	     sw_optarg_t *optarg_p, sw_multi_output_t *output_p,
	     sw_mem_t *mem);

//------------------------------------------------------------------------------------

#endif // SW_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
