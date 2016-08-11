#ifndef SW_COMMONS_I32_H
#define SW_COMMONS_I32_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sw_commons.h"

typedef struct sw_mem_i32 {
  int max_q_len;
  int max_r_len;
  int q_size;
  int r_size;
  int q_lens[SIMD_DEPTH];
  int r_lens[SIMD_DEPTH];
  int H_size;
  int F_size;
  int aux_size;
  int *qq;
  int *rr;
  int *H;
  int *F;
  int *C;
  char *q_aux;
  char *r_aux;
  char *q[SIMD_DEPTH];
  char *r[SIMD_DEPTH];
} sw_mem_i32_t;

//-------------------------------------------------------------------------

sw_mem_i32_t *sw_mem_i32_new();
void sw_mem_i32_free(sw_mem_i32_t *mem);
void sw_mem_i32_alloc(char **q, char **r, int num_seqs, sw_mem_i32_t *mem);

//-------------------------------------------------------------------------


#endif // end of SW_COMMONS_I32_H
