#ifndef SW_COMMONS_I8_H
#define SW_COMMONS_I8_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sw_commons.h"

typedef struct sw_mem_i8 {
  int max_q_len;
  int max_r_len;
  int q_size;
  int r_size;
  int q_lens[SIMD_DEPTH];
  int r_lens[SIMD_DEPTH];
  int H_size;
  int F_size;
  int aux_size;
  char *qq;
  char *rr;
  char *H;
  char *F;
  char *C;
  char *q_aux;
  char *r_aux;
  char *q[SIMD_DEPTH];
  char *r[SIMD_DEPTH];
} sw_mem_i8_t;

//-------------------------------------------------------------------------

sw_mem_i8_t *sw_mem_i8_new();
void sw_mem_i8_free(sw_mem_i8_t *mem);
void sw_mem_i8_alloc(char **q, char **r, int num_seqs, sw_mem_i8_t *mem);

//-------------------------------------------------------------------------


#endif // end of SW_COMMONS_I8_H
