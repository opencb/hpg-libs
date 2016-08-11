#include "sw_commons_i8.h"

//-------------------------------------------------------------------------
// sw_mem_i8_t functions
//-------------------------------------------------------------------------

sw_mem_i8_t *sw_mem_i8_new() {
  const int len = 250;
  sw_mem_i8_t *mem = (sw_mem_i8_t *) malloc(sizeof(sw_mem_i8_t));
  mem->q_size = len;
  mem->qq = (char *) _mm_malloc(SIMD_DEPTH * len * sizeof(char), SIMD_ALIGNED);

  mem->r_size = len;
  mem->rr = (char *) _mm_malloc(SIMD_DEPTH * len * sizeof(char), SIMD_ALIGNED);

  mem->H_size = len * len;
  mem->H = (char *) _mm_malloc(SIMD_DEPTH * mem->H_size * sizeof(char), SIMD_ALIGNED);
  mem->C = (char *) _mm_malloc(SIMD_DEPTH * mem->H_size * sizeof(char), SIMD_ALIGNED);

  mem->F_size = len;
  mem->F = (char *) _mm_malloc(SIMD_DEPTH * len * sizeof(char), SIMD_ALIGNED);

  mem->aux_size = len;
  mem->q_aux = (char *) calloc(len * 2, sizeof(char));
  mem->r_aux = (char *) calloc(len * 2, sizeof(char));

  for(int i=0; i < SIMD_DEPTH; i++) {
    mem->q_lens[i] = 0;
    mem->r_lens[i] = 0;
  }

  return mem;
}

//-------------------------------------------------------------------------

void sw_mem_i8_free(sw_mem_i8_t *mem) {
  _mm_free(mem->qq);
  _mm_free(mem->rr);
  _mm_free(mem->H);
  _mm_free(mem->C);
  _mm_free(mem->F);
  free(mem->q_aux);
  free(mem->r_aux);
  free(mem);
}

//-------------------------------------------------------------------------

void sw_mem_i8_alloc(char **q, char **r, int num_seqs,
		      sw_mem_i8_t *mem) {
  int len, i, j;
  int max_q_len = 0;
  int max_r_len = 0;
  for(int i=0; i < num_seqs; i++) {
    len = strlen(q[i]);
    mem->q_lens[i] = len;
    if (len > max_q_len) max_q_len = len;
    mem->q[i] = q[i];

    len = strlen(r[i]);
    mem->r_lens[i] = len;
    if (len > max_r_len) max_r_len = len;
    mem->r[i] = r[i];
  }

  // check if we have room, otherwise we have to allocate more memory
  mem->max_q_len = max_q_len;
  if (max_q_len > mem->q_size) {
    _mm_free(mem->qq);
    mem->qq = (char *) _mm_malloc(SIMD_DEPTH * max_q_len * sizeof(char), SIMD_ALIGNED);
    mem->q_size = max_q_len;
  }

  mem->max_r_len = max_r_len;
  if (max_r_len > mem->r_size) {
    _mm_free(mem->rr);
    mem->rr = (char *) _mm_malloc(SIMD_DEPTH * max_r_len * sizeof(char), SIMD_ALIGNED);
    mem->r_size = max_r_len;
  }

  // copy seqs
  for(i=0; i < SIMD_DEPTH; i++) {

    for(j=0; j < mem->q_lens[i]; j++) {
      mem->qq[SIMD_DEPTH * j + i] = (char) q[i][j];
    }
    for(; j < max_q_len; j++) {
      mem->qq[SIMD_DEPTH * j + i] = ';';
    }
    for(j=0; j < mem->r_lens[i]; j++) {
      mem->rr[SIMD_DEPTH * j + i] = (char) r[i][j];
    }
    for(; j < max_r_len; j++) {
      mem->rr[SIMD_DEPTH * j + i] = ':';
    }
  }

  // check room for H and C matrices
  int matrix_size = max_q_len * max_r_len;
  if (matrix_size > mem->H_size) {
    _mm_free(mem->H);
    _mm_free(mem->C);

    int size = SIMD_DEPTH * matrix_size;
    mem->H = (char *) _mm_malloc(size * sizeof(char), SIMD_ALIGNED);
    mem->C = (char *) _mm_malloc(size * sizeof(char), SIMD_ALIGNED);
    //printf("new H %x, C %x\n", *H, *C);
    mem->H_size = matrix_size;
  }
  //  memset(*H, 0, size_h);
  //  memset(*C, 0, size_c);

  // check room for F vector
  if (max_q_len > mem->F_size) {
    _mm_free(mem->F);
    mem->F = (char *) _mm_malloc(SIMD_DEPTH * max_q_len * sizeof(char), SIMD_ALIGNED);
    mem->F_size = max_q_len;
  }
  //  memset(*F, 0, size_f);

  // check room for alignment sizes 
  int max_size = (max_r_len > max_q_len ? max_r_len : max_q_len);
  if (max_size > mem->aux_size) {
    free(mem->q_aux);
    free(mem->r_aux);

    mem->q_aux = (char *) calloc(max_size * 2, sizeof(char));
    mem->r_aux = (char *) calloc(max_size * 2, sizeof(char));
    //printf("new q_aux %x, r_aux %x\n", *q_aux, *r_aux);
    mem->aux_size = max_size;
  }
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


