#include "sw_i16.h"

#ifdef __AVX512__
  const unsigned int SIMD_DEPTH = 32;
  const unsigned int SIMD_ALIGNED = 64;
#else
#ifdef __AVX2__
  const unsigned int SIMD_DEPTH = 16;
  const unsigned int SIMD_ALIGNED = 32;
#else
  const unsigned int SIMD_DEPTH = 8;
  const unsigned int SIMD_ALIGNED = 16;
#endif
#endif

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void sw_mqmr_i16(char **query_p, char **ref_p, unsigned int num_queries,
		 sw_optarg_t *optarg_p, sw_multi_output_t *output_p) {
#ifdef __AVX512__
  const unsigned int simd_depth = 32;
  const unsigned int aligned_mem = 64;
#else
#ifdef __AVX2__
  const unsigned int simd_depth = 16;
  const unsigned int aligned_mem = 32;
#else
  const unsigned int simd_depth = 8;
  const unsigned int aligned_mem = 16;
#endif
#endif

  if (output_p == NULL) {
    printf("Error: output buffer is null.\n");
    exit(-1);
  }
  if (output_p->num_queries < num_queries) {
    printf("Error: num. sequencias (%i) and output buffer length (%i) mismatch !\n",
	   num_queries, output_p->num_queries);
    exit(-1);
  }
  
  int tid = omp_get_thread_num();

#ifdef TIMING
  double partial_t;
#endif // TIMING
  char *q_aux = NULL, *r_aux = NULL;
  int depth, aux_size = 0, H_size = 0, F_size = 0, max_q_len = 0, max_r_len = 0;
  int *q, *r;
  int len, index, q_lens[simd_depth], r_lens[simd_depth], alig_lens[simd_depth];
  short *H = NULL, *F = NULL, *C = NULL;

  short match = optarg_p->match;
  short mismatch = optarg_p->mismatch;
  short gap_open = optarg_p->gap_open;
  short gap_extend = optarg_p->gap_extend;
  float *score_p = output_p->score_p;

  //printf("num queries = %i\n", num_queries);
  //for(int k = 0; k < 8; k++) printf("%0.2f ", score_p[k]);
  //printf("\n");

  depth = 0;
  for (int i = 0; i < num_queries; i++) {
    //printf("smith_waterman.c: query: #%i\n", i);
    len = strlen(query_p[i]);
    if (len > max_q_len) max_q_len = len;
    q_lens[depth] = len;
    q[depth] = query_p[i];
    
    len = strlen(ref_p[i]);
    if (len > max_r_len) max_r_len = len;
    r_lens[depth] = len;
    r[depth] = ref_p[i];
    
    depth++;
    if (depth == simd_depth) {
      index = i - depth + 1;

      reallocate_memory_i16(max_q_len, max_r_len, simd_depth, 
			    &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
      
      // generating score matrix
#ifdef TIMING
      partial_t = sw_tic();
#endif // TIMING

#ifdef __AVX2__
      avx2_matrix_i16(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
		      match, mismatch, gap_open, gap_extend, H, F, C, &score_p[index]);
#else
      printf("16-bit integer version: not yet implemented!!\n");
      printf("For the moment, only for AVX2.\n");
      exit(-1);
#endif // __AVX2__

#ifdef TIMING
      sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING

#ifdef TIMING
      partial_t = sw_tic();
#endif // TIMING

      // tracebacking
      simd_traceback_i16(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
	gap_open, gap_extend, H, C, &score_p[index],
	&output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
	&output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
	q_aux, r_aux);
#ifdef TIMING
      sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING
      depth = 0;
      max_q_len = 0;
      max_r_len = 0;
    }
  }
}
      /*
      //    printf("depth = %i\n", depth);

        if (depth > 0) {

            //float max_score[simd_depth];
            float *max_score = (float *) _mm_malloc(simd_depth * sizeof(float), aligned_mem);

            for (unsigned int i = depth; i < simd_depth; i++) {
                q[i] = q[0];
                q_lens[i] = q_lens[0];

                r[i] = r[0];
                r_lens[i] = r_lens[0];
            }

            reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);

            index = num_queries - depth;

            // generating score matrix
#ifdef TIMING
            partial_t = sw_tic();
#endif // TIMING

ma#ifdef __AVX512__
	    avx512_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
		           match, mismatch, gap_open, gap_extend, H, F, C, &score_p[index]);
#else
#ifdef __AVX2__
#ifdef SW_VERSION_1
            avx2_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                        optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#else
            avx2_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                         match, mismatch, gap_open, gap_extend, H, F, C, max_score);
#endif // SW_VERSION_1
#else
#ifdef SW_VERSION_1
            sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                       optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#else
            sse_matrix2(depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                        match, mismatch, gap_open, gap_extend, H, F, C, max_score);
#endif // SW_VERSION_1
#endif // __AVX2__
#endif // __AVX512__

#ifdef TIMING
            sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING

#ifdef TIMING
            partial_t = sw_tic();
#endif // TIMING
            // tracebacking
            simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
                           gap_open, gap_extend, H, C, max_score,
                           &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
                           &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
                           q_aux, r_aux);
#ifdef TIMING
            sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING

            for (unsigned int i = 0; i < depth; i++) {
                score_p[index +  i] = max_score[i];
            }
            _mm_free(max_score);

        }
        //    printf("Free 1\n");
        // free memory
        if (H != NULL) _mm_free(H);
        if (C != NULL) _mm_free(C);
        if (F != NULL) _mm_free(F);
        if (q_aux != NULL) free(q_aux);
        if (r_aux != NULL) free(r_aux);
      */


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

typedef struct sw_mem_i16 {
  int max_q_len;
  int max_r_len;
  int depth;
  int *H_size;
  float **H;
  int **C;
  int *F_size;
  float **F;
  int *aux_size; 
  char **q_aux; 
  char **r_aux;
} sw_mem_i16_t;

sw_mem_i16_t *sw_mem_new(int depth) {
  sw_mem_i16_t *mem = (sw_mem_i16_t *) malloc(sizeof(sw_mem_i16_t));

  mem->depth = depth;

  return mem;
}

//------------------------------------------------------------------------------------

void sw_mem_reallocate(char **q, char **r, int num_seqs, sw_mem_i16_t *mem) {

    int *qq = (int *) _mm_malloc(depth * max_q_len * sizeof(int), 32);
    int *rr = (int *) _mm_malloc(depth * max_r_len * sizeof(int), 32);

    int length, i, j;
    for(int i=0; i < num_seqs; i++) {
      length = strlen(q[i]);
      if (length > max_q_len) max_q_len = length;
      length = strlen(r[i]);
      if (length > max_r_len) max_r_len = length;
    }

    int size, matrix_size = max_q_len * max_r_len;
    if (matrix_size > mem->matrix_size) {
        if (mem->matrix_size > 0) {
            _mm_free(mem->H);
            _mm_free(mem->C);
        }

	size = matrix_size * SIMD_DEPTH * sizeof(short);
        mem->H = (short *) _mm_malloc(size, SIMD_ALIGNED);
        mem->C = (short *) _mm_malloc(size, SIMD_ALIGNED);
        //printf("new H %x, C %x\n", *H, *C);
        mem->matrix_size = matrix_size;
    }

    size = max_q_len * SIMD_DEPTH * sizeof(short);
    if (max_q_len > mem->max_q_len) {
      if (max_q_len > 0) {
	_mm_free(mem->F);
      }
      mem->F = (short *) _mm_malloc(size, SIMD_ALIGNED);
      //printf("new F %x\n", *F);
      mem->max_q_len = max_q_len;
    }

    max_length = (max_r_len > max_q_len ? max_r_len : max_q_len);
    if (max_length > mem->max_length) {
      if (mem->max_length > 0) {
	free(mem->q_aux);
	free(mem->r_aux);
      }
      mem->q_aux = (char *) calloc(max_length * 2, sizeof(char));
      mem->r_aux = (char *) calloc(max_length * 2, sizeof(char));
      //printf("new q_aux %x, r_aux %x\n", *q_aux, *r_aux);
      mem->max_length = max_length;
    }




    }
    size_h = simd_depth * matrix_size * sizeof(short);

    for(int i=0; i < num_seqs; i++) {
        for(jj=0; jj < q_len[ii]; jj++) {
            qq[mem->depth * jj + ii] = (int) q[ii][jj];
        }
        for(; jj < max_q_len; jj++) {
            qq[depth * jj + ii] = ';';
        }
        for(jj=0; jj < r_len[ii]; jj++) {
            rr[depth * jj + ii] = (int) r[ii][jj];
        }
        for(; jj < max_r_len; jj++) {
            rr[depth * jj + ii] = ':';
        }
    }



void reallocate_memory(int max_q_len, int max_r_len, int simd_depth,
                       int *H_size, float **H, int **C, int *F_size, float **F,
                       int *aux_size, char **q_aux, char **r_aux) {

    int size_h, size_c, size_f;
    unsigned int matrix_size = max_q_len * max_r_len;

    size_h = simd_depth * matrix_size * sizeof(float);
    size_c = simd_depth * matrix_size * sizeof(int);
    if (matrix_size > *H_size) {
        if (*H_size > 0) {
            _mm_free(*H);
            _mm_free(*C);
        }

        *H = (float *) _mm_malloc(size_h, aligned_mem);
        *C = (int *) _mm_malloc(size_c, aligned_mem);
        //printf("new H %x, C %x\n", *H, *C);
        *H_size = matrix_size;
    }
    //  memset(*H, 0, size_h);
    //  memset(*C, 0, size_c);

    size_f = simd_depth * max_q_len * sizeof(float);
    if (max_q_len > *F_size) {
        if (*F_size > 0) _mm_free(*F);
        *F = (float *) _mm_malloc(size_f, aligned_mem);
        //printf("new F %x\n", *F);
        *F_size = max_q_len;
    }
    //  memset(*F, 0, size_f);

    int max_size = (max_r_len > max_q_len ? max_r_len : max_q_len);
    if (max_size > *aux_size) {
        if (*aux_size > 0) {
            free(*q_aux);
            free(*r_aux);
        }
        *q_aux = (char *) calloc(max_size * 2, sizeof(char));
        *r_aux = (char *) calloc(max_size * 2, sizeof(char));
        //printf("new q_aux %x, r_aux %x\n", *q_aux, *r_aux);
        *aux_size = max_size;
    }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
