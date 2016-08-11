#include "sw.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

typedef struct pair {
  int index;
  int query_len;
  int ref_len;
  char *query;
  char *ref;
} pair_t;

void pair_set(int index, char *q, char *r, pair_t *pair) {
  pair->index = index;
  pair->query_len = strlen(q);
  pair->ref_len = strlen(r);
  pair->query = q;
  pair->ref = r;
}

//-------------------------------------------------------------------------

int pair_cmp(const void *a, const void *b) {
  return ((pair_t *)a)->query_len -  ((pair_t *)b)->query_len;
}  

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int lower_bound(pair_t *pairs, int num_pairs, int val) {
  int lo = 0, hi = num_pairs - 1;
  while (lo < hi) {
    int mid = lo + (hi - lo)/2;
    if (pairs[mid].query_len < val)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void sw_mem_alloc(pair_t *pairs, int num_seqs, int depth, sw_mem_t *mem);
void sw_mem_print_seq(char *msg, int query, int index, int depth, sw_mem_t *mem);

//-------------------------------------------------------------------------

void sw_mqmr(char **query, char **ref, unsigned int num_queries,
	     sw_optarg_t *optarg, sw_multi_output_t *output,
	     sw_mem_t *mem) {

  if (output == NULL) {
    printf("Error: output buffer is null.\n");
    exit(-1);
  }
  if (output->num_queries < num_queries) {
    printf("Error: num. sequencias (%i) and output buffer length (%i) mismatch !\n",
	   num_queries, output->num_queries);
    exit(-1);
  }

  pair_t *pairs = (pair_t *) malloc(num_queries * sizeof(pair_t));
  printf("before sorting...\n");
  for(int i=0; i<num_queries; i++) {
    pair_set(i, query[i], ref[i], &pairs[i]);
    printf("%i -> q len = %i %s\n", i, pairs[i].query_len, pairs[i].query);
  }
  printf("sorting...\n");
  qsort(pairs, num_queries, sizeof(pair_t), pair_cmp);
  printf("after sorting...\n");
  for(int i=0; i<num_queries; i++) {
    printf("%i -> q len = %i %s\n", pairs[i].index, pairs[i].query_len, pairs[i].query);
  }

  int bound = lower_bound(pairs, num_queries, (256 / (int) optarg->match));
  printf("bound = %i\n", bound);
  
  int tid = omp_get_thread_num();

#ifdef TIMING
  double partial_t;
#endif // TIMING

  int alig_lens[32];
  
  short match = (short) optarg->match;
  short mismatch = (short) optarg->mismatch;
  short gap_open = (short) optarg->gap_open;
  short gap_extend = (short) optarg->gap_extend;
  short *score = (short*) output->score_p;

  
  //printf("num queries = %i\n", num_queries);
  //for(int k = 0; k < 8; k++) printf("%0.2f ", score_p[k]);
  //printf("\n");

  int i, j, k;
  int num_seqs, depth;

  // first, short sequences, computing 32 8-bit matrices
  char score_i8[32];
  depth = 32;
  printf("first, short reads\n");
  for(i = 0; i < bound; i += depth) {
    num_seqs = bound - i;
    if (num_seqs > depth) num_seqs = depth; 
    printf("\tprocessing %i: from %i to %i\n", num_seqs, i, i + num_seqs - 1);
    sw_mem_alloc(&pairs[i], num_seqs, 32, mem);

    //printf("-----> after alloc\n");

    // generating score matrix
#ifdef TIMING
    partial_t = sw_tic();
#endif // TIMING
    
    matrix_avx2_i8(num_seqs, (char *) mem->qq, mem->q_lens, mem->max_q_len, 
		   (char *) mem->rr, mem->r_lens, mem->max_r_len,
		   match, mismatch, gap_open, gap_extend, 
		   (char *) mem->H, (char *) mem->F, (char *) mem->C, score_i8);

    //printf("start avx2_matrix (index %i)\n", index);
    //printf("end avx2_matrix\n");
    // update scores
    for(j = 0, k = i; j < num_seqs; j++, k++) {
      score[k] = (short) score_i8[j];
      printf("score %i: %i\n", k, score_i8[j]);
    }

#ifdef TIMING
    sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING
    
#ifdef TIMING
    partial_t = sw_tic();
#endif // TIMING

    // backtracking
    backtracking_i8(depth, num_seqs, 
		    mem->q, mem->q_lens, mem->max_q_len, 
		    mem->r, mem->r_lens, mem->max_r_len,
		    gap_open, gap_extend, 
		    (char *) mem->H, (char *) mem->C, score_i8,
		    &output->query_map_p[i], (int *)&output->query_start_p[i],
		    &output->ref_map_p[i], (int *)&output->ref_start_p[i], 
		    alig_lens,
		    mem->q_aux, mem->r_aux);

#ifdef TIMING
    sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING
  }
  /*
  // second, long sequences, computing 16 16-bit matrices
  depth = 16;
  printf("second, long reads\n");
  for(int i = bound; i < num_queries; i += depth) {
    num_seqs = num_queries - i;
    if (num_seqs > depth) num_seqs = depth; 
    printf("\tprocessing %i: from %i to %i\n", num_seqs, i, i + num_seqs - 1);
    sw_mem_alloc(&pairs[i], num_seqs, 16, mem);

    //    printf("-----> after alloc\n");

    // generating score matrix
#ifdef TIMING
    partial_t = sw_tic();
#endif // TIMING
    matrix_avx2_i16(num_seqs, mem->qq, mem->q_lens, mem->max_q_len, 
		    mem->rr, mem->r_lens, mem->max_r_len,
		    match, mismatch, gap_open, gap_extend, 
		    mem->H, mem->F, mem->C, &score[i]);
    //printf("start avx2_matrix (index %i)\n", index);
    //printf("end avx2_matrix\n");
    
#ifdef TIMING
    sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING
    
#ifdef TIMING
    partial_t = sw_tic();
#endif // TIMING
    // backtracking
    backtracking_i16(depth, num_seqs, 
		     mem->q, mem->q_lens, mem->max_q_len, 
		     mem->r, mem->r_lens, mem->max_r_len,
		     gap_open, gap_extend, 
		     mem->H, mem->C, &score[i],
		     &output->query_map_p[i], (int *)&output->query_start_p[i],
		     &output->ref_map_p[i], (int *)&output->ref_start_p[i], 
		     alig_lens,
		     mem->q_aux, mem->r_aux);
#ifdef TIMING
    sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING
  }
  */
}

//-------------------------------------------------------------------------
// sw_mem functions
//-------------------------------------------------------------------------

sw_mem_t *sw_mem_new() {
  const int depth = 32;
  const int align = 32;
  const int len = 250;

  sw_mem_t *mem = (sw_mem_t *) malloc(sizeof(sw_mem_t));
  mem->q_size = len;
  mem->qq = (short *) _mm_malloc(depth * len * sizeof(short), align);

  mem->r_size = len;
  mem->rr = (short *) _mm_malloc(depth * len * sizeof(short), align);

  mem->H_size = len * len;
  mem->H = (short *) _mm_malloc(depth * mem->H_size * sizeof(short), align);
  mem->C = (short *) _mm_malloc(depth * mem->H_size * sizeof(short), align);

  mem->F_size = len;
  mem->F = (short *) _mm_malloc(depth * len * sizeof(short), align);

  mem->aux_size = len;
  mem->q_aux = (char *) calloc(len * 2, sizeof(char));
  mem->r_aux = (char *) calloc(len * 2, sizeof(char));

  for(int i=0; i < depth; i++) {
    mem->q_lens[i] = 0;
    mem->r_lens[i] = 0;
  }

  return mem;
}

//-------------------------------------------------------------------------

void sw_mem_free(sw_mem_t *mem) {
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

void sw_mem_alloc(pair_t *pairs, int num_seqs, int depth, sw_mem_t *mem) {
  const int align = 32;
  const int max_depth = 32;

  printf("begin: mem_alloc depth = %i\n", depth);

  int len, i, j;
  int max_q_len = 0;
  int max_r_len = 0;
  for(int i=0; i < num_seqs; i++) {
    len = pairs[i].query_len;
    mem->q_lens[i] = len;
    if (len > max_q_len) max_q_len = len;
    mem->q[i] = pairs[i].query;

    len = pairs[i].ref_len;
    mem->r_lens[i] = len;
    if (len > max_r_len) max_r_len = len;
    mem->r[i] = pairs[i].ref;
  }

  printf("00\n");

  // check if we have room, otherwise we have to allocate more memory
  mem->max_q_len = max_q_len;
  if (max_q_len > mem->q_size) {
    _mm_free(mem->qq);
    mem->qq = (short *) _mm_malloc(max_depth * max_q_len * sizeof(short), align);
    mem->q_size = max_q_len;
  }

  mem->max_r_len = max_r_len;
  if (max_r_len > mem->r_size) {
    _mm_free(mem->rr);
    mem->rr = (short *) _mm_malloc(max_depth * max_r_len * sizeof(short), align);
    mem->r_size = max_r_len;
  }

  printf("11\n");

  // copy seqs
  if (depth == 32) {
    printf("q 0: %s\n", pairs[0].query);
    printf("r 0: %s\n", pairs[0].ref);
    printf("q 1: %s\n", pairs[1].query);
    printf("r 1: %s\n", pairs[1].ref);
    printf("max q len: %i\n", mem->max_q_len);
    printf("max r len: %i\n", mem->max_r_len);
    char *s = (char *) mem->qq;
    for(i=0; i < depth; i++) {
      for(j=0; j < mem->q_lens[i]; j++) {
	s[depth * j + i] = (char) pairs[i].query[j];
      }
      for(; j < max_q_len; j++) {
	s[depth * j + i] = ';';
      }
      s = (char *) mem->rr;
      for(j=0; j < mem->r_lens[i]; j++) {
        s[depth * j + i] = (char) pairs[i].ref[j];
      }
      for(; j < max_r_len; j++) {
	s[depth * j + i] = ':';
      }
    }

    sw_mem_print_seq("q 0: ", 1, 0, depth, mem);
    sw_mem_print_seq("q 1: ", 1, 1, depth, mem);
    sw_mem_print_seq("r 0: ", 0, 0, depth, mem);
    sw_mem_print_seq("r 1: ", 0, 1, depth, mem);

  } else {
    // depth = 16 (short)
    for(i=0; i < depth; i++) {
      for(j=0; j < mem->q_lens[i]; j++) {
	mem->qq[depth * j + i] = (short) pairs[i].query[j];
      }
      for(; j < max_q_len; j++) {
	mem->qq[depth * j + i] = ';';
      }
      for(j=0; j < mem->r_lens[i]; j++) {
	mem->rr[depth * j + i] = (short) pairs[i].ref[j];
      }
      for(; j < max_r_len; j++) {
	mem->rr[depth * j + i] = ':';
      }
    }
  }

  printf("22\n");

  // check room for H and C matrices
  int matrix_size = max_q_len * max_r_len;
  if (matrix_size > mem->H_size) {
    _mm_free(mem->H);
    _mm_free(mem->C);

    int size = max_depth * matrix_size;
    mem->H = (short *) _mm_malloc(size * sizeof(short), align);
    mem->C = (short *) _mm_malloc(size * sizeof(short), align);
    //printf("new H %x, C %x\n", *H, *C);
    mem->H_size = matrix_size;
  }
  //  memset(*H, 0, size_h);
  //  memset(*C, 0, size_c);

  // check room for F vector
  if (max_q_len > mem->F_size) {
    _mm_free(mem->F);
    mem->F = (short *) _mm_malloc(max_depth * max_q_len * sizeof(short), align);
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

  printf("end: mem_alloc depth = %i\n", depth);
}

//-------------------------------------------------------------------------

void sw_mem_print_seq(char *msg, int query, int index, int depth, sw_mem_t *mem) {
  printf("%s", msg);
  if (query == 1) {
    for (int i = 0; i < mem->max_q_len; i++) {
      printf("%c", ((char *)mem->qq)[SIMD_DEPTH_32 * i + index]);
    }
  } else {
    for (int i = 0; i < mem->max_r_len; i++) {
      printf("%c", ((char *)mem->rr)[SIMD_DEPTH_32 * i + index]);
    }
  }
  printf("\n");
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

