#include "sw_f32.h"

//extern double sse_matrix_t, sse_tracking_t;
//extern double sse1_matrix_t, sse1_tracking_t;

sw_mem_f32_t *sw_mems[1024];

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void sw_mqmr_f32(char **query, char **ref, unsigned int num_queries,
		 sw_optarg_t *optarg, sw_multi_output_t *output,
		 sw_mem_f32_t *mem) {

  if (output == NULL) {
    printf("Error: output buffer is null.\n");
    exit(-1);
  }
  if (output->num_queries < num_queries) {
    printf("Error: num. sequencias (%i) and output buffer length (%i) mismatch !\n",
	   num_queries, output->num_queries);
    exit(-1);
  }
  
  int tid = omp_get_thread_num();

#ifdef TIMING
  double partial_t;
#endif // TIMING

  //  char *q_aux = NULL, *r_aux = NULL;
  //  int depth, aux_size = 0, H_size = 0, F_size = 0, max_q_len = 0, max_r_len = 0;
  //  char *q[SIMD_DEPTH], *r[SIMD_DEPTH];
  //  int len, index, q_lens[SIMD_DEPTH], r_lens[SIMD_DEPTH], alig_lens[SIMD_DEPTH];
  //  float *H = NULL, *F = NULL;
  //  int *C = NULL;
  int alig_lens[SIMD_DEPTH];
  
  float match = optarg->match;
  float mismatch = optarg->mismatch;
  float gap_open = optarg->gap_open;
  float gap_extend = optarg->gap_extend;
  float *score = output->score_p;
  
  //printf("num queries = %i\n", num_queries);
  //for(int k = 0; k < 8; k++) printf("%0.2f ", score_p[k]);
  //printf("\n");
  
  int num_seqs;
  for(int i = 0; i < num_queries; i += SIMD_DEPTH) {
    num_seqs = num_queries - i;
    if (num_seqs > SIMD_DEPTH) num_seqs = SIMD_DEPTH; 
    //    printf("after sw_mem_f32_new(): num_seqs = %i\n", num_seqs);
    sw_mem_f32_alloc(&query[i], &ref[i], num_seqs, mem);

    //    printf("-----> after alloc\n");

    // generating score matrix
#ifdef TIMING
    partial_t = sw_tic();
#endif // TIMING
      
#ifdef __AVX2__
    matrix_avx2_f32(num_seqs, mem->qq, mem->q_lens, mem->max_q_len, 
		    mem->rr, mem->r_lens, mem->max_r_len,
		    match, mismatch, gap_open, gap_extend, 
		    mem->H, mem->F, mem->C, &score[i]);
#else
    printf("Only 32-bit float AVX2 is implemented!");
    exit(-1);
#endif // __AVX2__
      
    //printf("start avx2_matrix (index %i)\n", index);
    //printf("end avx2_matrix\n");
    
#ifdef TIMING
    sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING
    
#ifdef TIMING
    partial_t = sw_tic();
#endif // TIMING
    // backtracking
    backtracking_f32(SIMD_DEPTH, num_seqs, 
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
}

//------------------------------------------------------------------------------------

void matrix_avx2_f32(int num_seqs,
		     int *qq, int *q_len, int max_q_len,
		     int *rr, int *r_len, int max_r_len,
		     float match, float mismatch, float gap_open, float gap_extend,
		     float *H, float *F, int *C, float *max_score) {
    __m256i q1, r1;

    __m256 match_simd = _mm256_set1_ps(match);
    __m256 mismatch_simd = _mm256_set1_ps(mismatch);

    __m256 h_simd, e_simd, f_simd, diagonal_simd;
    __m256 temp_simd, subst_simd, subst_simd1;

    __m256i zeroi = _mm256_setzero_si256();

    __m256 score_simd = _mm256_setzero_ps();
    __m256 zero_simd = _mm256_setzero_ps();
    __m256 one_simd = _mm256_set1_ps(1);
    __m256 gap_open_simd = _mm256_set1_ps(gap_open);
    __m256 gap_extend_simd = _mm256_set1_ps(gap_extend);

    __m256 max_de, max_fz;
    __m256 cmp_de, cmp_fz, cmp_de_fz;
    __m256i c;

    int offset, idx, j_depth;
    int q_len_depth = SIMD_DEPTH * max_q_len;

    h_simd = zero_simd;
    e_simd = zero_simd;

    // initial loop
    r1 = _mm256_load_si256((__m256i *) rr);
    for (int j = 0; j < max_q_len; j++) {
        j_depth = SIMD_DEPTH * j;

        q1 = _mm256_load_si256((__m256i *) (qq + j_depth));

        // left value: gap in reference
        e_simd = _mm256_max_ps(_mm256_sub_ps(e_simd, gap_extend_simd),
                               _mm256_sub_ps(h_simd, gap_open_simd));

        // diagonal value: match or mismatch
        subst_simd = _mm256_blendv_ps(mismatch_simd, match_simd, _mm256_castsi256_ps(_mm256_cmpeq_epi32(q1, r1)));
        diagonal_simd = _mm256_add_ps(zero_simd, subst_simd);

        cmp_de = _mm256_min_ps(_mm256_cmp_ps(diagonal_simd, e_simd, _CMP_GE_OQ), one_simd);
        max_de = _mm256_max_ps(diagonal_simd, e_simd);

        // up value: gap in query
        f_simd = _mm256_max_ps(_mm256_sub_ps(zero_simd, gap_extend_simd),
                               _mm256_sub_ps(zero_simd, gap_open_simd));

        cmp_fz = _mm256_min_ps(_mm256_cmp_ps(f_simd, zero_simd, _CMP_GE_OQ), one_simd);
        max_fz = _mm256_max_ps(f_simd, zero_simd);

        // get max. value and save it
        cmp_de_fz = _mm256_min_ps(_mm256_cmp_ps(max_de, max_fz, _CMP_GE_OQ), one_simd);
        h_simd = _mm256_max_ps(max_de, max_fz);

        score_simd = _mm256_max_ps(score_simd, h_simd);

        // compass (save left, diagonal, up or zero?)
        c = _mm256_slli_epi32(_mm256_or_si256(zeroi, _mm256_cvtps_epi32(cmp_de)), 1);
        c = _mm256_slli_epi32(_mm256_or_si256(c, _mm256_cvtps_epi32(cmp_fz)), 1);
        c = _mm256_or_si256(c, _mm256_cvtps_epi32(cmp_de_fz));

        // update matrices
        _mm256_store_ps(&H[j_depth], h_simd);
        _mm256_store_ps(&F[j_depth], f_simd);
        _mm256_store_si256((__m256i *)&C[j_depth], c);
    }

    // main loop
    for (int i = 1; i < max_r_len; i++) {

        h_simd = zero_simd;
        e_simd = zero_simd;
        temp_simd = zero_simd;

        idx = i * q_len_depth;

        r1 = _mm256_load_si256((__m256i *) (rr + SIMD_DEPTH * i));
        for (int j = 0; j < max_q_len; j++) {
            j_depth = SIMD_DEPTH * j;
            offset = idx + j_depth;

            // left value: gap in reference
            e_simd = _mm256_max_ps(_mm256_sub_ps(e_simd, gap_extend_simd),
                                   _mm256_sub_ps(h_simd, gap_open_simd));

            // diagonal value: match or mismatch
            q1 = _mm256_load_si256((__m256i *) (qq + j_depth));
            subst_simd = _mm256_blendv_ps(mismatch_simd, match_simd, _mm256_castsi256_ps(_mm256_cmpeq_epi32(q1, r1)));
            diagonal_simd = _mm256_add_ps(temp_simd, subst_simd);

            cmp_de = _mm256_min_ps(_mm256_cmp_ps(diagonal_simd, e_simd, _CMP_GE_OQ), one_simd);
            max_de = _mm256_max_ps(diagonal_simd, e_simd);

            // up value: gap in query
            temp_simd = _mm256_load_ps(&H[offset - q_len_depth]);

            f_simd = _mm256_load_ps(&F[j_depth]);
            f_simd = _mm256_max_ps(_mm256_sub_ps(f_simd, gap_extend_simd),
                                   _mm256_sub_ps(temp_simd, gap_open_simd));

            cmp_fz = _mm256_min_ps(_mm256_cmp_ps(f_simd, zero_simd, _CMP_GE_OQ), one_simd);
            max_fz = _mm256_max_ps(f_simd, zero_simd);

            // get max. value
            cmp_de_fz = _mm256_min_ps(_mm256_cmp_ps(max_de, max_fz, _CMP_GE_OQ), one_simd);
            h_simd = _mm256_max_ps(max_de, max_fz);

            score_simd = _mm256_max_ps(score_simd, h_simd);

            // compass (save left, diagonal, up or zero?)
            c = _mm256_slli_epi32(_mm256_or_si256(zeroi, _mm256_cvtps_epi32(cmp_de)), 1);
            c = _mm256_slli_epi32(_mm256_or_si256(c, _mm256_cvtps_epi32(cmp_fz)), 1);
            c = _mm256_or_si256(c, _mm256_cvtps_epi32(cmp_de_fz));

            // update matrices
            _mm256_store_ps(&H[offset], h_simd);
            _mm256_store_ps(&F[j_depth], f_simd);
            _mm256_store_si256((__m256i *)&C[offset], c);

        }
    }

    // return back max scores
    _mm256_store_ps(max_score, score_simd);
}

//------------------------------------------------------------------------

void backtracking_f32(int depth, int num_seqs,
		      char **q, int *q_len, int max_q_len,
		      char **r, int *r_len, int max_r_len,
		      float gap_open, float gap_extend,
		      float *H, int *C,
		      float *max_score,
		      char **q_alig, int *q_start,
		      char **r_alig, int *r_start,
		      int *len_alig,
		      char *q_aux, char *r_aux) {

    char *qq, *qq_alig;
    char *rr, *rr_alig;

    int jj, kk;

    int ix, iy, gapcnt, index;
    double bimble;
    double errbounds = (double) 0.01;

    int qq_len, rr_len;
    int len = 0;
    char c;
    float score;

    /*
    //max_len = max_q_len * max_r_len;
    max_len = 2 * max_r_len;
    printf("allocating max_len %i...\n", max_len);
    char *q_aux = (char *) calloc(max_len, sizeof(char));
    char *r_aux = (char *) calloc(max_len, sizeof(char));
    printf("done\n");
    */

    for (int i = 0; i < num_seqs ; i++) {
        qq = q[i];
        rr = r[i];
        qq_len = q_len[i];
        rr_len = r_len[i];

        //    printf("%i of %i\n", i, num_seqs);
        find_position_f32(depth, i, qq, qq_len, rr, rr_len, H, max_q_len, max_r_len, max_score[i], &kk, &jj);
        //printf("index %i: kk = %i, jj = %i\n", i, kk, jj);

        len = 0;
        //    while ((c = C[(jj * max_q_len * 4) + (kk * 4) + i])) {
        while (1) {
            index = (jj * max_q_len * depth) + (kk * depth) + i;
            score = H[index];
            if (score == 0.0f) {
                break;
            }

            c = C[index];

            if (c == 5 || c == 7) { // diagonal
                q_aux[len] = qq[kk];
                r_aux[len] = rr[jj];
                len++;
                jj--;
                kk--;

                if (jj < 0) break;
                if (kk < 0) break;

                continue;
            } else if (c == 1 || c == 3) { // left
                gapcnt = 0;
                ix     = kk - 1;

                while (1) {
                    bimble = H[(jj * max_q_len * depth) + (ix * depth) + i] - gap_open - (gapcnt * gap_extend);

                    if (!ix || fabs((double)score - (double)bimble) < errbounds)
                        break;

                    --ix;
                    if (ix < 0) {
                        printf("Error reconstructing alignment (walking left)\n");
                        exit(-1);
                    }
                    ++gapcnt;
                }

                if(bimble <= 0.0)
                    break;

                for(int ic = 0; ic <= gapcnt; ++ic) {
                    /*if (len > (max_r_len * 2)) {
                      printf("SEQ(%i):%s\n", qq_len, qq);
                      printf("REF(%i):%s\n", rr_len, rr);
                      printf("%i > %i\n", len, max_r_len*2 );
                      exit(-1);
                      }*/
                    q_aux[len] = qq[kk--];
                    r_aux[len] = '-';
                    len++;
                    //if (jj < 0) break;
                }

                continue;
            } else if (c == 2 || c == 6) { // down
                gapcnt = 0;
                iy     = jj - 1;

                while (1) {
                    bimble = H[(iy * max_q_len * depth) + (kk * depth) + i] - gap_open - (gapcnt * gap_extend);
                    //printf("bimble = %0.2f\n", bimble);

                    if (!iy || fabs((double)score - (double)bimble) < errbounds)
                        break;

                    --iy;
                    if (iy < 0) {
                        printf("Error reconstructing alignment (walking down)\n");
                        exit(-1);
                    }
                    ++gapcnt;
                }

                if(bimble <= 0.0)
                    break;

                //printf("\tDOWN: gapcnt = %i\n", gapcnt);
                for(int ic = 0; ic <= gapcnt; ++ic) {
                    q_aux[len] = '-';
                    r_aux[len] = rr[jj--];
                    len++;
                    //printf("\tlen = %i, q_aux: %s\n", len, q_aux);
                    //printf("\tr_aux:%s\n", r_aux);
                    //if (kk < 0) break;
                }

                continue;
            } else {
                break;
                //printf("Error reconstructing alignment\n");
                //exit(-1);
            }
            //printf("pos (j,k) = (%i, %i), len = %i\n", jj, kk, len);
        }


        q_start[i] = kk + 1;
        r_start[i] = jj + 1;
        //q_start[i] = (kk < 0 ? 1 : kk + 1);
        //r_start[i] = (jj < 0 ? 1 : jj + 1);

        q_aux[len] = 0;
        r_aux[len] = 0;

        //    printf("index %i\tq = %s\n", i, q_aux);
        //    printf("index %i\tr = %s\n", i, r_aux);

        qq_alig = (char *) calloc(len + 1, sizeof(char));
        rr_alig = (char *) calloc(len + 1, sizeof(char));
        jj = len - 1;

        for (int i = 0; i < len; i++) {
            qq_alig[i] = q_aux[jj];
            rr_alig[i] = r_aux[jj];
            jj--;
        }
        qq_alig[len] = 0;
        rr_alig[len] = 0;

        q_alig[i] = qq_alig;
        r_alig[i] = rr_alig;

        len_alig[i] = len;
        //    printf("rev. index %i\tq = %s\n", i, qq_alig);
        //    printf("rev. index %i\tr = %s\n", i, rr_alig);
    }
    /*
    printf("free %i...\n", max_len);
    free(q_aux);
    free(r_aux);
    printf("done\n");
    //  printf("%s\n%s\n", q_alig, r_alig);
    */
}

//-------------------------------------------------------------------------

void find_position_f32(int depth, int index, char *q, int q_len, char *r, int r_len,
		       float *H, int cols, int rows, float score,
		       int *q_pos, int *r_pos) {
    *r_pos = 0;
    *q_pos = 0;

    //  score = 473;
    //  score = H[100];
    //  printf("max. score = %0.2f\n", score);
    //  printf("\tcols = %i, rows = %i, depth = %i, index = %i, r_len = %i, q_len = %i\n",
    //  	 score, cols, rows, depth, index, r_len, q_len);
    int i, ii = 0;
    for (int j = 0; j < r_len; j++) {
        i = j * cols * depth;
        for (int k = 0; k < q_len; k++) {
            ii = i + (k * depth) + index;
            //      printf("\tii (%i)= i (%i) + [k  (%i) * depth (%i)] + index (%i)\n", ii, i, k, depth, index);
            if (H[ii] == score) {
                *r_pos = j;
                *q_pos = k;
                return;
            }
        }
    }
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

