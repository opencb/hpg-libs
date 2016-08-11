#ifdef __AVX2__

#include "avx2_i32.h"

//extern double sse_matrix_t, sse_tracking_t;
//extern double sse1_matrix_t, sse1_tracking_t;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void avx2_matrix2_i32(int num_seqs,
		      char **q, int *q_len, int max_q_len,
		      char **r, int *r_len, int max_r_len,
		      int match, int mismatch, int gap_open, int gap_extend,
		      int *H, int *F, int *C, int *max_score) {

    const int depth = 8;

    int *qq = (int *) _mm_malloc(depth * max_q_len * sizeof(int), 32);
    int *rr = (int *) _mm_malloc(depth * max_r_len * sizeof(int), 32);

    int jj;
    for(int ii=0; ii < depth; ii++) {
        for(jj=0; jj < q_len[ii]; jj++) {
            qq[depth * jj + ii] = (int) q[ii][jj];
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

    __m256i q1, r1;

    __m256i match_simd = _mm256_set1_epi32(match);
    __m256i mismatch_simd = _mm256_set1_epi32(mismatch);

    __m256i h_simd, e_simd, f_simd, diagonal_simd;
    __m256i temp_simd, subst_simd, subst_simd1;

    __m256i zeroi = _mm256_setzero_si256();

    __m256i score_simd = _mm256_setzero_si256();
    __m256i zero_simd = _mm256_setzero_si256();
    __m256i one_simd = _mm256_set1_epi32(1);
    __m256i gap_open_simd = _mm256_set1_epi32(gap_open);
    __m256i gap_extend_simd = _mm256_set1_epi32(gap_extend);

    __m256i max_de, max_fz;
    __m256i cmp_de, cmp_fz, cmp_de_fz;
    __m256i c;

    int offset, idx, j_depth;
    int q_len_depth = depth * max_q_len;

    h_simd = zero_simd;
    e_simd = zero_simd;

    // initial loop
    r1 = _mm256_load_si256((__m256i *) rr);
    for (int j = 0; j < max_q_len; j++) {
        j_depth = depth * j;

        q1 = _mm256_load_si256((__m256i *) (qq + j_depth));

        // left value: gap in reference
        e_simd = _mm256_max_epi32(_mm256_sub_epi32(e_simd, gap_extend_simd),
				  _mm256_sub_epi32(h_simd, gap_open_simd));

        // diagonal value: match or mismatch
	//        subst_simd = _mm256_blend_epi32(mismatch_simd, match_simd, _mm256_cmpeq_epi32(q1, r1));
        diagonal_simd = _mm256_add_epi32(zero_simd, subst_simd);

	//        cmp_de = _mm256_min_epi32(_mm256_cmp_epi32(diagonal_simd, e_simd, _CMP_GE_OQ), one_simd);
        max_de = _mm256_max_epi32(diagonal_simd, e_simd);

        // up value: gap in query
        f_simd = _mm256_max_epi32(_mm256_sub_epi32(zero_simd, gap_extend_simd),
				  _mm256_sub_epi32(zero_simd, gap_open_simd));

	cmp_fz = _mm256_min_epi32(_mm256_cmpge_epi32(f_simd, zero_simd), one_simd);
        max_fz = _mm256_max_epi32(f_simd, zero_simd);

        // get max. value and save it
	//        cmp_de_fz = _mm256_min_epi32(_mm256_cmp_epi32(max_de, max_fz, _CMP_GE_OQ), one_simd);
        h_simd = _mm256_max_epi32(max_de, max_fz);

        score_simd = _mm256_max_epi32(score_simd, h_simd);

        // compass (save left, diagonal, up or zero?)
        c = _mm256_slli_epi32(_mm256_or_si256(zeroi, cmp_de), 1);
        c = _mm256_slli_epi32(_mm256_or_si256(c, cmp_fz), 1);
        c = _mm256_or_si256(c, cmp_de_fz);

        // update matrices
        _mm256_store_si256((__m256i *) &H[j_depth], h_simd);
        _mm256_store_si256((__m256i *) &F[j_depth], f_simd);
        _mm256_store_si256((__m256i *) &C[j_depth], c);
    }

    // main loop
    for (int i = 1; i < max_r_len; i++) {

        h_simd = zero_simd;
        e_simd = zero_simd;
        temp_simd = zero_simd;

        idx = i * q_len_depth;

        r1 = _mm256_load_si256((__m256i *) (rr + depth * i));
        for (int j = 0; j < max_q_len; j++) {
            j_depth = depth * j;
            offset = idx + j_depth;

            // left value: gap in reference
            e_simd = _mm256_max_epi32(_mm256_sub_epi32(e_simd, gap_extend_simd),
				      _mm256_sub_epi32(h_simd, gap_open_simd));

            // diagonal value: match or mismatch
            q1 = _mm256_load_si256((__m256i *) (qq + j_depth));
	    //            subst_simd = _mm256_blend_epi32(mismatch_simd, match_simd, _mm256_cmpeq_epi32(q1, r1));
            diagonal_simd = _mm256_add_epi32(temp_simd, subst_simd);

	    //            cmp_de = _mm256_min_epi32(_mm256_cmp_epi32(diagonal_simd, e_simd, _CMP_GE_OQ), one_simd);
            max_de = _mm256_max_epi32(diagonal_simd, e_simd);

            // up value: gap in query
            temp_simd = _mm256_load_si256((__m256i *) &H[offset - q_len_depth]);
	    
            f_simd = _mm256_load_si256((__m256i *) &F[j_depth]);
            f_simd = _mm256_max_epi32(_mm256_sub_epi32(f_simd, gap_extend_simd),
				      _mm256_sub_epi32(temp_simd, gap_open_simd));

	    //            cmp_fz = _mm256_min_epi32(_mm256_cmp_epi32(f_simd, zero_simd, _CMP_GE_OQ), one_simd);
            max_fz = _mm256_max_epi32(f_simd, zero_simd);

            // get max. value
	    //            cmp_de_fz = _mm256_min_epi32(_mm256_cmp_epi32(max_de, max_fz, _CMP_GE_OQ), one_simd);
            h_simd = _mm256_max_epi32(max_de, max_fz);

            score_simd = _mm256_max_epi32(score_simd, h_simd);

            // compass (save left, diagonal, up or zero?)
            c = _mm256_slli_epi32(_mm256_or_si256(zeroi, cmp_de), 1);
            c = _mm256_slli_epi32(_mm256_or_si256(c, cmp_fz), 1);
            c = _mm256_or_si256(c, cmp_de_fz);

            // update matrices
            _mm256_store_si256((__m256i *) &H[offset], h_simd);
            _mm256_store_si256((__m256i *) &F[j_depth], f_simd);
            _mm256_store_si256((__m256i *) &C[offset], c);

        }
    }

    // return back max scores
    _mm256_store_si256((__m256i *) max_score, score_simd);

    // free memory
    _mm_free(qq);
    _mm_free(rr);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
#endif // __AVX2__