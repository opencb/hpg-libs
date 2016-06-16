#ifdef __AVX512__

#include "avx512.h"

//extern double sse_matrix_t, sse_tracking_t;
//extern double sse1_matrix_t, sse1_tracking_t;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void avx512_matrix2(int num_seqs,
		    char **q, int *q_len, int max_q_len,
		    char **r, int *r_len, int max_r_len,
		    float match, float mismatch, float gap_open, float gap_extend,
		    float *H, float *F, int *C, float *max_score) {

    const int depth = 16;

    int *qq = (int *) _mm_malloc(depth * max_q_len * sizeof(int), 64);
    int *rr = (int *) _mm_malloc(depth * max_r_len * sizeof(int), 64);

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

    __m512i q1, r1;

    __m512 match_simd = _mm512_set1_ps(match);
    __m512 mismatch_simd = _mm512_set1_ps(mismatch);

    __m512 h_simd, e_simd, f_simd, diagonal_simd;
    __m512 temp_simd, subst_simd, subst_simd1;

    __m512i zeroi = _mm512_setzero_si512();

    __m512 score_simd = _mm512_setzero_ps();
    __m512 zero_simd = _mm512_setzero_ps();
    __m512 one_simd = _mm512_set1_ps(1);
    __m512 gap_open_simd = _mm512_set1_ps(gap_open);
    __m512 gap_extend_simd = _mm512_set1_ps(gap_extend);

    __m512 max_de, max_fz;
    __m512 cmp_de, cmp_fz, cmp_de_fz;
    __m512i c;

    __mmask16 mask;

    int offset, idx, j_depth;
    int q_len_depth = depth * max_q_len;

    h_simd = zero_simd;
    e_simd = zero_simd;

    // initial loop
    r1 = _mm512_load_si512((__m512i *) rr);
    for (int j = 0; j < max_q_len; j++) {
        j_depth = depth * j;

        q1 = _mm512_load_si512((__m512i *) (qq + j_depth));

        // left value: gap in reference
        e_simd = _mm512_max_ps(_mm512_sub_ps(e_simd, gap_extend_simd),
                               _mm512_sub_ps(h_simd, gap_open_simd));

        // diagonal value: match or mismatch
        mask = _mm512_cmpeq_epi32_mask(q1, r1);
        subst_simd = _mm512_mask_blend_ps(mask, mismatch_simd, match_simd); //, _mm512_castsi512_ps(_mm512_mask_blend_epi32(mask, q1, r1)));
        diagonal_simd = _mm512_add_ps(zero_simd, subst_simd);

	mask = _mm512_cmp_ps_mask(diagonal_simd, e_simd, _CMP_GE_OQ);
	cmp_de = _mm512_min_ps(_mm512_mask_blend_ps(mask, diagonal_simd, e_simd), one_simd);
        max_de = _mm512_max_ps(diagonal_simd, e_simd);

        // up value: gap in query
        f_simd = _mm512_max_ps(_mm512_sub_ps(zero_simd, gap_extend_simd),
                               _mm512_sub_ps(zero_simd, gap_open_simd));

        mask = _mm512_cmp_ps_mask(f_simd, zero_simd, _CMP_GE_OQ);
        cmp_fz = _mm512_min_ps(_mm512_mask_blend_ps(mask, f_simd, zero_simd), one_simd);
        max_fz = _mm512_max_ps(f_simd, zero_simd);

        // get max. value and save it
        mask = _mm512_cmp_ps_mask(max_de, max_fz, _CMP_GE_OQ);
        cmp_de_fz = _mm512_min_ps(_mm512_mask_blend_ps(mask, max_de, max_fz), one_simd);
        h_simd = _mm512_max_ps(max_de, max_fz);

        score_simd = _mm512_max_ps(score_simd, h_simd);

        // compass (save left, diagonal, up or zero?)
        c = _mm512_slli_epi32(_mm512_or_si512(zeroi, _mm512_cvtps_epi32(cmp_de)), 1);
        c = _mm512_slli_epi32(_mm512_or_si512(c, _mm512_cvtps_epi32(cmp_fz)), 1);
        c = _mm512_or_si512(c, _mm512_cvtps_epi32(cmp_de_fz));
        // update matrices
        _mm512_store_ps(&H[j_depth], h_simd);
        _mm512_store_ps(&F[j_depth], f_simd);
        _mm512_store_si512((__m512i *)&C[j_depth], c);
    }

    // main loop
    for (int i = 1; i < max_r_len; i++) {

        h_simd = zero_simd;
        e_simd = zero_simd;
        temp_simd = zero_simd;

        idx = i * q_len_depth;

        r1 = _mm512_load_si512((__m512i *) (rr + depth * i));
        for (int j = 0; j < max_q_len; j++) {
            j_depth = depth * j;
            offset = idx + j_depth;

            // left value: gap in reference
            e_simd = _mm512_max_ps(_mm512_sub_ps(e_simd, gap_extend_simd),
                                   _mm512_sub_ps(h_simd, gap_open_simd));

            // diagonal value: match or mismatch
            q1 = _mm512_load_si512((__m512i *) (qq + j_depth));
            mask = _mm512_cmpeq_epi32_mask(q1, r1);
	    subst_simd = _mm512_mask_blend_ps(mask, mismatch_simd, match_simd); //, _mm512_castsi512_ps(_mm512_mask_blend_epi32(mask, q1, r1)));
            diagonal_simd = _mm512_add_ps(temp_simd, subst_simd);

            mask = _mm512_cmp_ps_mask(diagonal_simd, e_simd, _CMP_GE_OQ);
            cmp_de = _mm512_min_ps(_mm512_mask_blend_ps(mask, diagonal_simd, e_simd), one_simd);
            max_de = _mm512_max_ps(diagonal_simd, e_simd);

            // up value: gap in query
            temp_simd = _mm512_load_ps(&H[offset - q_len_depth]);

            f_simd = _mm512_load_ps(&F[j_depth]);
            f_simd = _mm512_max_ps(_mm512_sub_ps(f_simd, gap_extend_simd),
                                   _mm512_sub_ps(temp_simd, gap_open_simd));

            mask = _mm512_cmp_ps_mask(f_simd, zero_simd, _CMP_GE_OQ);
            cmp_fz = _mm512_min_ps(_mm512_mask_blend_ps(mask, f_simd, zero_simd), one_simd);
            max_fz = _mm512_max_ps(f_simd, zero_simd);

            // get max. value
            mask = _mm512_cmp_ps_mask(max_de, max_fz, _CMP_GE_OQ);
            cmp_de_fz = _mm512_min_ps(_mm512_mask_blend_ps(mask, max_de, max_fz), one_simd);
            h_simd = _mm512_max_ps(max_de, max_fz);

            score_simd = _mm512_max_ps(score_simd, h_simd);

            // compass (save left, diagonal, up or zero?)
            c = _mm512_slli_epi32(_mm512_or_si512(zeroi, _mm512_cvtps_epi32(cmp_de)), 1);
            c = _mm512_slli_epi32(_mm512_or_si512(c, _mm512_cvtps_epi32(cmp_fz)), 1);
            c = _mm512_or_si512(c, _mm512_cvtps_epi32(cmp_de_fz));

            // update matrices
            _mm512_store_ps(&H[offset], h_simd);
            _mm512_store_ps(&F[j_depth], f_simd);
            _mm512_store_si512((__m512i *)&C[offset], c);

        }
    }

    // return back max scores
    _mm512_store_ps(max_score, score_simd);

    // free memory
    _mm_free(qq);
    _mm_free(rr);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
#endif // __AVX512__
