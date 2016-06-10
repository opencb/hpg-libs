#include "avx2.h"

//extern double sse_matrix_t, sse_tracking_t;
//extern double sse1_matrix_t, sse1_tracking_t;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void avx2_matrix(int num_seqs,
                 char **q, int *q_len, int max_q_len,
                 char **r, int *r_len, int max_r_len,
                 float profile[128][128], float gap_open, float gap_extend,
                 float *H, float *F, int *C, float *max_score) {

    const int depth = 8;

    __m256 h_simd, e_simd, f_simd, diagonal_simd;
    __m256 temp_simd, subst_simd;

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
    int q_len_depth = depth * max_q_len;

    h_simd = zero_simd;
    e_simd = zero_simd;

    // initial loop
    for (int j = 0; j < max_q_len; j++) {

        j_depth = depth * j;

        // left value: gap in reference
        e_simd = _mm256_max_ps(_mm256_sub_ps(e_simd, gap_extend_simd),
                               _mm256_sub_ps(h_simd, gap_open_simd));

        // diagonal value: match or mismatch
        subst_simd = _mm256_set_ps((q_len[7] > j) ? profile[(unsigned char)q[7][j]][(unsigned char)r[7][0]] : -1000.0f,
                                   (q_len[6] > j) ? profile[(unsigned char)q[6][j]][(unsigned char)r[6][0]] : -1000.0f,
                                   (q_len[5] > j) ? profile[(unsigned char)q[5][j]][(unsigned char)r[5][0]] : -1000.0f,
                                   (q_len[4] > j) ? profile[(unsigned char)q[4][j]][(unsigned char)r[4][0]] : -1000.0f,
                                   (q_len[3] > j) ? profile[(unsigned char)q[3][j]][(unsigned char)r[3][0]] : -1000.0f,
                                   (q_len[2] > j) ? profile[(unsigned char)q[2][j]][(unsigned char)r[2][0]] : -1000.0f,
                                   (q_len[1] > j) ? profile[(unsigned char)q[1][j]][(unsigned char)r[1][0]] : -1000.0f,
                                   (q_len[0] > j) ? profile[(unsigned char)q[0][j]][(unsigned char)r[0][0]] : -1000.0f);

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

        for (int j = 0; j < max_q_len; j++) {
            j_depth = depth * j;
            offset = idx + j_depth;

            // left value: gap in reference
            e_simd = _mm256_max_ps(_mm256_sub_ps(e_simd, gap_extend_simd),
                                   _mm256_sub_ps(h_simd, gap_open_simd));

            // diagonal value: match or mismatch
            diagonal_simd = _mm256_add_ps(temp_simd,
                                          _mm256_set_ps((q_len[7] > j && r_len[7] > i) ? profile[(unsigned char)q[7][j]][(unsigned char)r[7][i]] : -1000.0f,
                                                        (q_len[6] > j && r_len[6] > i) ? profile[(unsigned char)q[6][j]][(unsigned char)r[6][i]] : -1000.0f,
                                                        (q_len[5] > j && r_len[5] > i) ? profile[(unsigned char)q[5][j]][(unsigned char)r[5][i]] : -1000.0f,
                                                        (q_len[4] > j && r_len[4] > i) ? profile[(unsigned char)q[4][j]][(unsigned char)r[4][i]] : -1000.0f,
                                                        (q_len[3] > j && r_len[3] > i) ? profile[(unsigned char)q[3][j]][(unsigned char)r[3][i]] : -1000.0f,
                                                        (q_len[2] > j && r_len[2] > i) ? profile[(unsigned char)q[2][j]][(unsigned char)r[2][i]] : -1000.0f,
                                                        (q_len[1] > j && r_len[1] > i) ? profile[(unsigned char)q[1][j]][(unsigned char)r[1][i]] : -1000.0f,
                                                        (q_len[0] > j && r_len[0] > i) ? profile[(unsigned char)q[0][j]][(unsigned char)r[0][i]] : -1000.0f)
            );

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

//-------------------------------------------------------------------------

void avx2_matrix2(int num_seqs,
                  char **q, int *q_len, int max_q_len,
                  char **r, int *r_len, int max_r_len,
                  float match, float mismatch, float gap_open, float gap_extend,
                  float *H, float *F, int *C, float *max_score) {

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
    int q_len_depth = depth * max_q_len;

    h_simd = zero_simd;
    e_simd = zero_simd;

    // initial loop
    r1 = _mm256_load_si256((__m256i *) rr);
    for (int j = 0; j < max_q_len; j++) {
        j_depth = depth * j;

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

        r1 = _mm256_load_si256((__m256i *) (rr + depth * i));
        for (int j = 0; j < max_q_len; j++) {
            j_depth = depth * j;
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

    // free memory
    _mm_free(qq);
    _mm_free(rr);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
