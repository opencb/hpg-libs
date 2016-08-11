#ifdef __AVX2_I16__
#include "avx2_i16.h"

//-------------------------------------------------------------------------

void avx2_matrix_i16(int num_seqs,
		     short *qq, int max_q_len,
		     short *rr, int max_r_len,
		     short match, short mismatch, short gap_open, short gap_extend,
		     short *H, short *F, short *C, short *max_score) {

    const int depth = 16;

    __m256i q1, r1;

    __m256i match_simd = _mm256_set1_epi16(match);
    __m256i mismatch_simd = _mm256_set1_epi16(mismatch);

    __m256i h_simd, e_simd, f_simd, diagonal_simd, mask;
    __m256i temp_simd, subst_simd, subst_simd1;

    __m256i zeroi = _mm256_setzero_si256();

    __m256i score_simd = _mm256_setzero_si256();
    __m256i zero_simd = _mm256_setzero_si256();
    __m256i one_simd = _mm256_set1_epi16(1);
    __m256i gap_open_simd = _mm256_set1_epi16(gap_open);
    __m256i gap_extend_simd = _mm256_set1_epi16(gap_extend);

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
        e_simd = _mm256_max_epi16(_mm256_sub_epi16(e_simd, gap_extend_simd),
				  _mm256_sub_epi16(h_simd, gap_open_simd));

        // diagonal value: match or mismatch
	mask = _mm256_cmpeq_epi16(q1, r1);
        subst_simd = _mm256_or_si256(_mm256_andnot_si256(mask, mismatch_simd), _mm256_and_si256(mask, match_simd));
        diagonal_simd = _mm256_add_epi16(zero_simd, subst_simd);

        cmp_de = _mm256_min_epi16(_mm256_cmp_epi16_mask(diagonal_simd, e_simd, _CMP_GE_OQ), one_simd);
        max_de = _mm256_max_epi16(diagonal_simd, e_simd);

        // up value: gap in query
        f_simd = _mm256_max_epi16(_mm256_sub_epi16(zero_simd, gap_extend_simd),
				  _mm256_sub_epi16(zero_simd, gap_open_simd));

        cmp_fz = _mm256_min_epi16(_mm256_cmp_epi16_mask(f_simd, zero_simd, _CMP_GE_OQ), one_simd);
        max_fz = _mm256_max_epi16(f_simd, zero_simd);

        // get max. value and save it
        cmp_de_fz = _mm256_min_epi16(_mm256_cmp_epi16_mask(max_de, max_fz, _CMP_GE_OQ), one_simd);
        h_simd = _mm256_max_epi16(max_de, max_fz);

        score_simd = _mm256_max_epi16(score_simd, h_simd);

        // compass (save left, diagonal, up or zero?)
        c = _mm256_slli_epi16(_mm256_or_si256(zeroi, cmp_de), 1);
        c = _mm256_slli_epi16(_mm256_or_si256(c, cmp_fz), 1);
        c = _mm256_or_si256(c, cmp_de_fz);

        // update matrices
        _mm256_store_si256(&H[j_depth], h_simd);
        _mm256_store_si256(&F[j_depth], f_simd);
        _mm256_store_si256(&C[j_depth], c);
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
            e_simd = _mm256_max_epi16(_mm256_sub_epi16(e_simd, gap_extend_simd),
				      _mm256_sub_epi16(h_simd, gap_open_simd));

            // diagonal value: match or mismatch
            q1 = _mm256_load_si256((__m256i *) (qq + j_depth));
	    mask = _mm256_cmpeq_epi16(q1, r1);
	    subst_simd = _mm256_or_si256(_mm256_andnot_si256(mask, mismatch_simd), _mm256_and_si256(mask, match_simd));
            diagonal_simd = _mm256_add_epi16(temp_simd, subst_simd);

            cmp_de = _mm256_min_epi16(_mm256_cmp_epi16_mask(diagonal_simd, e_simd, _CMP_GE_OQ), one_simd);
            max_de = _mm256_max_epi16(diagonal_simd, e_simd);

            // up value: gap in query
            temp_simd = _mm256_load_epi16(&H[offset - q_len_depth]);

            f_simd = _mm256_load_epi16(&F[j_depth]);
            f_simd = _mm256_max_epi16(_mm256_sub_epi16(f_simd, gap_extend_simd),
				      _mm256_sub_epi16(temp_simd, gap_open_simd));

            cmp_fz = _mm256_min_epi16(_mm256_cmp_epi16_mask(f_simd, zero_simd, _CMP_GE_OQ), one_simd);
            max_fz = _mm256_max_epi16(f_simd, zero_simd);

            // get max. value
            cmp_de_fz = _mm256_min_epi16(_mm256_cmp_epi16_mask(max_de, max_fz, _CMP_GE_OQ), one_simd);
            h_simd = _mm256_max_epi16(max_de, max_fz);

            score_simd = _mm256_max_epi16(score_simd, h_simd);

            // compass (save left, diagonal, up or zero?)
            c = _mm256_slli_epi16(_mm256_or_si256(zeroi, cmp_de), 1);
            c = _mm256_slli_epi16(_mm256_or_si256(c, cmp_fz), 1);
            c = _mm256_or_si256(c, cmp_de_fz);

            // update matrices
            _mm256_store_si256(&H[offset], h_simd);
            _mm256_store_si256(&F[j_depth], f_simd);
            _mm256_store_si256(&C[offset], c);

        }
    }

    // return back max scores
    _mm256_store_si256(max_score, score_simd);

}

//------------------------------------------------------------------------

void simd_traceback_i16(int depth, int num_seqs,
			char **q, int *q_len, int max_q_len,
			char **r, int *r_len, int max_r_len,
			short gap_open, short gap_extend,
			short *H, short *C,
			short *max_score,
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
        simd_find_position_i16(depth, i, qq, qq_len, rr, rr_len, H, max_q_len, max_r_len, max_score[i], &kk, &jj);
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

void simd_find_position_i16(int depth, int index, char *q, int q_len, char *r, int r_len,
			    short *H, int cols, int rows, short score,
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
#endif // __AVX2_I16__

