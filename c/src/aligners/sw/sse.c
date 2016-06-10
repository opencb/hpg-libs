#include "sse.h"

extern double sse_matrix_t, sse_tracking_t;
extern double sse1_matrix_t, sse1_tracking_t;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void sse_matrix(int num_seqs, 
		char **q, int *q_len, int max_q_len,
		char **r, int *r_len, int max_r_len,
		float profile[128][128], float gap_open, float gap_extend,
		float *H, float *F, int *C, float *max_score) {
  
  const int depth = 4;

  __m128 h_simd, e_simd, f_simd, diagonal_simd;
  __m128 temp_simd, subst_simd;

  __m128i zeroi   = _mm_setzero_si128();

  __m128 score_simd = _mm_setzero_ps();
  __m128 zero_simd = _mm_setzero_ps();
  __m128 one_simd = _mm_set1_ps(1);
  __m128 gap_open_simd = _mm_set1_ps(gap_open);
  __m128 gap_extend_simd = _mm_set1_ps(gap_extend);

  __m128 max_de, max_fz;
  __m128 cmp_de, cmp_fz, cmp_de_fz;
  __m128i c;

  int offset, idx, j_depth;
  int q_len_depth = depth * max_q_len;

  h_simd = zero_simd;
  e_simd = zero_simd;

  // initial loop
  for (int j = 0; j < max_q_len; j++) {

    j_depth = depth * j;
    
    // left value: gap in reference
    e_simd = _mm_max_ps(_mm_sub_ps(e_simd, gap_extend_simd), 
			_mm_sub_ps(h_simd, gap_open_simd));


    // diagonal value: match or mismatch
    subst_simd = _mm_set_ps((q_len[3] > j) ? profile[(unsigned char)q[3][j]][(unsigned char)r[3][0]] : -1000.0f,
                            (q_len[2] > j) ? profile[(unsigned char)q[2][j]][(unsigned char)r[2][0]] : -1000.0f,
                            (q_len[1] > j) ? profile[(unsigned char)q[1][j]][(unsigned char)r[1][0]] : -1000.0f,
                            (q_len[0] > j) ? profile[(unsigned char)q[0][j]][(unsigned char)r[0][0]] : -1000.0f);

    diagonal_simd = _mm_add_ps(zero_simd, subst_simd);

    cmp_de = _mm_min_ps(_mm_cmpge_ps(diagonal_simd, e_simd), one_simd);
    max_de = _mm_max_ps(diagonal_simd, e_simd);

    // up value: gap in query
    f_simd = _mm_max_ps(_mm_sub_ps(zero_simd, gap_extend_simd), 
			_mm_sub_ps(zero_simd, gap_open_simd));

    cmp_fz = _mm_min_ps(_mm_cmpge_ps(f_simd, zero_simd), one_simd);
    max_fz = _mm_max_ps(f_simd, zero_simd);
    
    // get max. value and save it
    cmp_de_fz = _mm_min_ps(_mm_cmpge_ps(max_de, max_fz), one_simd);
    h_simd = _mm_max_ps(max_de, max_fz);

    score_simd = _mm_max_ps(score_simd, h_simd);

    // compass (save left, diagonal, up or zero?)
    c = _mm_slli_epi32(_mm_or_si128(zeroi, _mm_cvtps_epi32(cmp_de)), 1);
    c = _mm_slli_epi32(_mm_or_si128(c, _mm_cvtps_epi32(cmp_fz)), 1);
    c = _mm_or_si128(c, _mm_cvtps_epi32(cmp_de_fz));

    // update matrices
    _mm_store_ps(&H[j_depth], h_simd);
    _mm_store_ps(&F[j_depth], f_simd);
    _mm_store_si128((__m128i *) &C[j_depth], c);
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
      e_simd = _mm_max_ps(_mm_sub_ps(e_simd, gap_extend_simd), 
			  _mm_sub_ps(h_simd, gap_open_simd));
 
      // diagonal value: match or mismatch
      diagonal_simd = _mm_add_ps(temp_simd,
				 _mm_set_ps((q_len[3] > j && r_len[3] > i) ? profile[(unsigned char)q[3][j]][(unsigned char)r[3][i]] : -1000.0f, 
					    (q_len[2] > j && r_len[2] > i) ? profile[(unsigned char)q[2][j]][(unsigned char)r[2][i]] : -1000.0f, 
					    (q_len[1] > j && r_len[1] > i) ? profile[(unsigned char)q[1][j]][(unsigned char)r[1][i]] : -1000.0f, 
					    (q_len[0] > j && r_len[0] > i) ? profile[(unsigned char)q[0][j]][(unsigned char)r[0][i]] : -1000.0f)
				 );

      cmp_de = _mm_min_ps(_mm_cmpge_ps(diagonal_simd, e_simd), one_simd);
      max_de = _mm_max_ps(diagonal_simd, e_simd);

      // up value: gap in query
      temp_simd = _mm_load_ps(&H[offset - q_len_depth]);

      f_simd = _mm_load_ps(&F[j_depth]);
      f_simd = _mm_max_ps(_mm_sub_ps(f_simd, gap_extend_simd), 
			  _mm_sub_ps(temp_simd, gap_open_simd));

      cmp_fz = _mm_min_ps(_mm_cmpge_ps(f_simd, zero_simd), one_simd);
      max_fz = _mm_max_ps(f_simd, zero_simd);

      // get max. value
      cmp_de_fz = _mm_min_ps(_mm_cmpge_ps(max_de, max_fz), one_simd);
      h_simd = _mm_max_ps(max_de, max_fz);

      score_simd = _mm_max_ps(score_simd, h_simd);

      // compass (save left, diagonal, up or zero?)
      c = _mm_slli_epi32(_mm_or_si128(zeroi, _mm_cvtps_epi32(cmp_de)), 1);
      c = _mm_slli_epi32(_mm_or_si128(c, _mm_cvtps_epi32(cmp_fz)), 1);
      c = _mm_or_si128(c, _mm_cvtps_epi32(cmp_de_fz));

      // update matrices
      _mm_store_ps(&H[offset], h_simd);
      _mm_store_ps(&F[j_depth], f_simd); 
      _mm_store_si128((__m128i *) &C[offset], c);

    }
  }
  _mm_store_ps(max_score, score_simd);
}

//-------------------------------------------------------------------------

void sse_matrix2(int num_seqs,
                char **q, int *q_len, int max_q_len,
                char **r, int *r_len, int max_r_len,
                float match, float mismatch, float gap_open, float gap_extend,
                float *H, float *F, int *C, float *max_score) {

  const int depth = 4;

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

  __m128i q1, r1;

  __m128 match_simd = _mm_set1_ps(match);
  __m128 mismatch_simd = _mm_set1_ps(mismatch);

  __m128 h_simd, e_simd, f_simd, diagonal_simd;
  __m128 temp_simd, subst_simd;

  __m128i zeroi   = _mm_setzero_si128();

  __m128 score_simd = _mm_setzero_ps();
  __m128 zero_simd = _mm_setzero_ps();
  __m128 one_simd = _mm_set1_ps(1);
  __m128 gap_open_simd = _mm_set1_ps(gap_open);
  __m128 gap_extend_simd = _mm_set1_ps(gap_extend);

  __m128 max_de, max_fz;
  __m128 cmp_de, cmp_fz, cmp_de_fz;
  __m128i c;

  int offset, idx, j_depth;
  int q_len_depth = depth * max_q_len;

  h_simd = zero_simd;
  e_simd = zero_simd;

  // initial loop
  r1 = _mm_load_si128((__m128i *) rr);
  for (int j = 0; j < max_q_len; j++) {

    j_depth = depth * j;

    q1 = _mm_load_si128((__m128i *) (qq + j_depth));

    // left value: gap in reference
    e_simd = _mm_max_ps(_mm_sub_ps(e_simd, gap_extend_simd),
                        _mm_sub_ps(h_simd, gap_open_simd));


    // diagonal value: match or mismatch
    subst_simd = _mm_blendv_ps(mismatch_simd, match_simd, _mm_castsi128_ps(_mm_cmpeq_epi32(q1, r1)));
    diagonal_simd = _mm_add_ps(zero_simd, subst_simd);

    cmp_de = _mm_min_ps(_mm_cmpge_ps(diagonal_simd, e_simd), one_simd);
    max_de = _mm_max_ps(diagonal_simd, e_simd);

    // up value: gap in query
    f_simd = _mm_max_ps(_mm_sub_ps(zero_simd, gap_extend_simd),
                        _mm_sub_ps(zero_simd, gap_open_simd));

    cmp_fz = _mm_min_ps(_mm_cmpge_ps(f_simd, zero_simd), one_simd);
    max_fz = _mm_max_ps(f_simd, zero_simd);

    // get max. value and save it
    cmp_de_fz = _mm_min_ps(_mm_cmpge_ps(max_de, max_fz), one_simd);
    h_simd = _mm_max_ps(max_de, max_fz);

    score_simd = _mm_max_ps(score_simd, h_simd);

    // compass (save left, diagonal, up or zero?)
    c = _mm_slli_epi32(_mm_or_si128(zeroi, _mm_cvtps_epi32(cmp_de)), 1);
    c = _mm_slli_epi32(_mm_or_si128(c, _mm_cvtps_epi32(cmp_fz)), 1);
    c = _mm_or_si128(c, _mm_cvtps_epi32(cmp_de_fz));

    // update matrices
    _mm_store_ps(&H[j_depth], h_simd);
    _mm_store_ps(&F[j_depth], f_simd);
    _mm_store_si128((__m128i *) &C[j_depth], c);
  }

  // main loop
  for (int i = 1; i < max_r_len; i++) {

    h_simd = zero_simd;
    e_simd = zero_simd;
    temp_simd = zero_simd;

    idx = i * q_len_depth;

    r1 = _mm_load_si128((__m128i *) (rr + depth * i));
    for (int j = 0; j < max_q_len; j++) {
      j_depth = depth * j;
      offset = idx + j_depth;

      // left value: gap in reference
      e_simd = _mm_max_ps(_mm_sub_ps(e_simd, gap_extend_simd),
                          _mm_sub_ps(h_simd, gap_open_simd));

      // diagonal value: match or mismatch
      q1 = _mm_load_si128((__m128i *) (qq + j_depth));
      subst_simd = _mm_blendv_ps(mismatch_simd, match_simd, _mm_castsi128_ps(_mm_cmpeq_epi32(q1, r1)));
      diagonal_simd = _mm_add_ps(temp_simd, subst_simd);

      cmp_de = _mm_min_ps(_mm_cmpge_ps(diagonal_simd, e_simd), one_simd);
      max_de = _mm_max_ps(diagonal_simd, e_simd);

      // up value: gap in query
      temp_simd = _mm_load_ps(&H[offset - q_len_depth]);

      f_simd = _mm_load_ps(&F[j_depth]);
      f_simd = _mm_max_ps(_mm_sub_ps(f_simd, gap_extend_simd),
                          _mm_sub_ps(temp_simd, gap_open_simd));

      cmp_fz = _mm_min_ps(_mm_cmpge_ps(f_simd, zero_simd), one_simd);
      max_fz = _mm_max_ps(f_simd, zero_simd);

      // get max. value
      cmp_de_fz = _mm_min_ps(_mm_cmpge_ps(max_de, max_fz), one_simd);
      h_simd = _mm_max_ps(max_de, max_fz);

      score_simd = _mm_max_ps(score_simd, h_simd);

      // compass (save left, diagonal, up or zero?)
      c = _mm_slli_epi32(_mm_or_si128(zeroi, _mm_cvtps_epi32(cmp_de)), 1);
      c = _mm_slli_epi32(_mm_or_si128(c, _mm_cvtps_epi32(cmp_fz)), 1);
      c = _mm_or_si128(c, _mm_cvtps_epi32(cmp_de_fz));

      // update matrices
      _mm_store_ps(&H[offset], h_simd);
      _mm_store_ps(&F[j_depth], f_simd);
      _mm_store_si128((__m128i *) &C[offset], c);

    }
  }
  _mm_store_ps(max_score, score_simd);

  // free memory
  _mm_free(qq);
  _mm_free(rr);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
