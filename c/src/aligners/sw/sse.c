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

  __m128i zeroi = _mm_set_epi32(0, 0, 0, 0);

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
  /*
      for (int i = 0; i < 4; i++) {
        printf("query %i:%s\nref.  %i:%s\n\n", i, q[i], i, r[i]);
      }
  */
  h_simd = zero_simd;
  e_simd = zero_simd;
  
  for (int j = 0; j < max_q_len; j++) {

    j_depth = depth * j;
    
    // left value: gap in reference
    e_simd = _mm_max_ps(_mm_sub_ps(e_simd, gap_extend_simd), 
			_mm_sub_ps(h_simd, gap_open_simd));

    //    printf("from left: %0.2f\n", ((float *)&e_simd)[0]);
    
    // diagonal value: match or mismatch
    subst_simd = _mm_set_ps((q_len[3] > j) ? profile[(unsigned char)q[3][j]][(unsigned char)r[3][0]] : -1000.0f,
                            (q_len[2] > j) ? profile[(unsigned char)q[2][j]][(unsigned char)r[2][0]] : -1000.0f,
                            (q_len[1] > j) ? profile[(unsigned char)q[1][j]][(unsigned char)r[1][0]] : -1000.0f,
                            (q_len[0] > j) ? profile[(unsigned char)q[0][j]][(unsigned char)r[0][0]] : -1000.0f);
    //printf("Init end\n");
    /*
    subst_simd = _mm_set_ps(profile[q[3][j]][r[3][0]], 
			    profile[q[2][j]][r[2][0]], 
			    profile[q[1][j]][r[1][0]], 
			    profile[q[0][j]][r[0][0]]);
    */

    diagonal_simd = _mm_add_ps(zero_simd, subst_simd);
    //    printf("from diagonal: temp = %0.2f %0.2f (%c, %c) -> %0.2f\n", ((float *)&temp_simd)[0], profile[q[0][j]][r[0][0]], q[0][j], r[0][0], ((float *)&diagonal_simd)[0]);

    cmp_de = _mm_min_ps(_mm_cmpge_ps(diagonal_simd, e_simd), one_simd);
    max_de = _mm_max_ps(diagonal_simd, e_simd);

    // up value: gap in query
    f_simd = _mm_max_ps(_mm_sub_ps(zero_simd, gap_extend_simd), 
			_mm_sub_ps(zero_simd, gap_open_simd));

    cmp_fz = _mm_min_ps(_mm_cmpge_ps(f_simd, zero_simd), one_simd);
    max_fz = _mm_max_ps(f_simd, zero_simd);
    
    //    printf("from up: %0.2f\n", ((float *)&f_simd)[0]);
    
    // get max. value and save it
    cmp_de_fz = _mm_min_ps(_mm_cmpge_ps(max_de, max_fz), one_simd);
    h_simd = _mm_max_ps(max_de, max_fz);

    score_simd = _mm_max_ps(score_simd, h_simd);
    //    printf("\t\t\t\t\tmax. score: %0.2f\n", ((float *)&h_simd)[0]);
    
    // compass (save left, diagonal, up or zero?)
    c = _mm_slli_epi32(_mm_or_si128(zeroi, _mm_cvtps_epi32(cmp_de)), 1);
    c = _mm_slli_epi32(_mm_or_si128(c, _mm_cvtps_epi32(cmp_fz)), 1);
    c = _mm_or_si128(c, _mm_cvtps_epi32(cmp_de_fz));

    //    printf("\t\t\t\t\tcompass: %i\n", ((int *)&c)[0]);

    // update matrices
    _mm_store_ps(&H[j_depth], h_simd);
    _mm_store_ps(&F[j_depth], f_simd);
    _mm_store_si128((__m128i *)&C[j_depth], c);

 
    //_mm_store_ps(&D[j_depth], diagonal_simd);

    /*
    offset = j_depth;
    printf("(row, col) = (%i, %i):\t \t%c-%c=%0.2f %c-%c=%0.2f %c-%c=%0.2f %c-%c=%0.2f\n", 0, j, q[0][j], r[0][0], profile[q[0][j]][r[0][0]], q[1][j], r[1][0], profile[q[1][j]][r[1][0]], q[2][j], r[2][0], profile[q[2][j]][r[2][0]], q[3][j], r[3][0], profile[q[3][j]][r[3][0]]);
    printf("(row, col) = (%i, %i):\tH\t%0.2f %0.2f %0.2f %0.2f\n", 0, j, H[offset], H[offset+1], H[offset+2], H[offset+3]);
    printf("(row, col) = (%i, %i):\tD\t%0.2f %0.2f %0.2f %0.2f\n", 0, j, D[offset], D[offset+1], D[offset+2], D[offset+3]);
    printf("(row, col) = (%i, %i):\td\t%0.2f %0.2f %0.2f %0.2f\n", 0, j, ((float *)&diagonal_simd)[0], ((float *)&diagonal_simd)[1], ((float *)&diagonal_simd)[2], ((float *)&diagonal_simd)[3]);

    printf("(row, col) = (%i, %i):\ts\t%0.2f %0.2f %0.2f %0.2f\n", 0, j, ((float *)&subst_simd)[0], ((float *)&subst_simd)[1], ((float *)&subst_simd)[2], ((float *)&subst_simd)[3]);
    */
  }
  //  printf("\n");

  //  exit(-1);

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
 
      //      if (i == 3 && j == 3) printf("from left: %0.2f\n", ((float *)&e_simd)[target]);

      // diagonal value: match or mismatch
      diagonal_simd = _mm_add_ps(temp_simd,
				 _mm_set_ps((q_len[3] > j && r_len[3] > i) ? profile[(unsigned char)q[3][j]][(unsigned char)r[3][i]] : -1000.0f, 
					    (q_len[2] > j && r_len[2] > i) ? profile[(unsigned char)q[2][j]][(unsigned char)r[2][i]] : -1000.0f, 
					    (q_len[1] > j && r_len[1] > i) ? profile[(unsigned char)q[1][j]][(unsigned char)r[1][i]] : -1000.0f, 
					    (q_len[0] > j && r_len[0] > i) ? profile[(unsigned char)q[0][j]][(unsigned char)r[0][i]] : -1000.0f)
				 );

      cmp_de = _mm_min_ps(_mm_cmpge_ps(diagonal_simd, e_simd), one_simd);
      max_de = _mm_max_ps(diagonal_simd, e_simd);

      //      if (i == 3 && j == 3)	printf("from diagonal: temp = %0.2f %0.2f (%c, %c) -> %0.2f\n", ((float *)&temp_simd)[target], profile[q[target][j]][r[target][i]], q[target][j], r[target][i], ((float *)&diagonal_simd)[target]);
      
      // up value: gap in query
      temp_simd = _mm_load_ps(&H[offset - q_len_depth]);

      f_simd = _mm_load_ps(&F[j_depth]);
      f_simd = _mm_max_ps(_mm_sub_ps(f_simd, gap_extend_simd), 
			  _mm_sub_ps(temp_simd, gap_open_simd));

      cmp_fz = _mm_min_ps(_mm_cmpge_ps(f_simd, zero_simd), one_simd);
      max_fz = _mm_max_ps(f_simd, zero_simd);

      //      if (i == 3 && j == 3) printf("from up: %0.2f\n", ((float *)&f_simd)[target]);

      // get max. value
      cmp_de_fz = _mm_min_ps(_mm_cmpge_ps(max_de, max_fz), one_simd);
      h_simd = _mm_max_ps(max_de, max_fz);

      score_simd = _mm_max_ps(score_simd, h_simd);

      //      if (i == 3 && j == 3) printf("\t\t\t\t\tmax. score: %0.2f\n", ((float *)&h_simd)[target]);

      // compass (save left, diagonal, up or zero?)
      c = _mm_slli_epi32(_mm_or_si128(zeroi, _mm_cvtps_epi32(cmp_de)), 1);
      c = _mm_slli_epi32(_mm_or_si128(c, _mm_cvtps_epi32(cmp_fz)), 1);
      c = _mm_or_si128(c, _mm_cvtps_epi32(cmp_de_fz));

      // update matrices
      _mm_store_ps(&H[offset], h_simd);
      _mm_store_ps(&F[j_depth], f_simd); 
      _mm_store_si128((__m128i *)&C[offset], c);

      /*
      if (j==0) {
	printf("(row, col) = (%i, %i):\tD\t%0.2f %0.2f %0.2f %0.2f\n", i, j, D[offset], D[offset+1], D[offset+2], D[offset+3]);
	printf("(row, col) = (%i, %i):\tH\t%0.2f %0.2f %0.2f %0.2f\n", i, j, H[offset], H[offset+1], H[offset+2], H[offset+3]);
      }
      */
      //      printf("(row, col) = (%i, %i):\t%0.2f %0.2f %0.2f %0.2f\n", i, j, H[offset], H[offset+1], H[offset+2], H[offset+3]);
    }
    //    printf("\n");
  }
  _mm_store_ps(max_score, score_simd);
  //  max_score[0] = 100.0f;
  //  max_score[1] = 100.0f;
  //  max_score[2] = 100.0f;
  //  max_score[3] = 100.0f;
  
  /*

  int rr_len = r_len[0];
  int qq_len = q_len[0];
  printf("r_len[0] = %i, q_len[0] = %i\n", rr_len, qq_len);

  printf("sse\n");
  for (int i = 0; i < rr_len; i++) {
    printf("\t");
    for (int j = 0; j < qq_len; j++) {
      printf("%0.2f\t", H[(i * max_q_len * 4) + (j * 4)]);
    }
    printf("\n");
  }
  */
  /*
  char filename[200];
  for (int i = 0; i < 4; i++) {
    sprintf(filename, "/tmp/sse1-%i.score", i);
    save_float_matrix(H, max_q_len, max_r_len, q[i], q_len[i], r[i], r_len[i], i, 4, filename);
  }
  */
  /*
      for (int i = 0; i < 4; i++) {
        printf("score %i:%0.2f\n\n", i, max_score[i]);
      }
  */
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
