#include "avx.h"

extern double avx_matrix_t, avx_tracking_t;
extern double avx1_matrix_t, avx1_tracking_t;


//float *H = NULL;//(float *) _mm_malloc(4 * max_q_len * max_r_len * sizeof(float), 16);
//int *C = NULL; //(int *) _mm_malloc(4 * max_q_len * max_r_len * sizeof(int), 16);

void sw_avx(char **q, char **r, float gap_open, float gap_extend, 
	    float *score, char **m, char **n, int *start1, int *start2) {
  /*  
    if (H == NULL) {
      H = (float *) _mm_malloc(4 * 100 * 800 * sizeof(float), 16);
      C = (int *) _mm_malloc(4 * 100 * 800 * sizeof(int), 16);
    }
  */
  float match = 5.0f, mismatch = -4.0f;

  int len, q_lens[8], r_lens[8], alig_lens[8];
  int max_q_len = 0;
  int max_r_len = 0;
  for (int i = 0; i < 8; i++) {
    len = strlen(q[i]);
    if (len > max_q_len) max_q_len = len;
    q_lens[i] = len;

    len = strlen(r[i]);
    if (len > max_r_len) max_r_len = len;
    r_lens[i] = len;

    //printf("length of q[%i] = %i\n", i, q_lens[i]);
    //printf("length of r[%i] = %i\n", i, r_lens[i]);
  }
  //  printf("max. length of query = %i, ref. = %i\n", max_q_len, max_r_len);

  float profile[128][128];
  profile['A']['A'] = 5.0f; profile['C']['A'] = -4.0f; profile['T']['A'] = -4.0f; profile['G']['A'] = -4.0f;
  profile['A']['C'] = -4.0f; profile['C']['C'] = 5.0f; profile['T']['C'] = -4.0f; profile['G']['C'] = -4.0f;
  profile['A']['G'] = -4.0f; profile['C']['T'] = -4.0f; profile['T']['T'] = 5.0f; profile['G']['T'] = -4.0f;
  profile['A']['T'] = -4.0f; profile['C']['G'] = -4.0f; profile['T']['G'] = -4.0f; profile['G']['G'] = 5.0f;

  int matrix_size = 8 * max_q_len * max_r_len * sizeof(float);

  float *H = (float *) _mm_malloc(matrix_size, 32);
  float *D = (float *) _mm_malloc(matrix_size, 32);
  float *E = (float *) _mm_malloc(matrix_size, 32);
  float *F = (float *) _mm_malloc(matrix_size, 32);

  double partial_t;
  partial_t = tic();
  avx_matrix(8, q, q_lens, max_q_len, r, r_lens, max_r_len, 
	     profile, gap_open, gap_extend, H, D, E, F, score);
  avx_matrix_t += toc(partial_t);

  /*
  char filename[200];
  for (int i = 0; i < 4; i++) {
    sprintf(filename, "/tmp/sse-%i.score", i);
    save_float_matrix(H, max_q_len, max_r_len,
		      q[i], q_lens[i], r[i], r_lens[i], 
		      i, 4, filename);
  }
  */

  partial_t = tic();
  simd_traceback(8, 8, q, q_lens, max_q_len, r, r_lens, max_r_len,
		 gap_open, gap_extend,
		 H, D, E, F, score,
		 m, start1,
		 n, start2, 
		 alig_lens);
  avx_tracking_t += toc(partial_t);
  
  _mm_free(H);
  _mm_free(D);
  _mm_free(E);
  _mm_free(F);
  //  _mm_free(M);
}

//-------------------------------------------------------------------------

void avx_matrix(int num_seqs, 
		char **q, int *q_len, int max_q_len,
		char **r, int *r_len, int max_r_len,
		float profile[128][128], float gap_open, float gap_extend,
		float *H, float *D, float *E, float *F, float *max_score) {
  
  const int depth = 8;

  __m256 h_simd, e_simd, f_simd, diagonal_simd;
  __m256 score_simd, temp_simd, subst_simd;

  __m256 zero_simd = _mm256_setzero_ps();
  __m256 gap_open_simd = _mm256_set1_ps(gap_open);
  __m256 gap_extend_simd = _mm256_set1_ps(gap_extend);

  int offset, index, idx, j_depth;
  int q_len_depth = depth * max_q_len;

  //    for (int i = 0; i < 4; i++) {
  //      printf("query %i:%s\nref.  %i:%s\n\n", i, q[i], i, r[i]);
  //    }

  h_simd = zero_simd;
  e_simd = zero_simd;
  
  for (int j = 0; j < max_q_len; j++) {

    j_depth = depth * j;
    
    // left value: gap in reference
    e_simd = _mm256_max_ps(_mm256_sub_ps(e_simd, gap_extend_simd), 
			  _mm256_sub_ps(h_simd, gap_open_simd));
    
    
    // diagonal value: match or mismatch
    subst_simd = _mm256_set_ps(5, 5, 5, 5, 5, 5, 5, 5);
      /*
    subst_simd = _mm256_set_ps(profile[q[7][j]][r[7][0]], 
			       profile[q[6][j]][r[6][0]], 
			       profile[q[5][j]][r[5][0]], 
			       profile[q[4][j]][r[4][0]], 
			       profile[q[3][j]][r[3][0]], 
			       profile[q[2][j]][r[2][0]], 
			       profile[q[1][j]][r[1][0]], 
			       profile[q[0][j]][r[0][0]]);
      */
    diagonal_simd = _mm256_add_ps(zero_simd, subst_simd);
    

    // up value: gap in query
    f_simd = _mm256_max_ps(_mm256_sub_ps(zero_simd, gap_extend_simd), 
			  _mm256_sub_ps(zero_simd, gap_open_simd));
    
    
    // get max. value
    h_simd = _mm256_max_ps(_mm256_max_ps(diagonal_simd, 
				       _mm256_max_ps(e_simd, f_simd)), 
			  zero_simd);
    score_simd = _mm256_max_ps(score_simd, h_simd);

    
    // update matrices
    _mm256_store_ps(&H[j_depth], h_simd);
    _mm256_store_ps(&E[j_depth], e_simd);
    _mm256_store_ps(&F[j_depth], f_simd); 
    _mm256_store_ps(&D[j_depth], diagonal_simd);

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
      
 
      //      if (i == 6 && j == 3) printf("from left: %0.2f\n", ((float *)&e_simd)[0]);

      // diagonal value: match or mismatch
      diagonal_simd = _mm256_add_ps(temp_simd,
				    _mm256_set_ps(5, 5, 5, 5, 5, 5, 5, 5) 
						  /*
				    _mm256_set_ps(profile[q[7][j]][r[7][i]], 
						  profile[q[6][j]][r[6][i]], 
						  profile[q[5][j]][r[5][i]], 
						  profile[q[4][j]][r[4][i]], 
						  profile[q[3][j]][r[3][i]], 
						  profile[q[2][j]][r[2][i]], 
						  profile[q[1][j]][r[1][i]], 
						  profile[q[0][j]][r[0][i]])
				    */
				    );

      //      if (i == 6 && j == 3) printf("from diagonal: temp = %0.2f %0.2f (%c, %c) -> %0.2f\n", ((float *)&temp_simd)[0], profile[q[0][j]][r[0][i]], q[0][j], r[0][i], ((float *)&diagonal_simd)[0]);
      
      // up value: gap in query
      temp_simd = _mm256_load_ps(&H[offset - q_len_depth]);
      
      f_simd = _mm256_load_ps(&F[offset - q_len_depth]);
      f_simd = _mm256_max_ps(_mm256_sub_ps(f_simd,gap_extend_simd), 
			     _mm256_sub_ps(temp_simd, gap_open_simd));
      

      //      if (i == 6 && j == 3) printf("from up: %0.2f\n", ((float *)&f_simd)[0]);

      // get max. value
      h_simd = _mm256_max_ps(_mm256_max_ps(diagonal_simd, 
					 _mm256_max_ps(e_simd, f_simd)), 
			    zero_simd);
      score_simd = _mm256_max_ps(score_simd, h_simd);

      //      if (i == 6 && j == 3) printf("max. score: %0.2f\n", ((float *)&h_simd)[0]);

      // update matrices
      _mm256_store_ps(&H[offset], h_simd);
      _mm256_store_ps(&E[offset], e_simd);
      _mm256_store_ps(&F[offset], f_simd); 
      _mm256_store_ps(&D[offset], diagonal_simd);
      
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
  _mm256_store_ps(max_score, score_simd);
  
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
  char filename[200] = "/tmp/sse-0.score";
  save_float_matrix(H, max_q_len, max_r_len, q[0], q_len[0], r[0], r_len[0], 0, 4, filename);
  */  
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void sw_avx1(char **q, char **r, float gap_open, float gap_extend, 
	     float *score, char **m, char **n, int *start1, int *start2) {
  /*  
    if (H == NULL) {
      H = (float *) _mm_malloc(4 * 100 * 800 * sizeof(float), 16);
      C = (int *) _mm_malloc(4 * 100 * 800 * sizeof(int), 16);
    }
  */
  float match = 5.0f, mismatch = -4.0f;
  
  int len, q_lens[8], r_lens[8], alig_lens[8];
  int max_q_len = 0;
  int max_r_len = 0;
  for (int i = 0; i < 8; i++) {
    len = strlen(q[i]);
    if (len > max_q_len) max_q_len = len;
    q_lens[i] = len;

    len = strlen(r[i]);
    if (len > max_r_len) max_r_len = len;
    r_lens[i] = len;

    //printf("length of q[%i] = %i\n", i, q_lens[i]);
    //printf("length of r[%i] = %i\n", i, r_lens[i]);
  }
  //  printf("max. length of query = %i, ref. = %i\n", max_q_len, max_r_len);

  float profile[128][128];
  profile['A']['A'] = 5.0f; profile['C']['A'] = -4.0f; profile['T']['A'] = -4.0f; profile['G']['A'] = -4.0f;
  profile['A']['C'] = -4.0f; profile['C']['C'] = 5.0f; profile['T']['C'] = -4.0f; profile['G']['C'] = -4.0f;
  profile['A']['G'] = -4.0f; profile['C']['T'] = -4.0f; profile['T']['T'] = 5.0f; profile['G']['T'] = -4.0f;
  profile['A']['T'] = -4.0f; profile['C']['G'] = -4.0f; profile['T']['G'] = -4.0f; profile['G']['G'] = 5.0f;

  int matrix_size = 8 * max_q_len * max_r_len * sizeof(float);

  float *H = (float *) _mm_malloc(matrix_size, 32);
  float *F = (float *) _mm_malloc(4 * max_q_len * sizeof(float), 32);
  int *C = (int *) _mm_malloc(matrix_size, 32);

  double partial_t;
  partial_t = tic();
  avx1_matrix(8, q, q_lens, max_q_len, r, r_lens, max_r_len, 
	      profile, gap_open, gap_extend, H, F, C, score);
  avx1_matrix_t += toc(partial_t);
  
  /*
  char filename[200];
  for (int i = 0; i < 4; i++) {
    sprintf(filename, "/tmp/sse-%i.score", i);
    save_float_matrix(H, max_q_len, max_r_len,
		      q[i], q_lens[i], r[i], r_lens[i], 
		      i, 4, filename);
  }
  */

  partial_t = tic();
  simd_traceback1(8, 8, q, q_lens, max_q_len, r, r_lens, max_r_len,
		  gap_open, gap_extend,
		  H, C, score,
		  m, start1,
		  n, start2, 
		  alig_lens);
  avx1_tracking_t += toc(partial_t);
  
  _mm_free(H);
  _mm_free(F);
  _mm_free(C);
}

//-------------------------------------------------------------------------

void avx1_matrix(int num_seqs, 
		 char **q, int *q_len, int max_q_len,
		 char **r, int *r_len, int max_r_len,
		 float profile[128][128], float gap_open, float gap_extend,
		 float *H, float *F, int *C, float *max_score) {
  
  const int depth = 8;

  __m256 h_simd, e_simd, f_simd, diagonal_simd;
  __m256 temp_simd, subst_simd;

  __m256i zeroi   = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 0);

  __m256 score_simd = _mm256_setzero_ps();
  __m256 zero_simd = _mm256_setzero_ps();
  __m256 one_simd = _mm256_set1_ps(1);
  __m256 gap_open_simd = _mm256_set1_ps(gap_open);
  __m256 gap_extend_simd = _mm256_set1_ps(gap_extend);

  __m256 max_de, max_fz;
  __m256 cmp_de, cmp_fz, cmp_de_fz;
  __m256i c;

  int offset, index, idx, j_depth;
  int q_len_depth = depth * max_q_len;

  //    for (int i = 0; i < 4; i++) {
  //      printf("query %i:%s\nref.  %i:%s\n\n", i, q[i], i, r[i]);
  //    }

  h_simd = zero_simd;
  e_simd = zero_simd;
  
  for (int j = 0; j < max_q_len; j++) {

    j_depth = depth * j;
    
    // left value: gap in reference
    e_simd = _mm256_max_ps(_mm256_sub_ps(e_simd, gap_extend_simd), 
			   _mm256_sub_ps(h_simd, gap_open_simd));

    //    printf("from left: %0.2f\n", ((float *)&e_simd)[0]);
    
    // diagonal value: match or mismatch
    subst_simd = _mm256_set_ps(profile[q[7][j]][r[7][0]], 
			       profile[q[6][j]][r[6][0]], 
			       profile[q[5][j]][r[5][0]], 
			       profile[q[4][j]][r[4][0]], 
			       profile[q[3][j]][r[3][0]], 
			       profile[q[2][j]][r[2][0]], 
			       profile[q[1][j]][r[1][0]], 
			       profile[q[0][j]][r[0][0]]);
    
    diagonal_simd = _mm256_add_ps(zero_simd, subst_simd);
    //    printf("from diagonal: temp = %0.2f %0.2f (%c, %c) -> %0.2f\n", ((float *)&temp_simd)[0], profile[q[0][j]][r[0][0]], q[0][j], r[0][0], ((float *)&diagonal_simd)[0]);

    cmp_de = _mm256_min_ps(_mm256_cmp_ps(diagonal_simd, e_simd, _CMP_GE_OQ), one_simd);
    max_de = _mm256_max_ps(diagonal_simd, e_simd);

    // up value: gap in query
    f_simd = _mm256_max_ps(_mm256_sub_ps(zero_simd, gap_extend_simd), 
			   _mm256_sub_ps(zero_simd, gap_open_simd));

    cmp_fz = _mm256_min_ps(_mm256_cmp_ps(f_simd, zero_simd, _CMP_GE_OQ), one_simd);
    max_fz = _mm256_max_ps(f_simd, zero_simd);
    
    //    printf("from up: %0.2f\n", ((float *)&f_simd)[0]);
    
    // get max. value and save it
    cmp_de_fz = _mm256_min_ps(_mm256_cmp_ps(max_de, max_fz, _CMP_GE_OQ), one_simd);
    h_simd = _mm256_max_ps(max_de, max_fz);

    score_simd = _mm256_max_ps(score_simd, h_simd);
    //    printf("\t\t\t\t\tmax. score: %0.2f\n", ((float *)&h_simd)[0]);
    
    // compass (save left, diagonal, up or zero?)
    c = _mm256_slli_epi32(_mm256_or_si256(zeroi, _mm256_cvtps_epi32(cmp_de)), 1);
    c = _mm256_slli_epi32(_mm256_or_si256(c, _mm256_cvtps_epi32(cmp_fz)), 1);
    c = _mm256_or_si256(c, _mm256_cvtps_epi32(cmp_de_fz));

    //    printf("\t\t\t\t\tcompass: %i\n", ((int *)&c)[0]);

    // update matrices
    _mm256_store_ps(&H[j_depth], h_simd);
    _mm256_store_ps(&F[j_depth], f_simd);
    _mm256_store_si256(&C[j_depth], c);

 
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
  int target = 0;
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
 
      //      if (i == 3 && j == 3) printf("from left: %0.2f\n", ((float *)&e_simd)[target]);

      // diagonal value: match or mismatch
      diagonal_simd = _mm256_add_ps(temp_simd,
				    _mm256_set_ps((q_len[7] > j && r_len[7] > i) ? profile[q[7][j]][r[7][i]] : -1000.0f, 
						  (q_len[6] > j && r_len[6] > i) ? profile[q[6][j]][r[6][i]] : -1000.0f, 
						  (q_len[5] > j && r_len[5] > i) ? profile[q[5][j]][r[5][i]] : -1000.0f, 
						  (q_len[4] > j && r_len[4] > i) ? profile[q[4][j]][r[4][i]] : -1000.0f, 
						  (q_len[3] > j && r_len[3] > i) ? profile[q[3][j]][r[3][i]] : -1000.0f, 
						  (q_len[2] > j && r_len[2] > i) ? profile[q[2][j]][r[2][i]] : -1000.0f, 
						  (q_len[1] > j && r_len[1] > i) ? profile[q[1][j]][r[1][i]] : -1000.0f, 
						  (q_len[0] > j && r_len[0] > i) ? profile[q[0][j]][r[0][i]] : -1000.0f)
				 );

      cmp_de = _mm256_min_ps(_mm256_cmp_ps(diagonal_simd, e_simd, _CMP_GE_OQ), one_simd);
      max_de = _mm256_max_ps(diagonal_simd, e_simd);

      //      if (i == 3 && j == 3)	printf("from diagonal: temp = %0.2f %0.2f (%c, %c) -> %0.2f\n", ((float *)&temp_simd)[target], profile[q[target][j]][r[target][i]], q[target][j], r[target][i], ((float *)&diagonal_simd)[target]);
      
      // up value: gap in query
      temp_simd = _mm256_load_ps(&H[offset - q_len_depth]);

      f_simd = _mm256_load_ps(&F[j_depth]);
      f_simd = _mm256_max_ps(_mm256_sub_ps(f_simd, gap_extend_simd), 
			     _mm256_sub_ps(temp_simd, gap_open_simd));

      cmp_fz = _mm256_min_ps(_mm256_cmp_ps(f_simd, zero_simd, _CMP_GE_OQ), one_simd);
      max_fz = _mm256_max_ps(f_simd, zero_simd);

      //      if (i == 3 && j == 3) printf("from up: %0.2f\n", ((float *)&f_simd)[target]);

      // get max. value
      cmp_de_fz = _mm256_min_ps(_mm256_cmp_ps(max_de, max_fz, _CMP_GE_OQ), one_simd);
      h_simd = _mm256_max_ps(max_de, max_fz);

      score_simd = _mm256_max_ps(score_simd, h_simd);

      //      if (i == 3 && j == 3) printf("\t\t\t\t\tmax. score: %0.2f\n", ((float *)&h_simd)[target]);

      // compass (save left, diagonal, up or zero?)
      c = _mm256_slli_epi32(_mm256_or_ps(zeroi, _mm256_cvtps_epi32(cmp_de)), 1);
      c = _mm256_slli_epi32(_mm256_or_ps(c, _mm256_cvtps_epi32(cmp_fz)), 1);
      c = _mm256_or_ps(c, _mm256_cvtps_epi32(cmp_de_fz));

      // update matrices
      _mm256_store_ps(&H[offset], h_simd);
      _mm256_store_ps(&F[j_depth], f_simd); 
      _mm256_store_si256(&C[offset], c);

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
  _mm256_store_ps(max_score, score_simd);
  
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
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
