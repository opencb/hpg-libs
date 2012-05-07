#include "smith_waterman_simd.h"

//------------------------------------------------------------------------------------

sw_simd_context_t *sw_simd_context_new(float gap_open, float gap_extend, subst_matrix_t *subst_matrix_p) {
    sw_simd_context_t* context_p = (sw_simd_context_t*) _mm_malloc(sizeof(sw_simd_context_t), SIMD_ALIGN);

    context_p->gap_open = gap_open;
    context_p->gap_extend = gap_extend;

#ifdef __AVX__
    context_p->zero_simd = _mm256_setzero_ps();
    context_p->gap_open_simd = _mm256_set1_ps(gap_open);
    context_p->gap_extend_simd = _mm256_set1_ps(gap_extend);
#else
    context_p->zero_simd = _mm_setzero_ps();
    context_p->gap_open_simd = _mm_set1_ps(gap_open);
    context_p->gap_extend_simd = _mm_set1_ps(gap_extend);
#endif // __AVX__     

    context_p->subst_matrix_p = subst_matrix_p;
          
    context_p->x_size = 0;
    context_p->y_size = 0;
    context_p->max_size = 0;
    
    context_p->E = NULL;
    context_p->F = NULL;
    context_p->H = NULL;
    
    context_p->compass_p = NULL;
    
    context_p->a_map = NULL;
    context_p->b_map = NULL;
    
    sw_simd_context_update(200, 800, context_p);
    
    return context_p;
}

//------------------------------------------------------------------------------------

void sw_simd_context_update(int x_size, int y_size, sw_simd_context_t* context_p) {
    int max_size = x_size * y_size;
    if (context_p->x_size < x_size || context_p->y_size < y_size) {
        if (context_p->x_size < x_size) {
	    context_p->x_size = x_size;
	}
	if (context_p->y_size < y_size) {
	    context_p->y_size = y_size;
	}
	context_p->max_size = x_size * y_size;
      
	if (context_p->E != NULL) _mm_free(context_p->E);
	context_p->E = (float*) _mm_malloc(context_p->x_size * SIMD_DEPTH * sizeof(float), SIMD_ALIGN);
      
	if (context_p->F != NULL) _mm_free(context_p->F);
	context_p->F = (float*) _mm_malloc(SIMD_DEPTH * sizeof(float), SIMD_ALIGN);
      
	if (context_p->H != NULL) _mm_free(context_p->H);
	context_p->H = (float*) _mm_malloc(context_p->max_size * SIMD_DEPTH * sizeof(float), SIMD_ALIGN);
	memset(context_p->H, 0, context_p->max_size * SIMD_DEPTH * sizeof(float));
	
	if (context_p->compass_p != NULL) free(context_p->compass_p);
	context_p->compass_p = (int*) calloc(context_p->max_size * SIMD_DEPTH, sizeof(int));
	
	if (context_p->a_map != NULL) free(context_p->a_map);
	context_p->a_map = (char*) calloc(context_p->max_size + 100, sizeof(char));
      
	if (context_p->b_map != NULL) free(context_p->b_map);
	context_p->b_map = (char*) calloc(context_p->max_size + 100, sizeof(char));
    }
    
    // we must set to zero E and F vectors, however
    // for the H matrix, we only need to set to zero the first row and column
    // and it's done when allocating memory (the first row and column are not
    // overwritten)
    memset(context_p->E, 0, context_p->x_size * SIMD_DEPTH * sizeof(float));
    memset(context_p->F, 0, SIMD_DEPTH * sizeof(float));
}

//------------------------------------------------------------------------------------

void sw_simd_context_free(sw_simd_context_t* context_p) {
    if (context_p == NULL) return;

    if (context_p->E != NULL) _mm_free(context_p->E);
    if (context_p->F != NULL) _mm_free(context_p->F);
    if (context_p->H != NULL) _mm_free(context_p->H);
    
    if (context_p->compass_p != NULL) free(context_p->compass_p);
    
    if (context_p->a_map != NULL) free(context_p->a_map);
    if (context_p->b_map != NULL) free(context_p->b_map);
    
    _mm_free(context_p);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void smith_waterman_simd(char** query_p, char** ref_p, unsigned int num_queries,
			 char **query_map_p, char **ref_map_p,
			 unsigned int *query_start_p, unsigned int *ref_start_p,
			 float *score_p, sw_simd_context_t* context_p) {

     unsigned int len, query_len[SIMD_DEPTH], ref_len[SIMD_DEPTH];
     int x_size = 0, y_size = 0;
     
     for (int i = 0; i < num_queries; i++) {
	  score_p[i] = 0;

	  len = strlen(query_p[i]);
	  if (len > x_size) x_size = len;
	  query_len[i] = len;

	  len = strlen(ref_p[i]);
	  if (len > y_size) y_size = len;
	  ref_len[i] = len;
     }
     x_size++;
     y_size++;


     sw_simd_context_update(x_size, y_size, context_p);
     
     int simd_size = (SIMD_DEPTH * sizeof(float));
     int x_size_depth = (SIMD_DEPTH * x_size);
     int x_size_1_depth = SIMD_DEPTH * (x_size + 1);
     
     int *seq_x_end_p = &context_p->seq_x_end;
     int *ref_y_end_p = &context_p->ref_y_end;
     
     float **h_end_p = (float**) context_p->h_end;
     char **seq_end_p = (char**) context_p->seq_end;
     char **ref_end_p = (char**) context_p->ref_end;
     char **compass_end_p = (char**) context_p->compass_end;
     
     float *e, *E = context_p->E;
     float *f, *F = context_p->F;
     float *h, *H = context_p->H;
     
     float final[SIMD_DEPTH];
     //float max_penalty = context_p->substitution[0] * 100.0;
     //float max_penalty = context_p->substitution[0] * 100.0;
     char seq, ref;

#ifdef __AVX__
     __m256 e_simd, f_simd, h_simd, e_value_simd, f_value_simd, h_value_simd, h_final_simd;
     __m256 he_simd, hf_simd, s_value_simd, substitution_simd;

     register __m256 zero_simd = context_p->zero_simd;
     register __m256 gap_open_simd = context_p->gap_open_simd;
     register __m256 gap_extend_simd = context_p->gap_extend_simd;
#else
     __m128 e_simd, f_simd, h_simd, e_value_simd, f_value_simd, h_value_simd, h_final_simd;
     __m128 he_simd, hf_simd, s_value_simd, substitution_simd;

     register __m128 zero_simd = context_p->zero_simd;
     register __m128 gap_open_simd = context_p->gap_open_simd;
     register __m128 gap_extend_simd = context_p->gap_extend_simd;
#endif // __AVX__
     
     //printf("compunting score matrix H\n");
     
     // computing score matrix (H)
     //
     int *compass_p = (int*) context_p->compass_p + ((x_size + 1) * SIMD_DEPTH);
     h = H + ((x_size + 1) * SIMD_DEPTH);
     for (int y = 1; y < y_size; y++) {
	  e = E;
	  f = F;
	  for (int x = 1; x < x_size; x++) {
	    //	      printf("(x, y) = (%i, %i) of (%i, %i)...\n", x, y, x_size-1, y_size-1);

	       // load SIMD registers
/*
	    gettimeofday(&t1, NULL);
*/
#ifdef __AVX__
	       e_simd  = _mm256_load_ps(e);
	       he_simd = _mm256_load_ps(h - x_size_depth);	       
	       f_simd  = _mm256_load_ps(f);
	       hf_simd = _mm256_load_ps(h - SIMD_DEPTH);
	       h_simd  = _mm256_load_ps(h - x_size_1_depth);
#else
	       e_simd  = _mm_load_ps(e);
	       he_simd = _mm_load_ps(h - x_size_depth);	       
	       f_simd  = _mm_load_ps(f);
	       hf_simd = _mm_load_ps(h - SIMD_DEPTH);
	       h_simd  = _mm_load_ps(h - x_size_1_depth);
#endif // __AVX__
/*
	       gettimeofday(&t2, NULL);
	       t_diff = ((double)t2.tv_sec*1000000 + (double)t2.tv_usec) - ((double)t1.tv_sec*1000000 + (double)t1.tv_usec);
	       t_total += t_diff;
*/
	       for (int i = 0; i < num_queries; i++) {
		    if (x>query_len[i] || y>ref_len[i]) {
		      ((float*) &substitution_simd)[i] = -1000.0f; //max_penalty;
		    } else {
			 seq = query_p[i][x-1];
			 ref = ref_p[i][y-1];
			 ((float*) &substitution_simd)[i] = (*context_p->subst_matrix_p)[seq][ref];
			 //((float*) &substitution_simd)[i] = context_p->substitution[seq == ref];
		    }
	       }
/*
	       gettimeofday(&t1, NULL);
*/
#ifdef __AVX__
	       //printf("e calculation...\n");
	       // E calculation
	       e_value_simd = _mm256_sub_ps(e_simd, gap_extend_simd);
	       h_value_simd = _mm256_sub_ps(he_simd, gap_open_simd);
	       e_value_simd = _mm256_max_ps(h_value_simd, e_value_simd);
	       
	       //printf("e calculation done!!\nf calculation...\n");
	       // F calculation
	       f_value_simd = _mm256_sub_ps(f_simd, gap_extend_simd);
	       h_value_simd = _mm256_sub_ps(hf_simd, gap_open_simd);
	       f_value_simd = _mm256_max_ps(h_value_simd, f_value_simd);
	       
	       //printf("f calculation done!!\nh calculation...\n");
	       // H calculation
	       h_value_simd = _mm256_add_ps(h_simd, substitution_simd);
	       h_final_simd = _mm256_max_ps(h_value_simd, f_value_simd);
	       h_final_simd = _mm256_max_ps(h_final_simd, e_value_simd);
	       h_final_simd = _mm256_max_ps(h_final_simd, zero_simd);
#else
	       // E calculation
	       e_value_simd = _mm_sub_ps(e_simd, gap_extend_simd);
	       h_value_simd = _mm_sub_ps(he_simd, gap_open_simd);
	       e_value_simd = _mm_max_ps(h_value_simd, e_value_simd);
	       
	       // F calculation
	       f_value_simd = _mm_sub_ps(f_simd, gap_extend_simd);
	       h_value_simd = _mm_sub_ps(hf_simd, gap_open_simd);
	       f_value_simd = _mm_max_ps(h_value_simd, f_value_simd);
	       
	       // H calculation
	       h_value_simd = _mm_add_ps(h_simd, substitution_simd);
	       h_final_simd = _mm_max_ps(h_value_simd, f_value_simd);
	       h_final_simd = _mm_max_ps(h_final_simd, e_value_simd);
	       h_final_simd = _mm_max_ps(h_final_simd, zero_simd);
#endif // __AVX__
/*
	       gettimeofday(&t2, NULL);
	       t_diff = ((double)t2.tv_sec*1000000 + (double)t2.tv_usec) - ((double)t1.tv_sec*1000000 + (double)t1.tv_usec);
	       t_total += t_diff;
*/

#ifdef __AVX__
	       _mm256_store_ps(h, h_final_simd);
	       _mm256_store_ps(e, e_value_simd);
	       _mm256_store_ps(f, f_value_simd);
#else
	       _mm_store_ps(h, h_final_simd);
	       _mm_store_ps(e, e_value_simd);
	       _mm_store_ps(f, f_value_simd);
#endif // __AVX__

	       //printf("h calculation done!!\ns calculation...\n");
	       // S (maximum score)
	       for (int i = 0; i < num_queries; i++) {
		    if (h[i] == ((float *) &h_value_simd)[i]) {
		      compass_p[i] = 0;
		    } else if (h[i] == f[i]) {
		      compass_p[i] = 1;
		    } else if (h[i] == e[i]) {
		      compass_p[i] = 2;
		    } else {
		      compass_p[i] = 3;
                    }

		    if (h[i] > score_p[i]) {
			 score_p[i] = h[i];
			 seq_end_p[i] = &query_p[i][x-1];
			 ref_end_p[i] = &ref_p[i][y-1];
			 h_end_p[i] = &h[i];
			 compass_end_p[i] = &compass_p[i];
			 
			 seq_x_end_p[i] = x-1;
			 ref_y_end_p[i] = y-1;
		    }
	       }
	       //printf("s calculation done!!\nupdate vectors...\n");
	       //printf("update vectors...\n");

	       // update score matrix from SIMD registers
/*
	       gettimeofday(&t1, NULL);
*/
/*
	       gettimeofday(&t2, NULL);
	       t_diff = ((double)t2.tv_sec*1000000 + (double)t2.tv_usec) - ((double)t1.tv_sec*1000000 + (double)t1.tv_usec);
	       t_total += t_diff;
*/

	       //printf("update vectors done!!\n");
	       // update pointers
	       e += SIMD_DEPTH;
	       h += SIMD_DEPTH;
	       compass_p += SIMD_DEPTH;
	  } // end for 1..x_size
	  h += SIMD_DEPTH;
	  compass_p += SIMD_DEPTH;
	  memset(F, 0, simd_size);
     } // end for 1..y_size

     /*
     char filename[1000];
     for (int i = 0; i < SIMD_DEPTH; i++) {
       sprintf(filename, "/tmp/score-matrix-depth-%i.txt", i);
       saveFloatMatrixByIndex(H, x_size, y_size, i, SIMD_DEPTH, input_p->seq_p[i], input_p->ref_p[i], filename);
       sprintf(filename, "/tmp/compass-matrix-depth-%i.txt", i);
       saveIntegerMatrixByIndex(context_p->compass_p, x_size, y_size, i, SIMD_DEPTH, input_p->seq_p[i], input_p->ref_p[i], filename);
       printf("*** %i -> max (x, y) = (%i, %i)\n", i, seq_x_end_p[i], ref_y_end_p[i]);
     }
     */
  
     char *a_mapp, *b_mapp;

     char* a_map = context_p->a_map; // = (char*) calloc(max_size, sizeof(char));
     char* b_map = context_p->b_map;; // = (char*) calloc(max_size, sizeof(char));
     int max_size = context_p->max_size;
     
     float sc, *sc_p, *aux_sc_p; // score variables
     int compass;
     
     int x_end, num_gaps; 
     float bimble, errbounds = 0.01f;
     float gap_open = context_p->gap_open;
     float gap_extend = context_p->gap_extend;

     // walking the score matrix (H)
     for (int i = 0; i < num_queries; i++) {
	  a_mapp = a_map + max_size;
	  b_mapp = b_map + max_size;
	  
	  *a_mapp = 0;
	  //a_mapp--; *a_mapp = *seq_end_p[i];
	  
	  *b_mapp = 0;
	  //b_mapp--; *b_mapp = *ref_end_p[i];
	  
	  //x_end = seq_x_end_p[i];
	  query_start_p[i] = seq_x_end_p[i] + 1;
	  ref_start_p[i] = ref_y_end_p[i] + 1;
	  
	  sc_p = h_end_p[i];
	  compass_p = compass_end_p[i];

	  int count = 0;
	  for (;;) {
	    compass = *compass_p;
	    
	    //printf("(score, compass) = (%.2f, %i, %x)\n", score, compass, compass_end_p[i]);

	    //if (score == 0 || compass == 3) break;
	    if (compass == 3) {
	      printf("compass = 3, I don't know what to do..exit !!\n");
	      exit(-1);
	    }

	    if (compass == 0) {
	      //printf("cipf: diagonal\n");
	      sc_p -= x_size_1_depth;
	      compass_p -= x_size_1_depth;
	      a_mapp--; *a_mapp = *(seq_end_p[i]--);
	      b_mapp--; *b_mapp = *(ref_end_p[i]--);
	      query_start_p[i]--;
	      ref_start_p[i]--;
	      
	      if (*sc_p <= 0) break;
	      //continue;
	      //if (--x_end <= 0) break;

	    } else if (compass == 2) {
	      //printf("cipf: down\n");
	      num_gaps = 0;
	      sc = *sc_p;
	      aux_sc_p = sc_p - x_size_depth;

	      while (1) {
		bimble = *aux_sc_p - gap_open - (num_gaps * gap_extend);
		//printf("down: fabs(score (%0.2f) -  bimble (%0.2f)) = %0.2f, errbounds = %0.2f\n", score, bimble, fabs(score - bimble), errbounds);

		if (aux_sc_p < &H[i] || fabs(sc - bimble) < errbounds) break;
		aux_sc_p -= x_size_depth;
		++num_gaps;
	      }
	      if (bimble <= 0.0) break;

	      //if (num_gaps > 0) printf("down: num gaps = %i\n", num_gaps);

	      for (int j = 0; j <= num_gaps ; j++) {
		sc_p -= x_size_depth;
		compass_p -= x_size_depth;
		a_mapp--; *a_mapp = '-';
		b_mapp--; *b_mapp = *(--ref_end_p[i]);
		ref_start_p[i]--;
	      }
	    } else {
	      //printf("cipf: left\n");
	      num_gaps = 0;
	      sc = *sc_p;
	      aux_sc_p = sc_p - SIMD_DEPTH;
	      //printf("left: aux_score_p = %li, &H[%i] = %li\n", aux_score_p, i, &H[i]);

	      while (1) {
		bimble = *aux_sc_p - gap_open - (num_gaps * gap_extend);
		//printf("left: aux_score_p = %f, fabs(score (%f) -  bimble (%f)) = %f, errbounds = %f, condition = %i\n", *aux_score_p, score, bimble, fabs(score - bimble), errbounds, fabs(score - bimble) < errbounds);

		//printf("aux_score_p = %li, &H[%i] = %li\n", aux_score_p, i, &H[i]);
		if (aux_sc_p < &H[i] || fabs(sc - bimble) < errbounds) break;

		aux_sc_p -= SIMD_DEPTH;
		++num_gaps;
	      }
	      if (bimble <= 0.0) break;

	      //if (num_gaps > 0) 
		//printf("----------left: num gaps = %i\n", num_gaps);
	      //exit(-1);

	      for (int j = 0; j <= num_gaps ; j++) {
		sc_p -= SIMD_DEPTH;
		compass_p -= SIMD_DEPTH;
		a_mapp--; *a_mapp = *(seq_end_p[i]--);
		b_mapp--; *b_mapp = '-';
		query_start_p[i]--;
	      }
	      //if (count++ > 1) break;
	      //break;
	      //continue;
	      //if (--x_end <= 0) break;
	    }
	  } // end of for ;;
	  
	  len = a_map + max_size - a_mapp;
	  query_map_p[i] = (char*) calloc(len, sizeof(char));
	  memcpy(query_map_p[i], a_mapp, len);
	  
	  len = b_map + max_size - b_mapp;
	  ref_map_p[i] = (char*) calloc(len, sizeof(char));
	  memcpy(ref_map_p[i], b_mapp, len);
     } // end of for 0..SIMD_DEPTH
     //exit(-1);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void saveFloatMatrixByIndex(float* H, unsigned int x_size, unsigned int y_size, unsigned int index, unsigned int depth, char* A, char* B, char* filename) {
  int x, y;
  float *h;

  printf("---> saving matrix: %s\n", filename);
  FILE* fd = fopen(filename, "w");

  fprintf(fd, "\t-\t");
  for(x=0 ; x<x_size-1 ; x++) {
    fprintf(fd, "%c\t", A[x]);
  }
  fprintf(fd, "\n");

  h = H;
  for(y=0 ; y<y_size ; y++) {
    for(x=0 ; x<x_size ; x++) {
      if (x == 0) {
        if (y == 0) {
          fprintf(fd, "-\t");
        } else {
          fprintf(fd, "%c\t", B[y-1]);
        }
      }
      fprintf(fd, "%02.2f\t", h[index]);
      h += depth;
    }
    fprintf(fd, "\n");
  }

  fclose(fd);
}

void saveFloatMatrixByIndex0(float* H, unsigned int x_size, unsigned int y_size, unsigned int index, unsigned int depth, char* A, char* B, char* filename) {
  int x, y;
  float *h;

  printf("---> saving matrix: %s\n", filename);
  FILE* fd = fopen(filename, "w");

  fprintf(fd, "\t");
  for(y=0 ; y<y_size ; y++) {
    fprintf(fd, "%c\t", B[y]);
  }
  fprintf(fd, "\n");

  h = H;
  for(x=0 ; x<x_size ; x++) {
    for(y=0 ; y<y_size ; y++) {
        if (y == 0) {
	  if (x == 0) {
	    fprintf(fd, "\t");
        } else {
          fprintf(fd, "%c\t", A[x]);
        }
      }
      fprintf(fd, "%02.2f\t", h[index]);
      h += depth;
    }
    fprintf(fd, "\n");
  }

  fclose(fd);
}

//-------------------------------------------------------------

void saveIntegerMatrixByIndex(int* H, unsigned int x_size, unsigned int y_size, unsigned int index, unsigned int depth, char* A, char* B, char* filename) {
  int x, y;
  int *h;

  printf("---> saving matrix: %s\n", filename);
  FILE* fd = fopen(filename, "w");

  fprintf(fd, "\t-\t");
  for(x=0 ; x<x_size-1 ; x++) {
    fprintf(fd, "%c\t", A[x]);
  }
  fprintf(fd, "\n");

  h = H;
  for(y=0 ; y<y_size ; y++) {
    for(x=0 ; x<x_size ; x++) {
      if (x == 0) {
        if (y == 0) {
          fprintf(fd, "-\t");
        } else {
          fprintf(fd, "%c\t", B[y-1]);
        }
      }
      fprintf(fd, "%i\t", h[index]);
      h += depth;
    }
    fprintf(fd, "\n");
  }

  fclose(fd);
}

void saveIntegerMatrixByIndex0(int* H, unsigned int x_size, unsigned int y_size, unsigned int index, unsigned int depth, char* A, char* B, char* filename) {
  int x, y;
  int *h;

  printf("---> saving matrix: %s\n", filename);
  FILE* fd = fopen(filename, "w");

  fprintf(fd, "\t");
  for(y=0 ; y<y_size ; y++) {
    fprintf(fd, "%c\t", B[y]);
  }
  fprintf(fd, "\n");

  h = H;
  for(x=0 ; x<x_size ; x++) {
    for(y=0 ; y<y_size ; y++) {
      if (y == 0) {
	if (x == 0) {
          fprintf(fd, "\t");
        } else {
          fprintf(fd, "%c\t", A[x]);
        }
      }
      fprintf(fd, "%i\t", h[index]);
      h += depth;
    }
    fprintf(fd, "\n");
  }

  fclose(fd);
}

//-------------------------------------------------------------

void saveCharMatrixByIndex(char* H, unsigned int x_size, unsigned int y_size, unsigned int index, unsigned int depth, char* A, char* B, char* filename) {
  int x, y;
  int *h;

  printf("---> saving matrix: %s\n", filename);
  FILE* fd = fopen(filename, "w");

  fprintf(fd, "\t-\t");
  for(x=0 ; x<x_size-1 ; x++) {
    fprintf(fd, "%c\t", A[x]);
  }
  fprintf(fd, "\n");

  h = H;
  for(y=0 ; y<y_size ; y++) {
    for(x=0 ; x<x_size ; x++) {
      if (x == 0) {
        if (y == 0) {
          fprintf(fd, "-\t");
        } else {
          fprintf(fd, "%c\t", B[y-1]);
        }
      }
      fprintf(fd, "%i\t", h[index]);
      h += depth;
    }
    fprintf(fd, "\n");
  }

  fclose(fd);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
