#include "sw.h"

extern struct timeval t1, t2;
extern double t_diff, t_total; 

//====================================================================================
// Smith-Waterman structures and functions (SIMD version)
//====================================================================================

sw_simd_input_t* sw_simd_input_new(unsigned int depth) {
     sw_simd_input_t* input_p = calloc(1, sizeof(sw_simd_input_t));
     
     input_p->depth = depth;
     input_p->seq_p = (char**) calloc(depth, sizeof(char*));
     input_p->ref_p = (char**) calloc(depth, sizeof(char*));
     input_p->seq_len_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     input_p->ref_len_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     
     return input_p;
}

//------------------------------------------------------------------------------------

void sw_simd_input_free(sw_simd_input_t* input_p) {
     if (input_p == NULL) return;

     if (input_p->seq_p != NULL) free(input_p->seq_p);
     if (input_p->ref_p != NULL) free(input_p->ref_p);
     if (input_p->seq_len_p != NULL) free(input_p->seq_len_p);
     if (input_p->ref_len_p != NULL) free(input_p->ref_len_p);
     
     free(input_p);
}

//------------------------------------------------------------------------------------

void sw_simd_input_add(char* seq_p, unsigned int seq_len, char* ref_p,
		       unsigned int ref_len, unsigned int index, sw_simd_input_t* input_p) {
     if (index >= input_p->depth) {
	  printf("ERROR: out of range!\n");
	  return;
     }
     
     input_p->seq_p[index] = seq_p;
     input_p->ref_p[index] = ref_p;
     input_p->seq_len_p[index] = seq_len;
     input_p->ref_len_p[index] = ref_len;
}

//------------------------------------------------------------------------------------

void sw_simd_input_display(unsigned int depth, sw_simd_input_t* input_p) {
     char str[2048];
     
     printf("sw_simd_input (depth = %i of %i):\n", depth, input_p->depth);
     for(int i = 0; i < depth; i++) {	  
	  memcpy(str, input_p->seq_p[i], input_p->seq_len_p[i]);
	  str[input_p->seq_len_p[i]] = '\0';
	  printf("%i:\n", i);
	  printf("seq: %s (%i)\n", str, input_p->seq_len_p[i]);
	  
	  memcpy(str, input_p->ref_p[i], input_p->ref_len_p[i]);
	  str[input_p->ref_len_p[i]] = '\0';
	  printf("ref: %s (%i)\n", str, input_p->ref_len_p[i]);
	  
	  printf("\n");
     }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

sw_simd_output_t* sw_simd_output_new(unsigned int depth) {
     sw_simd_output_t* output_p = calloc(1, sizeof(sw_simd_output_t));
     
     output_p->depth = depth;
     output_p->mapped_seq_p = (char**) calloc(depth, sizeof(char*));
     output_p->mapped_ref_p = (char**) calloc(depth, sizeof(char*));
     output_p->mapped_len_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     //output_p->mapped_ref_len_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     
     output_p->start_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     output_p->start_seq_p = (unsigned int*) calloc(depth, sizeof(unsigned int));
     
     output_p->score_p = (float*) calloc(depth, sizeof(float));
     output_p->norm_score_p = (float*) calloc(depth, sizeof(float));
     
     return output_p;
}

//------------------------------------------------------------------------------------

void sw_simd_output_free(sw_simd_output_t* output_p) {
     if (output_p == NULL) return;
     
     if (output_p->mapped_seq_p != NULL) free(output_p->mapped_seq_p);
     if (output_p->mapped_ref_p != NULL) free(output_p->mapped_ref_p);
     if (output_p->mapped_len_p != NULL) free(output_p->mapped_len_p);
     //if (output_p->mapped_ref_len_p != NULL) free(output_p->mapped_ref_len_p);
     
     if (output_p->start_p != NULL) free(output_p->start_p);
     if (output_p->start_seq_p != NULL) free(output_p->start_seq_p);
     if (output_p->score_p != NULL) free(output_p->score_p);
     if (output_p->norm_score_p != NULL) free(output_p->norm_score_p);
     
     free(output_p);
}
//------------------------------------------------------------------------------------

void sw_simd_output_display(unsigned int depth, sw_simd_output_t* output_p) {
     char str[2048];
     
     printf("sw_simd_output (depth = %i of %i):\n", depth, output_p->depth);
     for(int i = 0; i < depth; i++) {
	  
	  memcpy(str, output_p->mapped_seq_p[i], output_p->mapped_len_p[i]);
	  str[output_p->mapped_len_p[i]] = '\0';
	  printf("%i:\n", i);
	  printf("seq: %s (len = %i)\n", str, output_p->mapped_len_p[i]);
	  
	  memcpy(str, output_p->mapped_ref_p[i], output_p->mapped_len_p[i]);
	  str[output_p->mapped_len_p[i]] = '\0';
	  printf("ref: %s (len = %i, start = %i)\n", str, output_p->mapped_len_p[i], output_p->start_p[i]);
	  
	  printf("dif: ");
	  for(int j = 0; j < output_p->mapped_len_p[i]; j++) {
	       printf(output_p->mapped_seq_p[i][j] == output_p->mapped_ref_p[i][j] ? " " : "x");
	  }
	  
	  printf(" (score = %.2f, norm. score = %.2f)\n", output_p->score_p[i], output_p->norm_score_p[i]);
	  
	  printf("\n");
     }
}

//------------------------------------------------------------------------------------

sw_simd_context_t* sw_simd_context_new(float match, float mismatch, float gap_open, float gap_extend) {
  sw_simd_context_t* context_p = (sw_simd_context_t*) malloc(sizeof(sw_simd_context_t));
  //sw_simd_context_t* context_p = (sw_simd_context_t*) calloc(1, sizeof(sw_simd_context_t));

  context_p->gap_open = gap_open;
  context_p->gap_extend = gap_extend;
  /*
#ifdef __AVX__
  //printf("avx....size of zero (%i)\n", sizeof(context_p->zero));
  context_p->zero_simd = _mm256_setzero_ps();
  context_p->gap_open_simd = _mm256_set1_ps(gap_open);
  context_p->gap_extend_simd = _mm256_set1_ps(gap_extend);
  
  /*
  context_p->zero_p = (__m256 *) _mm_malloc(sizeof(__m256), SIMD_ALIGN);
  context_p->gap_open_p = (__m256 *) _mm_malloc(sizeof(__m256), SIMD_ALIGN);
  context_p->gap_extend_p = (__m256 *) _mm_malloc(sizeof(__m256), SIMD_ALIGN);

  *(context_p->zero_p) = _mm256_setzero_ps();
  *(context_p->gap_open_p) = _mm256_set1_ps(gap_open);
  *(context_p->gap_extend_p) = _mm256_set1_ps(gap_extend);
  
#else
  //printf("sse....size of zero (%i)\n", sizeof(context_p->zero));

     context_p->zero_simd = _mm_setzero_ps();
     context_p->gap_open_simd = _mm_set1_ps(gap_open);
     context_p->gap_extend_simd = _mm_set1_ps(gap_extend);
#endif // __AVX__     
*/
     //printf("mismatch = %0.2f, match = %0.2f\n", mismatch, match);

     context_p->substitution[0] = mismatch;
     context_p->substitution[1] = match;
          
     context_p->x_size = 0;
     context_p->y_size = 0;
     context_p->max_size = 0;
     
     context_p->E = NULL;
     context_p->F = NULL;
     context_p->H = NULL;
     context_p->C = NULL;

     context_p->compass_p = NULL;
     
     context_p->a_map = NULL;
     context_p->b_map = NULL;
     
     context_p->q_aux = NULL;
     context_p->r_aux = NULL;
     context_p->aux_size = 0;
     context_p->H_size = 0;
     context_p->F_size = 0; 

       // init matrix to -1000.0f
     for (int i = 0; i<128; i++) {
       for (int j = 0; j<128; j++) {
	 context_p->matrix[i][j] = -1000.0f;
       }
     }

     context_p->matrix['A']['A'] = 5.0f; context_p->matrix['C']['A'] = -4.0f; context_p->matrix['T']['A'] = -4.0f; context_p->matrix['G']['A'] = -4.0f;

     context_p->matrix['A']['C'] = -4.0f; context_p->matrix['C']['C'] = 5.0f; context_p->matrix['T']['C'] = -4.0f; context_p->matrix['G']['C'] = -4.0f;
     context_p->matrix['A']['G'] = -4.0f; context_p->matrix['C']['T'] = -4.0f; context_p->matrix['T']['T'] = 5.0f; context_p->matrix['G']['T'] = -4.0f;
     context_p->matrix['A']['T'] = -4.0f; context_p->matrix['C']['G'] = -4.0f; context_p->matrix['T']['G'] = -4.0f; context_p->matrix['G']['G'] = 5.0f;

     //     sw_simd_context_update(200, 800, context_p);


     return context_p;
}

//------------------------------------------------------------------------------------

void sw_simd_context_update(int x_size, int y_size, sw_simd_context_t* context_p) {
  //printf("sw_simd_context_update...\n");
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
     //printf("sw_simd_context_update done !!\n");
}

//------------------------------------------------------------------------------------

void sw_simd_context_free(sw_simd_context_t* context_p) {
  if (context_p == NULL) return;

  if (context_p->E != NULL) _mm_free(context_p->E);
  if (context_p->F != NULL) _mm_free(context_p->F);
  if (context_p->H != NULL) _mm_free(context_p->H);
  if (context_p->C != NULL) _mm_free(context_p->C);

  if (context_p->compass_p != NULL) free(context_p->compass_p);

  if (context_p->a_map != NULL) free(context_p->a_map);
  if (context_p->b_map != NULL) free(context_p->b_map);

  if (context_p->q_aux != NULL) free(context_p->q_aux);
  if (context_p->r_aux != NULL) free(context_p->r_aux);
  //free(context_p);
  free(context_p);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void smith_waterman_simd2(sw_simd_input_t* input_p, sw_simd_output_t* output_p, sw_simd_context_t* context_p) {
     int x_size = 0, y_size = 0;
     
     for (int i = 0; i < SIMD_DEPTH; i++) {
	  output_p->score_p[i] = 0;
	  if (input_p->seq_len_p[i] > x_size) x_size = input_p->seq_len_p[i];
	  if (input_p->ref_len_p[i] > y_size) y_size = input_p->ref_len_p[i];
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
     
     float max_penalty = context_p->substitution[0] * 100.0;
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
	       for (int i = 0; i < SIMD_DEPTH; i++) {
		    if (x>input_p->seq_len_p[i] || y>input_p->ref_len_p[i]) {
			 ((float*) &substitution_simd)[i] = max_penalty;
		    } else {
			 seq = input_p->seq_p[i][x-1];
			 ref = input_p->ref_p[i][y-1];
			 ((float*) &substitution_simd)[i] = context_p->substitution[seq == ref];
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
	       //printf("h calculation done!!\ns calculation...\n");
	       // S (maximum score)
	       for (int i = 0; i < SIMD_DEPTH; i++) {
		    compass_p[i] = 3;
		    if (((float *) &h_final_simd)[i] == ((float *) &h_value_simd)[i]) {
		      compass_p[i] = 0;
		    } else if (((float *) &h_final_simd)[i] == ((float *) &f_value_simd)[i]) {
		      compass_p[i] = 1;
		    } else if (((float *) &h_final_simd)[i] == ((float *) &e_value_simd)[i]) {
		      compass_p[i] = 2;
		    }

		    if (((float *) &h_final_simd)[i] > output_p->score_p[i]) {
			 output_p->score_p[i] = ((float *) &h_final_simd)[i];;
			 seq_end_p[i] = &input_p->seq_p[i][x-1];
			 ref_end_p[i] = &input_p->ref_p[i][y-1];
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
#ifdef __AVX__
	       _mm256_store_ps(h, h_final_simd);
	       _mm256_store_ps(e, e_value_simd);
	       _mm256_store_ps(f, f_value_simd);
#else
	       _mm_store_ps(h, h_final_simd);
	       _mm_store_ps(e, e_value_simd);
	       _mm_store_ps(f, f_value_simd);
#endif // __AVX__
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
     
     float score, *score_p, *aux_score_p; //diag, up, left;
     int compass;
     
     int x_end, num_gaps; 
     float bimble, errbounds = 0.01f;
     float gap_open = context_p->gap_open;
     float gap_extend = context_p->gap_extend;
     //printf("walking through the score matrix H\n");
     //printf("address: compass_p %x, end compass %x, max_size end compass = %x\n", compass_p, compass_end_p[0], (int*) &(context_p->compass_p + (max_size * SIMD_DEPTH))[0]);
     //printf("value: compass_p %i, end compass %i, max_size end compass = %i\n", *compass_p, *compass_end_p[0], (context_p->compass_p + (max_size * SIMD_DEPTH))[0]);
     /*
     int *p = compass_end_p[0];
     printf("compass end: %x, %i\n", p, *p);
     p -= x_size_1_depth; // SIMD_DEPTH; // x_size_depth; //SIMD_DEPTH;
     printf("compass end: %x, %i\n", p, *p);
     p -= x_size_1_depth; // SIMD_DEPTH; // x_size_depth; //SIMD_DEPTH;
     printf("compass end: %x, %i\n", p, *p);
     p -= x_size_1_depth; // SIMD_DEPTH; // x_size_depth; //SIMD_DEPTH;
     printf("compass end: %x, %i\n", p, *p);
     p -= x_size_1_depth; // SIMD_DEPTH; // x_size_depth; //SIMD_DEPTH;
     printf("compass end: %x, %i\n", p, *p);
     p -= x_size_1_depth; // SIMD_DEPTH; // x_size_depth; //SIMD_DEPTH;
     printf("compass end: %x, %i\n", p, *p);
     p -= x_size_1_depth; // SIMD_DEPTH; // x_size_depth; //SIMD_DEPTH;
     printf("compass end: %x, %i\n", p, *p);
     p -= x_size_1_depth; // SIMD_DEPTH; // x_size_depth; //SIMD_DEPTH;
     printf("compass end: %x, %i\n", p, *p);
     p -= x_size_1_depth; // SIMD_DEPTH; // x_size_depth; //SIMD_DEPTH;
     printf("compass end: %x, %i\n", p, *p);
     //     printf("compass end (x-1): %x, %i\n", (compass_end_p[0]-SIMD_DEPTH), *((compass_end_p[0]-SIMD_DEPTH)[0]));
     */
     // walking the score matrix (H)
     for (int i = 0; i < SIMD_DEPTH; i++) {
	  a_mapp = a_map + max_size;
	  b_mapp = b_map + max_size;
	  
	  *a_mapp = 0;
	  //a_mapp--; *a_mapp = *seq_end_p[i];
	  
	  *b_mapp = 0;
	  //b_mapp--; *b_mapp = *ref_end_p[i];
	  
	  //x_end = seq_x_end_p[i];
	  output_p->start_p[i] = ref_y_end_p[i] + 1;
	  
	  score_p = h_end_p[i];
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
	      score_p -= x_size_1_depth;
	      compass_p -= x_size_1_depth;
	      a_mapp--; *a_mapp = *(seq_end_p[i]--);
	      b_mapp--; *b_mapp = *(ref_end_p[i]--);
	      output_p->start_p[i]--;
	      
	      if (*score_p <= 0) break;
	      //continue;
	      //if (--x_end <= 0) break;

	    } else if (compass == 2) {
	      //printf("cipf: down\n");
	      num_gaps = 0;
	      score = *score_p;
	      aux_score_p = score_p - x_size_depth;

	      while (1) {
		bimble = *aux_score_p - gap_open - (num_gaps * gap_extend);
		//printf("down: fabs(score (%0.2f) -  bimble (%0.2f)) = %0.2f, errbounds = %0.2f\n", score, bimble, fabs(score - bimble), errbounds);

		if (aux_score_p < &H[i] || fabs(score - bimble) < errbounds) break;
		aux_score_p -= x_size_depth;
		++num_gaps;
	      }
	      if (bimble <= 0.0) break;

	      //if (num_gaps > 0) printf("down: num gaps = %i\n", num_gaps);

	      for (int j = 0; j <= num_gaps ; j++) {
		score_p -= x_size_depth;
		compass_p -= x_size_depth;
		a_mapp--; *a_mapp = '-';
		b_mapp--; *b_mapp = *(--ref_end_p[i]);
		output_p->start_p[i]--;
	      }
	    } else {
	      //printf("cipf: left\n");
	      num_gaps = 0;
	      score = *score_p;
	      aux_score_p = score_p - SIMD_DEPTH;
	      //printf("left: aux_score_p = %li, &H[%i] = %li\n", aux_score_p, i, &H[i]);

	      while (1) {
		bimble = *aux_score_p - gap_open - (num_gaps * gap_extend);
		//printf("left: aux_score_p = %f, fabs(score (%f) -  bimble (%f)) = %f, errbounds = %f, condition = %i\n", *aux_score_p, score, bimble, fabs(score - bimble), errbounds, fabs(score - bimble) < errbounds);

		//printf("aux_score_p = %li, &H[%i] = %li\n", aux_score_p, i, &H[i]);
		if (aux_score_p < &H[i] || fabs(score - bimble) < errbounds) break;

		aux_score_p -= SIMD_DEPTH;
		++num_gaps;
	      }
	      if (bimble <= 0.0) break;

	      //if (num_gaps > 0) 
		//printf("----------left: num gaps = %i\n", num_gaps);
	      //exit(-1);

	      for (int j = 0; j <= num_gaps ; j++) {
		score_p -= SIMD_DEPTH;
		compass_p -= SIMD_DEPTH;
		a_mapp--; *a_mapp = *(seq_end_p[i]--);
		b_mapp--; *b_mapp = '-';
	      }
	      //if (count++ > 1) break;
	      //break;
	      //continue;
	      //if (--x_end <= 0) break;
	    }
	  } // end of for ;;
	  
	  //output_p->mapped_seq_p[i] = strdup(a_mapp);
	  output_p->mapped_len_p[i] = a_map + max_size - a_mapp;
	  output_p->mapped_seq_p[i] = (char*) calloc(output_p->mapped_len_p[i], sizeof(char));
	  memcpy(output_p->mapped_seq_p[i], a_mapp, output_p->mapped_len_p[i]);
	  
	  //output_p->mapped_ref_p[i] = strdup(b_mapp);
	  //output_p->mapped_ref_len_p[i] = b_map + max_size - b_mapp;
	  output_p->mapped_ref_p[i] = (char*) calloc(output_p->mapped_len_p[i], sizeof(char));
	  memcpy(output_p->mapped_ref_p[i], b_mapp, output_p->mapped_len_p[i]);
	  
	  output_p->norm_score_p[i] = output_p->score_p[i] / (output_p->mapped_len_p[i] * context_p->substitution[1]);
     } // end of for 0..SIMD_DEPTH
}

//-------------------------------------------------------------

void smith_waterman_simd(sw_simd_input_t* input, sw_simd_output_t* output, 
			 sw_simd_context_t* context) {

  int len, simd_depth = 4, num_queries = input->depth;
  int max_q_len = 0, max_r_len = 0;

  //float gap_open = context->gap_open;
  //float gap_extend = context->gap_extend;

  //  printf("Process %d reads\n", num_queries);

  for (size_t i = 0; i < num_queries; i++) {
    if (input->seq_len_p[i] > max_q_len) max_q_len = input->seq_len_p[i];
    if (input->ref_len_p[i] > max_r_len) max_r_len = input->ref_len_p[i];
  }

  reallocate_memory(max_q_len, max_r_len, simd_depth, 
		    &context->H_size, &context->H, &context->C, 
		    &context->F_size, &context->F, 
		    &context->aux_size, &context->q_aux, &context->r_aux);

  sse1_matrix(num_queries, 
	      input->seq_p, input->seq_len_p, max_q_len, 
	      input->ref_p, input->ref_len_p, max_r_len, 
	      context->matrix, context->gap_open, context->gap_extend, 
	      context->H, context->F, context->C, output->score_p);
  
  simd_traceback1(simd_depth, num_queries, 
		  input->seq_p, input->seq_len_p, max_q_len, 
		  input->ref_p, input->ref_len_p, max_r_len, 
		  context->gap_open, context->gap_extend, 
		  context->H, context->C, output->score_p,
		  output->mapped_seq_p, output->start_seq_p,
		  output->mapped_ref_p, output->start_p, 
		  output->mapped_len_p, 
		  context->q_aux, context->r_aux);
  
  for(int i = 0; i < simd_depth; i++){
    //printf("Out: %s\n", output->mapped_seq_p[i]);
    output->norm_score_p[i] = output->score_p[i] / (output->mapped_len_p[i] * context->substitution[1]);
    //printf("Output Score %d, %d, %d\n", output->norm_score_p[i], output->mapped_len_p[i], output->score_p[i]);
  }
}

//====================================================================================
// Smith-Waterman functions from EMBOSS package
//====================================================================================

float smith_waterman(char* seq_a, char* seq_b, float gapopen, float gapextend, 
		    char* m, char* n, int* start1, int* start2) {

  int len_a = strlen(seq_a);
  int len_b = strlen(seq_b);
  int total = (len_a + 1) * (len_b + 1);

  float score = 0;
  float* path = (float*) calloc(total, sizeof(float));
  int* compass = (int*) calloc(total, sizeof(int));

  score = AlignPathCalcSW(seq_a, seq_b, len_a, len_b, gapopen, gapextend, path, compass);


  //  saveFloatMatrixByIndex0(path, len_a, len_b, 0, 1, seq_a, seq_b, "/tmp/emboss-matrix.txt");
  //  saveIntegerMatrixByIndex0(compass, len_a, len_b, 0, 1, seq_a, seq_b, "/tmp/emboss-compass.txt");

  AlignWalkSWMatrix(path, compass, gapopen, gapextend, seq_a, seq_b,
		    (char*) m, (char *) n, strlen(seq_a), strlen(seq_b), start1, start2);

  free(path);
  free(compass);

  return score;
}

//-------------------------------------------------------------

/* @func embAlignPathCalcSW ***************************************************
**
** Create path matrix for Smith-Waterman
** Nucleotides or proteins as needed.
**
** @param [r] a [const char *] first sequence
** @param [r] b [const char *] second sequence
** @param [r] lena [ajint] length of first sequence
** @param [r] lenb [ajint] length of second sequence
** @param [r] gapopen [float] gap opening penalty
** @param [r] gapextend [float] gap extension penalty
** @param [w] path [float *] path matrix
** @param [w] compass [ajint *] Path direction pointer array
**
** @return [float] Maximum score
** @@
** Optimised to keep a maximum value to avoid looping down or left
** to find the maximum. (il 29/07/99)
******************************************************************************/

float AlignPathCalcSW(const char *a, const char *b, int lena, int lenb,
                      float gapopen, float gapextend, float* path,
                      int* compass)
{
    float ret;
    long xpos;
    long ypos;
    long i;
    long j;

    double match;
    double mscore;
    double result;
    double fnew;
    double* maxa;

    double bx;
    char compasschar;

    ret= -FLT_MAX;

    /* Create stores for the maximum values in a row or column */

    maxa = (double*) calloc(lena, sizeof(double));

    /* First initialise the first column and row */
    for(i=0;i<lena;++i)
    {
        result = (a[i]==b[0] ? 5.0 : -4.0);

	fnew = i==0 ? 0. :
		path[(i-1)*lenb] -(compass[(i-1)*lenb]==DOWN ?
			gapextend : gapopen);

	if (result > fnew && result>0)
	{
	  path[i*lenb] = (float) result;
	  compass[i*lenb] = 0;
	}
	else if (fnew>0)
	{
	  path[i*lenb] = (float) fnew;
	  compass[i*lenb] = DOWN;
	}
	else
	{
	    path[i*lenb] = 0.;
	    compass[i*lenb] = 3;
	}

	maxa[i] = i==0 ? path[i*lenb]-gapopen :
	path[i*lenb] - (compass[(i-1)*lenb]==DOWN ? gapextend : gapopen);
    }

    for(j=0;j<lenb;++j)
    {
        result = (a[0]==b[j] ? 5.0 : -4.0);

	fnew = j==0 ? 0. :
		path[j-1] -(compass[j-1]==LEFT ? gapextend : gapopen);

	if (result > fnew && result > 0)
	{
	  path[j] = (float) result;
	    compass[j] = 0;
	}
	else if (fnew >0)
	{
	  path[j] = (float) fnew;
	    compass[j] = LEFT;
	}
	else
	{
	    path[j] = 0.;
	    compass[j] = 3;
	}

    }


    /* xpos and ypos are the diagonal steps so start at 1 */
    xpos = 1;
    float aux;

    while(xpos!=lenb)
    {
	ypos  = 1;
	bx = path[xpos]-gapopen-gapextend;

	while(ypos < lena)
	{
	    /* get match for current xpos/ypos */
            match = (a[ypos]==b[xpos] ? 5.0 : -4.0);

	    /* Get diag score */
	    mscore = path[(ypos-1)*lenb+xpos-1] + match;
	    aux = mscore;

	    /* Set compass to diagonal value 0 */
	    compass[ypos*lenb+xpos] = 0;
	    path[ypos*lenb+xpos] = (float) mscore;


	    /* Now parade back along X axis */
            maxa[ypos] -= gapextend;
            fnew=path[(ypos)*lenb+xpos-1];
            fnew-=gapopen;

            if(fnew > maxa[ypos])
                maxa[ypos] = fnew;

            if( maxa[ypos] > mscore)
            {
                mscore = maxa[ypos];
                path[ypos*lenb+xpos] = (float) mscore;
                compass[ypos*lenb+xpos] = LEFT; /* Score comes from left */
            }

	    /* And then bimble down Y axis */
            bx -= gapextend;
            fnew = path[(ypos-1)*lenb+xpos];
            fnew-=gapopen;

            if(fnew > bx)
                bx = fnew;

            if(bx > mscore)
            {
                mscore = bx;
                path[ypos*lenb+xpos] = (float) mscore;
                compass[ypos*lenb+xpos] = DOWN; /* Score comes from bottom */
            }

	    /*
	    if (ypos == xpos) {
	      printf("(%i, %i):\tmscore = %0.2f (%0.2f)\tmaxa[%i] = %0.2f\tbx = %0.2f\t-> compass = %i\n", xpos, ypos, mscore, aux, ypos, maxa[ypos], bx, compass[ypos*lenb+xpos]);
	    }
	    */

            if(mscore > ret)
                ret = (float) mscore;

	    result = path[ypos*lenb+xpos];
	    if(result < 0.) {
		path[ypos*lenb+xpos] = 0.;
		compass[ypos*lenb+xpos] = 3;
	    }

	    ypos++;
	}
	++xpos;
    }

    //printf("before free maxa: %x\n", maxa);
    free(maxa);

    return ret;
}

/******************************************************************************/

void AlignWalkSWMatrix(const float* path, const int* compass,
		       float gapopen, float gapextend,
		       const char*  a, const char* b,
		       char* m, char* n,
		       int lena, int lenb,
		       int *start1, int *start2)
{
  long i;
  long j;
  long k;
  long gapcnt;
  double pmax;
  double score;
  double bimble;

  long ix;
  long iy;
  
  long xpos = 0;
  long ypos = 0;
  const char *p;
  const char *q;
  
  int ic;
  double errbounds;

  /* errbounds = gapextend; */
  errbounds = (double) 0.01;

  /* Get maximum path score and save position */
  pmax = -FLT_MAX;
  k = (long)lena*(long)lenb-1;

  for(i=lena-1; i>=0; --i)
    for(j=lenb-1; j>=0; --j)
      if((path[k--] > pmax) || E_FPEQ(path[k+1],pmax,U_FEPS))
	{
	    pmax = path[k+1];
	    xpos = j;
	    ypos = i;
	  }

  //ajStrAssignClear(m);
  //ajStrAssignClear(n);

  p = a;
  q = b;

  while(xpos>=0 && ypos>=0)
    {
      if(!compass[ypos*lenb+xpos])    /* diagonal */
        {
	  //printf("emboss: diagonal\n");
	  //ajStrAppendK(m,p[ypos--]);
	  //ajStrAppendK(n,q[xpos--]);

	  strncat(m, &p[ypos--], 1);
	  strncat(n, &q[xpos--], 1);

	  if(ypos >= 0 && xpos>=0 && path[(ypos)*lenb+xpos]<=0.)
	    break;

	  continue;
        }
      else if(compass[ypos*lenb+xpos]==LEFT) /* Left, gap(s) in vertical */
        {
	  //printf("emboss: left\n");
	  score  = path[ypos*lenb+xpos];
	  gapcnt = 0;
	  ix     = xpos-1;

	  while(1)
            {
	      bimble = path[ypos*lenb+ix]-gapopen-(gapcnt*gapextend);

	      if(!ix || fabs((double)score-(double)bimble)<errbounds)
		break;

	      --ix;
	      ++gapcnt;
            }

	  if(bimble<=0.0)
	    break;

	  //printf("LEFT: gapcnt = %i\n", gapcnt);
	  for(ic=0;ic<=gapcnt;++ic)
            {
	      //ajStrAppendK(m,'.');
	      //ajStrAppendK(n,q[xpos--]);
	      strcat(m, "-");
	      strncat(n, &q[xpos--], 1);
            }

	  continue;
        }
      else if(compass[ypos*lenb+xpos]==DOWN) /* Down, gap(s) in horizontal */
        {
	  //printf("emboss: down\n");
	  score  = path[ypos*lenb+xpos];
	  gapcnt = 0;
	  iy = ypos-1;

	  while(1)
            {
	      bimble=path[iy*lenb+xpos]-gapopen-(gapcnt*gapextend);

	      if(!iy || fabs((double)score-(double)bimble)<errbounds)
		break;

              --iy;

	      if(iy<0) {
		printf("SW: Error walking down");
		exit(-1);
	      }

	      ++gapcnt;
            }

	  if(bimble<=0.0)
	    break;

	  //printf("DOWN: gapcnt = %i\n", gapcnt);
	  for(ic=0;ic<=gapcnt;++ic)
            {
	      //ajStrAppendK(m,p[ypos--]);
	      //ajStrAppendK(n,'.');
	      strncat(m, &p[ypos--], 1);
	      strcat(n, "-");
            }
	  continue;
        }
      else {
	printf("Walk Error in SW");
	exit(-1);
      }
    }

  *start1 = (int) (ypos + 1); /* Potential lossy cast */
  *start2 = (int) (xpos + 1); /* Potential lossy cast */

  //ajStrReverse(m);            /* written with append, need to reverse */
  //ajStrReverse(n);
  
  revstr(m);
  revstr(n);

  return;
}


void revstr(char* str) {
  int i;

  int len = strlen(str);
  char cpstr[len+1];

  for(i=0; i < len ; i++) {
    cpstr[i] = str[len-i-1];
  }
  cpstr[i] = '\0';

  strcpy(str, cpstr);
}


//-------------------------------------------------------------
//-------------------------------------------------------------


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

//-------------------------------------------------------------
//-------------------------------------------------------------
