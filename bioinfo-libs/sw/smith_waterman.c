#include "smith_waterman.h"

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

sw_optarg_t* sw_optarg_new(float gap_open, float gap_extend, char *subst_matrix_name) {
    sw_optarg_t *optarg_p = calloc(1, sizeof(sw_optarg_t));

    optarg_p->gap_open = gap_open;
    optarg_p->gap_extend = gap_extend;
    
    init_subst_score_matrix(subst_matrix_name, optarg_p->subst_matrix);

    return optarg_p;
}

//------------------------------------------------------------------------------------

void sw_optarg_free(sw_optarg_t* optarg_p) {
    free(optarg_p);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

sw_multi_output_t* sw_multi_output_new(unsigned int num_queries) {
    sw_multi_output_t *output_p = calloc(1, sizeof(sw_multi_output_t));
  
    output_p->num_queries = num_queries;
    output_p->query_map_p = (char**) calloc(num_queries, sizeof(char*));
    output_p->ref_map_p = (char**) calloc(num_queries, sizeof(char*));

    output_p->query_map_len_p = (unsigned int*) calloc(num_queries, sizeof(unsigned int));
    output_p->ref_map_len_p = (unsigned int*) calloc(num_queries, sizeof(unsigned int));
    
    output_p->query_start_p = (unsigned int*) calloc(num_queries, sizeof(unsigned int));
    output_p->ref_start_p = (unsigned int*) calloc(num_queries, sizeof(unsigned int));
    
    output_p->score_p = (float*) calloc(num_queries, sizeof(float));
    
    return output_p;
}

//------------------------------------------------------------------------------------

void sw_multi_output_free(sw_multi_output_t* output_p) {
    if (output_p == NULL) return;
     
    if (output_p->query_map_p != NULL) free(output_p->query_map_p);
    if (output_p->ref_map_p != NULL) free(output_p->ref_map_p);

    if (output_p->query_map_len_p != NULL) free(output_p->query_map_len_p);
    if (output_p->ref_map_len_p != NULL) free(output_p->ref_map_len_p);
    
    if (output_p->query_start_p != NULL) free(output_p->query_start_p);
    if (output_p->ref_start_p != NULL) free(output_p->ref_start_p);
    
    if (output_p->score_p != NULL) free(output_p->score_p);
    
    free(output_p);
}

//------------------------------------------------------------------------------------

void sw_multi_output_save(sw_multi_output_t* output_p, FILE *file_p) {
    unsigned int num_queries = output_p->num_queries;
    unsigned int len, identity, gaps;

    if (file_p == NULL) {
      file_p = stdout;
    }
    
    for(int i = 0; i < num_queries; i++) {
        gaps = 0;
        identity = 0;
        len = strlen(output_p->query_map_p[i]);

	fprintf(file_p, "Query: %s\tStart at %i\n", output_p->query_map_p[i], output_p->query_start_p[i]);
	fprintf(file_p,"       ");
	for(int j = 0; j < len; j++) {
	  if (output_p->query_map_p[i][j] == '-' || output_p->ref_map_p[i][j] == '-') {
	    gaps++;
	  }
	  if (output_p->query_map_p[i][j] == output_p->ref_map_p[i][j]) {
	    fprintf(file_p, "|");
	    identity++;
	  } else {
	    fprintf(file_p, " ");
	  }
	}
	fprintf(file_p, "\n");
	fprintf(file_p, "Ref. : %s\tStart at %i\n", output_p->ref_map_p[i], output_p->ref_start_p[i]);
	fprintf(file_p, "Score: %.2f\tLength: %i\tIdentity: %.2f%\tGaps: %.2f%\n", 
		output_p->score_p[i], len, identity * 100.0f / len, gaps * 100.0f / len);
	       
	fprintf(file_p, "\n");
    }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void smith_waterman_mqmr(char **query_p, char **ref_p, unsigned int num_seqs, 
			 sw_optarg_t *optarg_p, unsigned int num_threads, 
			 sw_multi_output_t *output_p) {

    if (output_p == NULL) {
        printf("Error: output buffer is null.\n");
	exit(-1);
    }
    if (output_p->num_queries != num_seqs) {
        printf("Error: num. sequencias (%i) and output buffer length (%i) mismatch !\n", 
	       num_seqs, output_p->num_queries);
	exit(-1);
    }

    unsigned int num_packs = num_seqs / SIMD_DEPTH;
    if (num_seqs % SIMD_DEPTH) {
        num_packs++;
    }
    unsigned int packs_per_thread = num_packs / num_threads;
    if (num_packs % num_threads) {
        packs_per_thread++;
    }

    #pragma omp parallel num_threads(num_threads)
    {
        unsigned int tid = omp_get_thread_num();

	sw_simd_context_t *context_p = sw_simd_context_new(optarg_p->gap_open, optarg_p->gap_extend,
							   &(optarg_p->subst_matrix));
							   
	unsigned int i;
	unsigned int first_index = tid * packs_per_thread * SIMD_DEPTH;
	unsigned int last_index = first_index + ((packs_per_thread - 1) * SIMD_DEPTH);
	if (last_index > num_seqs) last_index = num_seqs;

	for (i = first_index; i < last_index; i += SIMD_DEPTH) {
  	    smith_waterman_simd(&query_p[i], &ref_p[i], SIMD_DEPTH,
				&(output_p->query_map_p[i]), &(output_p->ref_map_p[i]), 
				&(output_p->query_start_p[i]), &(output_p->ref_start_p[i]),
				&(output_p->score_p[i]), context_p);
	}

	int num_queries = (i + SIMD_DEPTH > num_seqs) ? (num_seqs - i) : SIMD_DEPTH;
	if (num_queries > 0) {
	    smith_waterman_simd(&query_p[i], &ref_p[i], num_queries,
				&(output_p->query_map_p[i]), &(output_p->ref_map_p[i]), 
				&(output_p->query_start_p[i]), &(output_p->ref_start_p[i]),
				&(output_p->score_p[i]), context_p);
	}

	free(context_p);
    }
}

//------------------------------------------------------------------------------------

void smith_waterman_mqsr(char **query_p, char *ref_p, unsigned int num_seqs, 
			 sw_optarg_t *optarg_p, unsigned int num_threads, 
			 sw_multi_output_t *output_p) {
  
    if (output_p == NULL) {
        printf("Error: output buffer is null.\n");
	exit(-1);
    }
    if (output_p->num_queries != num_seqs) {
        printf("Error: num. sequencias (%i) and output buffer length (%i) mismatch !\n", 
	       num_seqs, output_p->num_queries);
	exit(-1);
    }

    unsigned int num_packs = num_seqs / SIMD_DEPTH;
    if (num_seqs % SIMD_DEPTH) {
        num_packs++;
    }
    unsigned int packs_per_thread = num_packs / num_threads;
    if (num_packs % num_threads) {
        packs_per_thread++;
    }

    char* refs[SIMD_DEPTH];
    for (int i = 0; i < SIMD_DEPTH; i++) {
        refs[i] = ref_p;
    }

    #pragma omp parallel num_threads(num_threads)
    {
        unsigned int tid = omp_get_thread_num();

	sw_simd_context_t *context_p = sw_simd_context_new(optarg_p->gap_open, optarg_p->gap_extend,
							   &(optarg_p->subst_matrix));
							   
	unsigned int i;
	unsigned int first_index = tid * packs_per_thread * SIMD_DEPTH;
	unsigned int last_index = first_index + ((packs_per_thread - 1) * SIMD_DEPTH);
	if (last_index > num_seqs) last_index = num_seqs;

	for (i = first_index; i < last_index; i += SIMD_DEPTH) {
  	    smith_waterman_simd(&query_p[i], refs, SIMD_DEPTH,
				&(output_p->query_map_p[i]), &(output_p->ref_map_p[i]), 
				&(output_p->query_start_p[i]), &(output_p->ref_start_p[i]),
				&(output_p->score_p[i]), context_p);
	}

	int num_queries = (i + SIMD_DEPTH > num_seqs) ? (num_seqs - i) : SIMD_DEPTH;
	if (num_queries > 0) {
	    smith_waterman_simd(&query_p[i], refs, num_queries,
				&(output_p->query_map_p[i]), &(output_p->ref_map_p[i]), 
				&(output_p->query_start_p[i]), &(output_p->ref_start_p[i]),
				&(output_p->score_p[i]), context_p);
	}

	free(context_p);
    }
}

//------------------------------------------------------------------------------------

void init_subst_score_matrix(char *filename, subst_matrix_t matrix) {
  FILE *file = fopen(filename, "r");

  if (file == NULL) {
    printf("Error: substitution score matrix file (%s) not found\n", filename);
    exit(-1);
  }

  char header_line[4096], *header[256];
  char token_line[4096], *token[256];
  fgets(header_line, 4096, file);
  trim(header_line);

  // init matrix to -1000.0f
  for (int i = 0; i<128; i++) {
    for (int j = 0; j<128; j++) {
      matrix[i][j] = -1000.0f;
    }
  }

  // read header row
  unsigned int num_columns = 0;
  header[num_columns] = strtok(header_line, "\t");
  while (header[num_columns]!= NULL) {
    num_columns++;
    header[num_columns] = strtok(NULL, "\t");
  }

  // read the remain rows and update matrix
  unsigned int col = 0;
  while (fgets(token_line, 4096, file) != NULL) {
    trim(token_line);
    col = 0;
    token[col] = strtok(token_line, "\t");
    while (token[col]!= NULL) {
      col++;
      token[col] = strtok(NULL, "\t");
    }
    if (col != num_columns) {
      printf("Error: substitution score matrix invalid format\n");
      exit(-1);
    }

    for(col = 1; col < num_columns; col++) {
      matrix[header[col][0]][token[0][0]] = atof(token[col]);
    }
  }

  fclose(file);
}

//----------------------------------------------------------------                                               //----------------------------------------------------------------                                                 
