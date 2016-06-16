#include "sw_commons.h"

extern const unsigned int simd_depth;
extern const unsigned int aligned_mem;

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void init_subst_score_matrix(char *filename, subst_matrix_t matrix) {
    FILE *file = fopen(filename, "r");

    if (file == NULL) {
        printf("Error: substitution score matrix file (%s) not found\n", filename);
        exit(-1);
    }

    char *header[256],  *token[256];

    char *header_line = (char*) calloc(1, 4096);
    char *token_line = (char*) calloc(1, 4096);

    char *res = fgets(header_line, 4096, file);
    str_trim(header_line);


    res = NULL;

    // init matrix to -1000.0f
    for (int i = 0; i<128; i++) {
        for (int j = 0; j<128; j++) {
            matrix[i][j] = -1000.0f;
        }
    }
    /*
    matrix['A']['A'] = 5.0f; matrix['C']['A'] = -4.0f; matrix['T']['A'] = -4.0f; matrix['G']['A'] = -4.0f;
    matrix['A']['C'] = -4.0f; matrix['C']['C'] = 5.0f; matrix['T']['C'] = -4.0f; matrix['G']['C'] = -4.0f;
    matrix['A']['G'] = -4.0f; matrix['C']['T'] = -4.0f; matrix['T']['T'] = 5.0f; matrix['G']['T'] = -4.0f;
    matrix['A']['T'] = -4.0f; matrix['C']['G'] = -4.0f; matrix['T']['G'] = -4.0f; matrix['G']['G'] = 5.0f;
    */

    // read header row
    unsigned int num_columns = 0;

    header[num_columns] = strtok(header_line, "\t");
    //  res = strtok(header_line, "\t");
    while (header[num_columns]!= NULL) {
        num_columns++;
        header[num_columns] = strtok(NULL, "\t");
    }

    // read the remain rows and update matrix
    unsigned int col = 0;
    while (fgets(token_line, 4096, file) != NULL) {
        str_trim(token_line);
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
            matrix[(unsigned char)header[col][0]][(unsigned char)token[0][0]] = atof(token[col]);
        }
    }

    free(header_line);
    free(token_line);

    fclose(file);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

sw_optarg_t* sw_optarg_new(float match, float mismatch,
                           float gap_open, float gap_extend, char *subst_matrix_name) {
    sw_optarg_t *optarg_p = calloc(1, sizeof(sw_optarg_t));

    optarg_p->match = match;
    optarg_p->mismatch = mismatch;
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

    //output_p->score_p = (float*) calloc(num_queries, sizeof(float));
    output_p->score_p = (float *) _mm_malloc(num_queries * sizeof(float), aligned_mem);

    return output_p;
}

//------------------------------------------------------------------------------------

void sw_multi_output_free(sw_multi_output_t* output_p) {
    if (output_p == NULL) return;

    for (unsigned int i = 0; i < output_p->num_queries; i++) {
        if (output_p->query_map_p[i] != NULL) free(output_p->query_map_p[i]);
        if (output_p->ref_map_p[i] != NULL) free(output_p->ref_map_p[i]);
    }
    free(output_p->query_map_p);
    free(output_p->ref_map_p);

    if (output_p->query_map_len_p != NULL) free(output_p->query_map_len_p);
    if (output_p->ref_map_len_p != NULL) free(output_p->ref_map_len_p);

    if (output_p->query_start_p != NULL) free(output_p->query_start_p);
    if (output_p->ref_start_p != NULL) free(output_p->ref_start_p);

    //if (output_p->score_p != NULL) free(output_p->score_p);
    if (output_p->score_p != NULL) _mm_free(output_p->score_p);

    free(output_p);
}

//------------------------------------------------------------------------------------

void sw_multi_output_save(int num_alignments, sw_multi_output_t* output_p, FILE *file_p) {
    unsigned int len, identity, gaps;

    if (file_p == NULL) {
        file_p = stdout;
    }

    for(int i = 0; i < num_alignments; i++) {
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
                fprintf(file_p, "x");
            }
        }
        fprintf(file_p, "\n");
        fprintf(file_p, "Ref. : %s\tStart at %i\n", output_p->ref_map_p[i], output_p->ref_start_p[i]);

        fprintf(file_p, "Score: %.2f\tLength: %i\tIdentity: %.2f\tGaps: %.2f\n",
                output_p->score_p[i], len, identity * 100.0f / len, gaps * 100.0f / len);

        fprintf(file_p, "\n");
    }
}

//------------------------------------------------------------------------------------

int sw_multi_output_string(int num_alignments, sw_multi_output_t* output_p, char *buf) {
  int total, len, identity, gaps;
  total = 0;

  for(int i = 0; i < num_alignments; i++) {
    gaps = 0;
    identity = 0;
    len = strlen(output_p->query_map_p[i]);

    total += sprintf(buf + total, "Query: %s\tStart at %i\n", output_p->query_map_p[i], output_p->query_start_p[i]);
    total += sprintf(buf + total,"       ");
    for(int j = 0; j < len; j++) {
      if (output_p->query_map_p[i][j] == '-' || output_p->ref_map_p[i][j] == '-') {
	gaps++;
      }
      if (output_p->query_map_p[i][j] == output_p->ref_map_p[i][j]) {
	total += sprintf(buf + total, "|");
	identity++;
      } else {
	total += sprintf(buf + total, "x");
      }
    }
    total += sprintf(buf + total, "\n");
    total += sprintf(buf + total, "Ref. : %s\tStart at %i\n", output_p->ref_map_p[i], output_p->ref_start_p[i]);
    
    total += sprintf(buf + total, "Score: %.2f\tLength: %i\tIdentity: %.2f\tGaps: %.2f\n",
		     output_p->score_p[i], len, identity * 100.0f / len, gaps * 100.0f / len);
    
    total += sprintf(buf + total, "\n");
  }
  return total;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
