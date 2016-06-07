
#include "macros.h"
#include "sw_commons.h"
#include "smith_waterman.h"

#include "sw_omp.h"

//====================================================================
// O M P     F O R      H P G - S W
//====================================================================

void run_sw_omp(char *q_filename, char *r_filename,
		float match, float mismatch, float gap_open, float gap_extend,
		char *matrix_filename, int batch_size, int num_threads,
		char *out_filename) {

  printf("------------------------------------------\n");
  printf("Running function %s\n", __FUNCTION__);
  printf("------------------------------------------\n");

#ifdef SW_AVX2
  const unsigned int depth = 8;
#else
  const unsigned int depth = 4;
#endif

  const int max_length = 2048;

  FILE *q_file = fopen(q_filename, "r");
  FILE *r_file = fopen(r_filename, "r");
  FILE *out_file = fopen(out_filename, "w");

#ifdef TIMING
  double max_sse = 0.0f;
#endif

  #pragma omp parallel num_threads(num_threads)
  {
    double start_read = 0.0f, read = 0.0f;
    double start_write = 0.0f, write = 0.0f;
    double start_memory = 0.0f, memory = 0.0f;
    double sse_t = 0.0f, partial_t = 0.0f;

    sw_optarg_t *optarg_p = sw_optarg_new(match, mismatch, gap_open, gap_extend, matrix_filename);
    sw_multi_output_t *output_p;

    start_memory = sw_tic();
    char *q[batch_size], *r[batch_size];
    for (int i = 0; i < batch_size; i++) {
      q[i] = (char *) calloc(max_length, sizeof(char));
      r[i] = (char *) calloc(max_length, sizeof(char));
    }
    memory += sw_toc(start_memory);

    int count = 0, batches = 0, num_queries;

    while (1) {
      num_queries = 0;

      #pragma omp critical
      {
	start_read = sw_tic();
	// read queries
	for (int i = 0; i < batch_size; i++) {
	  if (fgets(q[i], max_length, q_file) == NULL) { break; }
	  str_trim(q[i]);
	  num_queries++;
	  count++;
	}

	// read references
	for (int i = 0; i < num_queries; i++) {
	  if (fgets(r[i], max_length, r_file) == NULL) { break; }
	  str_trim(r[i]);
	}
	read += sw_toc(start_read);
      }
      // end of #pragma omp critical

      // exit if no queries
      if (num_queries == 0) break;

      //printf("num queries = %i, count = %i\n", num_queries, count);
      start_memory = sw_tic();
      output_p = sw_multi_output_new(num_queries);
      memory += sw_toc(start_memory);

      // call smith-waterman
      partial_t = sw_tic();
      smith_waterman_mqmr(q, r, num_queries, optarg_p, 1, output_p);
      sse_t += sw_toc(partial_t);

      #pragma omp critical
      {
	// save results
	start_write = sw_tic();
	sw_multi_output_save(num_queries, output_p, out_file);
	write += sw_toc(start_write);
      } 
      // end of #pragma omp critical

      start_memory = sw_tic();
      sw_multi_output_free(output_p);
      memory += sw_toc(start_memory);

      batches++;
    }

    start_memory = sw_tic();
    // free memory and close files
    sw_optarg_free(optarg_p);

    for (int i = 0; i < batch_size; i++) {
      free(q[i]);
      free(r[i]);
    }
    memory += sw_toc(start_memory);

#ifdef TIMING
    //if (sse_matrix_t[thread_id] + sse_tracking_t[thread_id] > max_sse)
    //    max_sse = sse_matrix_t[thread_id] + sse_tracking_t[thread_id];
    //}
#endif

  }
  // end of #pragma omp parallel

  fclose(q_file);
  fclose(r_file);
  fclose(out_file);

#ifdef SW_AVX2
  printf("smith_waterman_mqmr function (AVX2)\n");
#else
  printf("smith_waterman_mqmr function (SSE)\n");
#endif

#ifdef TIMING
    printf("Thread: <matrix creation time> <backtracing time>\n");
    for(int i = 0; i < num_threads ; i++) {
      printf("\tThread %i:\t%0.3fs\t%0.3fs\n", i, sse_matrix_t[i], sse_tracking_t[i]);
    }
    printf("\tMax. time: %0.3fs\n", max_sse);
#endif
#ifdef TIMING
    //    printf("Memory mng. time    : %0.3f s\n", memory);
    //    printf("Read sequences time : %0.3f s\n", read);
    //    printf("Write results time  : %0.3f s\n", write);
#endif
    //printf("Alignment time      : %0.3f s (2 x %i seqs; %i threads)\n", sse_t, count, num_threads);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
