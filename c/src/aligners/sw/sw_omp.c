
#include "macros.h"
#include "sw_commons.h"
#include "smith_waterman.h"

#include "sw_omp.h"

#ifdef TIMING
extern double *sse_matrix_t, *sse_tracking_t;
#endif // TIMING

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

  double start = 0.0f, elapsed = 0.0f;
  start = sw_tic();

#ifdef SW_AVX2
  const unsigned int depth = 8;
#else
  const unsigned int depth = 4;
#endif

  const int max_length = 2048;

  FILE *q_file = fopen(q_filename, "r");
  FILE *r_file = fopen(r_filename, "r");

  FILE *out_file = NULL;
  if (out_filename) {
    out_file = fopen(out_filename, "w");
  }

#ifdef TIMING
  int batches[num_threads];
  double start_read[num_threads], read[num_threads];
  double start_write[num_threads], write[num_threads];
  double start_memory[num_threads], memory[num_threads];
  double start_convert[num_threads], convert[num_threads];
  double start_sw_back[num_threads], sw_back[num_threads];
#endif

  #pragma omp parallel num_threads(num_threads)
  {
    sw_optarg_t *optarg_p;
    sw_multi_output_t *output_p;
    int tid = omp_get_thread_num();

    int buffer_size;
    char *buffer = (char *) malloc(batch_size * 4096);

#ifdef TIMING
    batches[tid] = 0;
    read[tid] = 0;
    convert[tid] = 0;
    write[tid] = 0;
    memory[tid] = 0;
    sw_back[tid] = 0;
#endif // TIMING

    #pragma omp critical
    {
      optarg_p = sw_optarg_new(match, mismatch, gap_open, gap_extend, matrix_filename);
    }

#ifdef TIMING
    start_memory[tid] = sw_tic();
#endif
    char *q[batch_size], *r[batch_size];
    for (int i = 0; i < batch_size; i++) {
      q[i] = (char *) calloc(max_length, sizeof(char));
      r[i] = (char *) calloc(max_length, sizeof(char));
    }
#ifdef TIMING
    memory[tid] += sw_toc(start_memory[tid]);
#endif

    int count = 0, num_queries;

    while (1) {
      num_queries = 0;

      #pragma omp critical
      {
#ifdef TIMING
	start_read[tid] = sw_tic();
#endif
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
#ifdef TIMING
	read[tid] += sw_toc(start_read[tid]);
#endif
      }
      // end of #pragma omp critical

      // exit if no queries
      if (num_queries == 0) break;

      //printf("num queries = %i, count = %i\n", num_queries, count);
#ifdef TIMING
      start_memory[tid] = sw_tic();
#endif
      output_p = sw_multi_output_new(num_queries);
#ifdef TIMING
      memory[tid] += sw_toc(start_memory[tid]);
#endif
      
      // call smith-waterman
#ifdef TIMING
      start_sw_back[tid] = sw_tic();
#endif
      smith_waterman_mqmr(q, r, num_queries, optarg_p, 1, output_p);
#ifdef TIMING
      sw_back[tid] += sw_toc(start_sw_back[tid]);
#endif

      if (out_file) {
#ifdef TIMING
        start_convert[tid] = sw_tic();
#endif
        buffer[0] = '\0';
        buffer_size = sw_multi_output_string(num_queries, output_p, buffer);
#ifdef TIMING
        convert[tid] += sw_toc(start_convert[tid]);
#endif

        #pragma omp critical
	{
	  // save results
#ifdef TIMING
	  start_write[tid] = sw_tic();
#endif
	  fwrite(buffer, 1, buffer_size, out_file);
	  //	  sw_multi_output_save(num_queries, output_p, out_file);
#ifdef TIMING
	  write[tid] += sw_toc(start_write[tid]);
#endif
        }
        // end of #pragma omp critical
      }

#ifdef TIMING
      start_memory[tid] = sw_tic();
#endif
      sw_multi_output_free(output_p);
#ifdef TIMING
      memory[tid] += sw_toc(start_memory[tid]);
#endif

#ifdef TIMING
      batches[tid]++;
#endif // TIMING
  }

#ifdef TIMING
    start_memory[tid] = sw_tic();
#endif
    // free memory
    sw_optarg_free(optarg_p);
    for (int i = 0; i < batch_size; i++) {
      free(q[i]);
      free(r[i]);
    }
#ifdef TIMING
    memory[tid] += sw_toc(start_memory[tid]);
#endif
  }
  // end of #pragma omp parallel

  elapsed += sw_toc(start);

  fclose(q_file);
  fclose(r_file);
  if (out_file) {
    fclose(out_file);
  }

#if defined (SW_AVX512) || defined (__AVX512__)
    printf("AVX512 (time in seconds)\n");
#else
#if defined (SW_AVX2) || defined (__AVX2__)
  printf("AVX2 (time in seconds)\n");
#else
  printf("SSE (time in seconds)\n");
#endif
#endif

#ifdef TIMING
  double batches_total = 0, read_total = 0, write_total = 0;
  double max_batches = 0.0f, max_sw_back = 0.0f, max_convert = 0.0f, max_memory = 0.0f;
  double max_matrix = 0.0f, max_back = 0.0f;
  double total_sw_back = 0.0f, total_matrix = 0.0f, total_back = 0.0f, total_memory = 0.0f, total_convert = 0.0f;
  printf("Thread: <batches> <matrix + backtracing> <matrix> <backtracking> <read> <convert> <write> <memory>\n");
  for(int i = 0; i < num_threads ; i++) {
    printf("Thread_%i\t%i\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", i, batches[i], sw_back[i], sse_matrix_t[i], sse_tracking_t[i], read[i], convert[i], write[i], memory[i]);
    batches_total += batches[i];
    read_total += read[i];
    write_total += write[i];

    total_matrix += sse_matrix_t[i];
    total_back += sse_tracking_t[i];
    total_sw_back += sw_back[i];
    total_memory += memory[i];
    total_convert += convert[i];

    if (batches[i] > max_batches) max_batches = batches[i];
    if (sw_back[i] > max_sw_back) max_sw_back = sw_back[i];
    if (sse_matrix_t[i] > max_matrix) max_matrix = sse_matrix_t[i];
    if (sse_tracking_t[i] > max_back) max_back = sse_tracking_t[i];
    if (sw_back[i] > max_sw_back) max_sw_back = sw_back[i];
    if (convert[i] > max_convert) max_convert = convert[i];
    if (memory[i] > max_memory) max_memory = memory[i];
  }
  printf("Batches (max, mean)\t%0.1f\t%0.3f\n", max_batches, batches_total / num_threads);
  printf("Total reading time\t%0.3f\n", read_total);
  printf("Total writing time\t%0.3f\n", write_total);
  printf("Matrix + backtracking (max, mean)\t%0.3f\t%0.3f\n", max_sw_back, total_sw_back / num_threads);
  printf("\tMatrix (max, mean)\t%0.3f\t%0.3f\n", max_matrix, total_matrix / num_threads);
  printf("\tBacktracking (max, mean)\t%0.3f\t%0.3f\n", max_back, total_back / num_threads);
  printf("Memory management time (max, mean)\t%0.3f\t%0.3f\n", max_memory, total_memory / num_threads);
  printf("Conversion time (max, mean)\t%0.3f\t%0.3f\n", max_convert, total_convert / num_threads);

  //    printf("Alignment time      : %0.3f s (%i threads)\n", sw_back, num_threads);
#endif
    printf("\n");
    printf("Total hpg-sw time   : %0.3f\n", elapsed);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
