
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
  double max_sw_back = 0.0f;
  double start_read[num_threads], read[num_threads];
  double start_write[num_threads], write[num_threads];
  double start_memory[num_threads], memory[num_threads];
  double start_format[num_threads], format[num_threads];
  double start_sw_back[num_threads], sw_back[num_threads];
#endif

  #pragma omp parallel num_threads(num_threads)
  {
    sw_optarg_t *optarg_p;
    sw_multi_output_t *output_p;
    int tid = omp_get_thread_num();

    int buffer_size;
    char *buffer = (char *) malloc(batch_size * 4096);

    batches[tid] = 0;
    read[tid] = 0;
    format[tid] = 0;
    write[tid] = 0;
    memory[tid] = 0;
    sw_back[tid] = 0;

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
        start_format[tid] = sw_tic();
#endif
        buffer[0] = '\0';
        buffer_size = sw_multi_output_string(num_queries, output_p, buffer);
#ifdef TIMING
        format[tid] += sw_toc(start_format[tid]);
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

      batches[tid]++;
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

#ifdef TIMING
    for(int i = 0 ; i < num_threads; i++) {
      if (sw_back[i] > max_sw_back)
        max_sw_back = sw_back[i];
    }
#endif

  fclose(q_file);
  fclose(r_file);
  if (out_file) {
    fclose(out_file);
  }

#ifdef SW_AVX2
  printf("AVX2 (time in seconds)\n");
#else
  printf("SSE (time in seconds)\n");
#endif

#ifdef TIMING
  double batches_total = 0, read_total = 0, write_total = 0;
  printf("Thread: <batches> <matrix + backtracing> <matrix> <backtracking> <read> <format> <write> <memory>\n");
  for(int i = 0; i < num_threads ; i++) {
    printf("\tThread_%i\t%i\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", i, batches[i], sw_back[i], sse_matrix_t[i], sse_tracking_t[i], read[i], format[i], write[i], memory[i]);
    batches_total += batches[i];
    read_total += read[i];
    write_total += write[i];
  }
  printf("Max. SW + backtracking : %0.3f\n", max_sw_back);
  printf("Total reading time     : %0.3f\n", read_total);
  printf("Total writing time     : %0.3f\n", write_total);

  //    printf("Alignment time      : %0.3f s (%i threads)\n", sw_back, num_threads);
#endif
    printf("\n");
    printf("Total hpg-sw time   : %0.3f\n", elapsed);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
