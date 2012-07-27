#include "main.h"

double emboss_matrix_t = 0.0f, emboss_tracking_t = 0.0f;
double sse_matrix_t = 0.0f, sse_tracking_t = 0.0f;
double sse1_matrix_t = 0.0f, sse1_tracking_t = 0.0f;
double avx_matrix_t = 0.0f, avx_tracking_t = 0.0f;
double avx1_matrix_t = 0.0f, avx1_tracking_t = 0.0f;

double sse_t = 0.0f;

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
void test();

void main(int argc, char *argv[]) {

  float match = 5.0f;
  float mismatch = 5.0f;
  float gap_open = 10.0f;
  float gap_extend = 0.5f;
  char *q_filename, *r_filename;
  char *emboss_out_filename, *sse_out_filename, *sse1_out_filename, *avx_out_filename, *avx1_out_filename;
  int num_threads = 4;
  int batch_size = 2048;

  if (argc != 7) {
    printf("Error.\n");
    printf("Usage: %s queries-filename refs-filename emboss-out-filename sse-out-filename num-threads batch-size\n", argv[0]);
    exit(-1);
  }
  q_filename = argv[1]; //"../data/queries-50k.txt";
  r_filename = argv[2]; //"../data/refs-50k.txt";
  emboss_out_filename = argv[3]; // "./emboss.out";
  sse_out_filename = argv[4]; //"./sse.out";
  num_threads = atoi(argv[5]);
  batch_size = atoi(argv[6]);
  //  sse1_out_filename = argv[5]; //"./sse1.out";
  //  avx_out_filename = argv[6]; //"./avx.out";
  //  avx1_out_filename = argv[7]; //"./avx1.out";

  // EMBOSS
  run_emboss(q_filename, r_filename,
	     gap_open, gap_extend, match, mismatch,
	     emboss_out_filename);

  //emboss_matrix_t = 54.7; // 37.93; //29.45;;
  //  emboss_tracking_t = 7.42; // 5.188; // 4.35;
  printf("EMBOSS\n");
  printf("\tCalculating matrix: %0.5f s\n", emboss_matrix_t);
  printf("\tTracking back     : %0.5f s\n", emboss_tracking_t);
  printf("\tTotal             : %0.5f s\n", (emboss_matrix_t + emboss_tracking_t));

  // SSE
  run_sse(q_filename, r_filename,
	  gap_open, gap_extend, "./dnafull", batch_size, num_threads, sse_out_filename);

  
  printf("SSE\n");
  printf("\tCalculating matrix: %0.5f s\n", sse1_matrix_t);
  printf("\tTracking back     : %0.5f s\n", sse1_tracking_t);
  printf("\tTotal             : %0.5f s\n", (sse1_matrix_t + sse1_tracking_t));
  printf("\tTotal (function)  : %0.2f s\n", sse_t);
  
  printf("\nSpeed-up (EMBOS vs SSE)\n");
  printf("\tCalculating matrix: %0.2f\n", emboss_matrix_t / sse1_matrix_t);
  printf("\tTracking back     : %0.2f\n", emboss_tracking_t / sse1_tracking_t);
  printf("\tTotal             : %0.2f\n", (emboss_matrix_t + emboss_tracking_t) / (sse1_matrix_t + sse1_tracking_t));
  printf("\tTotal (function)  : %0.2f\n", (emboss_matrix_t + emboss_tracking_t) / (sse_t));
  printf("\n");
}

//-------------------------------------------------------------------------

void run_emboss(char *q_filename, char *r_filename, 
		  float gap_open, float gap_extend,
		  float match, float mismatch,
		  char *out_filename) {  

  const int max_length = 2048;
  double total_t = 0.0f, partial_t = 0.0f;

  int count = 0, start1, start2;

  float score;
  char *m, *n;

  FILE *q_file = fopen(q_filename, "r");
  FILE *r_file = fopen(r_filename, "r");
  FILE *out_file = fopen(out_filename, "w");

  char q[max_length + 100], r[max_length + 100];

  while (fgets(q, max_length, q_file) != NULL) {
    if (fgets(r, max_length, r_file) == NULL) {
      break;
    }

    trim(q);
    trim(r);

    m = (char *) calloc(4048, sizeof(char));
    n = (char *) calloc(4048, sizeof(char));
    
    //    partial_t = tic();
    score = sw_emboss(r, q, gap_open, gap_extend, m, n, &start1, &start2);
    //    emboss_t += toc(partial_t);

    fprintf(out_file, "%s\n%s\nScore = %0.2f, query start at %i, ref. start at %i\n\n", n, m, score,
	      start2, start1);

    free(m);
    free(n);

    count++;
    //if (count >= 1) break;
  }
  printf("\nEMBOSS done (%i reads)\n", count);

  fclose(q_file);
  fclose(r_file);
  fclose(out_file);

  //  return total_t;
}

//-------------------------------------------------------------------------

void run_sse(char *q_filename, char *r_filename, 
	     float gap_open, float gap_extend, char *matrix_filename,
	     int batch_size, int num_threads,
	     char *out_filename) {  

  const int max_length = 2048;
  double total_t = 0.0f, partial_t = 0.0f;


  sw_optarg_t *optarg_p = sw_optarg_new(gap_open, gap_extend, matrix_filename);
  sw_multi_output_t *output_p;


  FILE *q_file = fopen(q_filename, "r");
  FILE *r_file = fopen(r_filename, "r");
  FILE *out_file = fopen(out_filename, "w");

  char *q[batch_size], *r[batch_size];
  for (int i = 0; i < batch_size; i++) {
    q[i] = (char *) calloc(max_length, sizeof(char));
    r[i] = (char *) calloc(max_length, sizeof(char));
  }
  
  int count = 0, batches = 0, num_queries;

  while (1) {
    num_queries = 0;

    // read queries
    for (int i = 0; i < batch_size; i++) {
      if (fgets(q[i], max_length, q_file) == NULL) { break; }
      trim(q[i]);
      num_queries++;
      count++;
    }

    // exit if no queries
    if (num_queries == 0) break;

    // read references
    for (int i = 0; i < num_queries; i++) {
      if (fgets(r[i], max_length, r_file) == NULL) { break; }
      trim(r[i]);
    }

    //printf("num queries = %i, count = %i\n", num_queries, count);
    output_p = sw_multi_output_new(num_queries);

    // call smith-waterman
    partial_t = tic();
    smith_waterman_mqmr(q, r, num_queries, optarg_p, num_threads, output_p);
    sse_t += toc(partial_t);

    // save results
    sw_multi_output_save(num_queries, output_p, out_file);
    //break;

    sw_multi_output_free(output_p);

    batches++;
    //printf("%i batches\n", batches);
    //    if (batches == 7) break;
  }
  printf("\nSSE (new version) done (%i reads)\n", count);

  // free memory and close files
  sw_optarg_free(optarg_p);

  for (int i = 0; i < batch_size; i++) {
    free(q[i]);
    free(r[i]);
  }

  fclose(q_file);
  fclose(r_file);
  fclose(out_file);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
