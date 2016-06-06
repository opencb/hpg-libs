
#include "commons/workflow_scheduler.h"
#include "containers/array_list.h"

#include "macros.h"
#include "sw_commons.h"
#include "smith_waterman.h"

#include "sw_pipeline.h"

//====================================================================
// W O R K F L O W     F O R      H P G - S W
//====================================================================

#define CONSUMER_STAGE   -1

//--------------------------------------------------------------------
// workflow input
//--------------------------------------------------------------------

typedef struct sw_wf_input {
  FILE *query_file;
  FILE *ref_file;
  FILE *out_file;
  sw_optarg_t *options;
  int batch_size;
} sw_wf_input_t;

sw_wf_input_t *sw_wf_input_new(FILE *query_file,
			       FILE *ref_file,
			       FILE *out_file,
			       sw_optarg_t *opts,
			       int batch_size) {
  
  sw_wf_input_t *b = (sw_wf_input_t *) calloc(1, sizeof(sw_wf_input_t));

  b->query_file = query_file;
  b->ref_file = ref_file;
  b->out_file = out_file;
  b->options = opts;
  b->batch_size = batch_size;

  return b;
}

void sw_wf_input_free(sw_wf_input_t *b) {
  if (b) {
    free(b);
  }
}

//--------------------------------------------------------------------
// batch structure between the different workflow stages
//--------------------------------------------------------------------

typedef struct sw_wf_batch {
  int num_queries;
  int batch_size;
  sw_optarg_t *options;
  sw_multi_output_t *output;
  char **q;
  char **r;
  FILE *out_file;
} sw_wf_batch_t;

sw_wf_batch_t *sw_wf_batch_new(int num_queries, int batch_size, 
			       sw_optarg_t *opts,
			       char **q, char **r, FILE *out_file) {
  sw_wf_batch_t *b = (sw_wf_batch_t *) calloc(1, sizeof(sw_wf_batch_t));
  
  b->num_queries = num_queries;
  b->batch_size = batch_size;
  b->options = opts;
  b->output = NULL;
  b->q = q;
  b->r = r;
  b->out_file = out_file;
  
  return b;
}

void sw_wf_batch_free(sw_wf_batch_t *b) {
  if (b) {
    int size = b->batch_size;
    for (int i = 0; i < size; i++) {
        free(b->q[i]);
        free(b->r[i]);
    }
    free(b->q);
    free(b->r);

    sw_multi_output_free(b->output);

    free(b);
  }
}

//--------------------------------------------------------------------
// workflow producer : read sequences from query and ref. files
//--------------------------------------------------------------------

void *sw_producer(void *input) { 
  sw_wf_input_t *wf_input = (sw_wf_input_t *) input;
  sw_wf_batch_t *wf_batch = NULL;
  int batch_size = wf_input->batch_size;

  const int max_length = 2048;

  char **q = (char **) malloc(batch_size * sizeof(char *));
  char **r = (char **) malloc(batch_size * sizeof(char *));
  for (int i = 0; i < batch_size; i++) {
    q[i] = (char *) calloc(max_length, sizeof(char));
    r[i] = (char *) calloc(max_length, sizeof(char));
  }

  // read queries
  int num_queries = 0;
  for (int i = 0; i < batch_size; i++) {
    if (fgets(q[i], max_length, wf_input->query_file) == NULL) { break; }
    str_trim(q[i]);
    num_queries++;
  }

  // read references
  for (int i = 0; i < num_queries; i++) {
    if (fgets(r[i], max_length, wf_input->ref_file) == NULL) { break; }
    str_trim(r[i]);
  }
  
  if (num_queries > 0) {
    wf_batch = sw_wf_batch_new(num_queries, batch_size, wf_input->options,
			       q, r, wf_input->out_file);      
  } else {
    // free memory
    for (int i = 0; i < batch_size; i++) {
      free(q[i]);
      free(r[i]);
    }
    free(q);
    free(r);
  }
    
  return wf_batch;
}

//--------------------------------------------------------------------
// workflow worker : align sequences
//--------------------------------------------------------------------

int sw_worker(void *data) {
  sw_wf_batch_t *batch = (sw_wf_batch_t *) data;

  batch->output = sw_multi_output_new(batch->num_queries);

  smith_waterman_mqmr(batch->q, batch->r, 
		      batch->num_queries, batch->options, 1, batch->output);
    
  return CONSUMER_STAGE;
}

//--------------------------------------------------------------------
// workflow consumer : write alignments
//--------------------------------------------------------------------

int sw_consumer(void *data) {
  sw_wf_batch_t *batch = (sw_wf_batch_t *) data;

  if (batch->num_queries > 0) {
    sw_multi_output_save(batch->num_queries, batch->output, batch->out_file);
  }

  // free memory
  sw_wf_batch_free(batch);
}

//--------------------------------------------------------------------
// workflow description
//--------------------------------------------------------------------

void run_sw_pipeline(char *q_filename, char *r_filename,
		     float match, float mismatch, float gap_open, float gap_extend,
		     char *matrix_filename, int batch_size, int num_threads,
		     char *out_filename) {

    sw_optarg_t *options = sw_optarg_new(match, mismatch, gap_open, gap_extend, matrix_filename);

    FILE *q_file = fopen(q_filename, "r");
    FILE *r_file = fopen(r_filename, "r");
    FILE *out_file = fopen(out_filename, "w");

    //------------------------------------------------------------------
    // workflow management
    //
    sw_wf_input_t *input = sw_wf_input_new(q_file, r_file, out_file, options, batch_size);

    // create and initialize workflow
    workflow_t *wf = workflow_new();
  
    workflow_stage_function_t stage_functions[] = {sw_worker};
    char *stage_labels[] = {"SW worker"};
    workflow_set_stages(1, &stage_functions, stage_labels, wf);
  
    // optional producer and consumer functions
    workflow_set_producer(sw_producer, "SW producer", wf);
    workflow_set_consumer(sw_consumer, "SW consumer", wf);
  
    workflow_run_with(num_threads, input, wf);

#ifdef TIMING
    workflow_display_timing(wf);
#endif

    //
    // end of workflow management
    //------------------------------------------------------------------

    // free memory and close files
    workflow_free(wf);
    sw_wf_input_free(input);
    sw_optarg_free(options);
    fclose(q_file);
    fclose(r_file);
    fclose(out_file);
}

//-------------------------------------------------------------------------
/*
    for(int i = 0; i < num_threads ; i++) {
      printf("\tThread %i:\t%0.3fs\t%0.3fs\n", i, sse_matrix_t[i], sse_tracking_t[i]);
    }
    printf("\tMax. time: %0.3fs\n", max_sse);
#endif
#ifdef TIMING
    printf("Memory mng. time    : %0.3f s\n", memory);
    printf("Read sequences time : %0.3f s\n", read);
    printf("Write results time  : %0.3f s\n", write);
#endif
    printf("Alignment time      : %0.3f s (2 x %i seqs; %i threads)\n", sse_t, count, num_threads);
    printf("\n");
    printf("Total hpg-sw time   : %0.3fs\n", elapsed);
}
*/
