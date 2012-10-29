#include "batch_aligner.h"

//====================================================================================
//  main batch aligner function
//====================================================================================

void batch_aligner(batch_aligner_input_t *input) {

  //  printf("START: batch_aligner\n", omp_get_thread_num());

  size_t total_batches = 0;
		
  list_t *read_list = input->read_list;
  list_t *write_list = input->write_list;

  write_batch_t* write_batch = NULL;
  aligner_batch_t *aligner_batch = NULL;

  list_item_t *read_item = NULL, *write_item = NULL;

  unsigned int tid = omp_get_thread_num();

  struct timeval t1, t2;

  array_list_t *list1, *list2;

  // main loop
  while ( (read_item = list_remove_item(read_list)) != NULL ) {

    aligner_batch = aligner_batch_new((fastq_batch_t *) read_item->data_p);

    thr_batches[tid]++;

    //printf("********************** BATCH %d (batch aligner %d)\n", total_batches, omp_get_thread_num());

    // Burros-Wheeler transform
    gettimeofday(&t1, NULL);
    apply_bwt(input->bwt_input, aligner_batch);
    gettimeofday(&t2, NULL);
    bwt_time[tid] += ((t2.tv_sec - t1.tv_sec) * 1e6 + (t2.tv_usec - t1.tv_usec));
    //printf("---> %d, bwt, num targets = %d\n", tid, aligner_batch->num_targets);


    if (aligner_batch->num_targets > 0) {
      // seeding
      gettimeofday(&t1, NULL);
      apply_seeding(input->region_input, aligner_batch);
      gettimeofday(&t2, NULL);
      seeding_time[tid] += ((t2.tv_sec - t1.tv_sec) * 1e6 + (t2.tv_usec - t1.tv_usec));
      thr_seeding_items[tid] += aligner_batch->num_targets;
      //printf("---> %d, seeding, num targets = %d\n", tid, aligner_batch->num_targets);

      // seeking CALs
      gettimeofday(&t1, NULL);
      apply_caling(input->cal_input, aligner_batch);
      gettimeofday(&t2, NULL);
      cal_time[tid] += ((t2.tv_sec - t1.tv_sec) * 1e6 + (t2.tv_usec - t1.tv_usec));
      thr_cal_items[tid] += aligner_batch->num_targets;
      //printf("---> %d, cal, num targets = %d\n", tid, aligner_batch->num_targets);
    }

    // pair-mode managing
    if (input->pair_input != NULL) {      
      apply_pair(input->pair_input, aligner_batch);
      //      printf("---> %d, pair, num targets = %d\n", tid, aligner_batch->num_targets);
    }

    if (aligner_batch->num_targets > 0) {
      // Smith-Waterman
      gettimeofday(&t1, NULL);
      apply_sw(input->sw_input, aligner_batch);
      gettimeofday(&t2, NULL);
      sw_time[tid] += ((t2.tv_sec - t1.tv_sec) * 1e6 + (t2.tv_usec - t1.tv_usec));
      thr_sw_items[tid] += aligner_batch->num_targets;
      //printf("---> %d, sw, num targets = %d\n", tid, aligner_batch->num_targets);
    }

    if (aligner_batch->num_targets > 0) {
      // prepare alignments (converts sw-output to alignment, searches pairs...)
      prepare_alignments(input->pair_input, aligner_batch);
    }

    write_item = list_item_new(total_batches, 0, aligner_batch);
    list_insert_item(write_item, write_list);

    list_item_free(read_item);
    total_batches++;
  } // main loop

  /*
  printf("Thread %d: BWT time     = %0.4f s\n", tid, bwt_time / 1e6);
  printf("Thread %d: Seeding time = %0.4f s\n", tid, seeding_time / 1e6);
  printf("Thread %d: CAL time     = %0.4f s\n", tid, cal_time / 1e6);
  printf("Thread %d: SW time      = %0.4f s\n", tid, sw_time / 1e6);
  */

  // decreasing writers
  if (write_list != NULL) list_decr_writers(write_list);

  //  printf("END: batch_aligner (%d), (total batches %d): END\n", omp_get_thread_num(), total_batches);
}

//====================================================================================
//  structures and prototypes
//====================================================================================

void batch_aligner_input_init(list_t *read_list, list_t *write_list, 
			      bwt_server_input_t *bwt_input, 
			      region_seeker_input_t *region_input,
			      cal_seeker_input_t *cal_input,
			      pair_server_input_t *pair_input,
			      sw_server_input_t *sw_input,
			      batch_aligner_input_t *input) {
  input->read_list = read_list;
  input->write_list = write_list;

  input->bwt_input = bwt_input;
  input->region_input = region_input;
  input->cal_input = cal_input;
  input->pair_input = pair_input;
  input->sw_input = sw_input;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
