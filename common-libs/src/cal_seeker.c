#include "cal_seeker.h"

//------------------------------------------------------------------------------------
// cal_seeker_input functions: init
//------------------------------------------------------------------------------------
void cal_seeker_server(cal_seeker_input_t* input) {
  
  extern short int cal_type;
  
  unsigned int cal_id = omp_get_thread_num();

  LOG_DEBUG_F("cal_seeker_server(%d): START\n", cal_id); 
  
  list_item_t *item = NULL;
  list_item_t *write_item = NULL;
  list_item_t *sw_item = NULL;
  size_t num_reads;
  array_list_t *allocate_cals;
  mapping_batch_t *mapping_batch;
  size_t num_cals, select_cals, total_cals = 0;
  fastq_read_t *read;

  char *seq;
  char *header;
  char *quality;
  unsigned int seq_len, header_len;

  char *cigar;
  char *header_id;
  unsigned int write_size = input->batch_size;
  //write_batch_t* write_batch = write_batch_new(write_size, MATCH_FLAG);
  
  //list_t *write_list = input->write_list;
  sw_batch_t *sw_batch;
  fastq_read_t *allocate_reads; 
  size_t reads_with_cals;
  alignment_t *alignment;
  size_t num_batches = 0, num_reads_unmapped = 0, num_without_cals = 0;
  size_t total_reads = 0;
  size_t j;
  size_t min_seeds, max_seeds, num_targets, target_pos, total_targets;

  while ( (item = list_remove_item(input->regions_list)) != NULL ) {
    //printf("Cal Seeker %d processing batch...\n", omp_get_thread_num());
    num_batches++;
    if (time_on) { timing_start(CAL_SEEKER, cal_id, timing_p); }
    
    mapping_batch = (mapping_batch_t *)item->data_p;
    num_targets = mapping_batch->num_targets;
    target_pos = 0;
    total_targets = 0;
    //num_reads = cal_batch->unmapped_batch->num_reads;
    total_reads += num_targets;
    //allocate_cals = (array_list_t **)calloc(num_reads, sizeof(array_list_t *));
    //allocate_reads = (fastq_read_t **)calloc(num_reads, sizeof(fastq_read_t *));
    //reads_with_cals = 0;
    for (size_t i = 0; i < num_targets; i++) {
      allocate_cals = array_list_new(100, 
				     1.25f, 
				     COLLECTION_MODE_ASYNCHRONIZED);

      /*                    
      num_cals = bwt_generate_cal_list_rna_linkedlist(mapping_batch->mapping_lists[mapping_batch->targets[i]], 
						      input->cal_optarg, 
						      allocate_cals);  
      */
      
      num_cals = bwt_generate_cal_list_rna_linked_list(mapping_batch->mapping_lists[mapping_batch->targets[i]], 
						      input->cal_optarg, 
						      allocate_cals);  
      
      array_list_free(mapping_batch->mapping_lists[mapping_batch->targets[i]], region_bwt_free);
      mapping_batch->mapping_lists[mapping_batch->targets[i]] = allocate_cals;
      //printf("Num CALS %i\n", num_cals);
      if (num_cals > MAX_RNA_CALS) {
	//array_list_clear(allocate_cals, (void *)cal_free);
	select_cals = num_cals - MAX_RNA_CALS;

	for(j = num_cals - 1; j >= MAX_RNA_CALS; j--) {
	  cal_free(array_list_remove_at(j, mapping_batch->mapping_lists[mapping_batch->targets[i]]));
	}
	mapping_batch->targets[target_pos++] = mapping_batch->targets[i];
	num_reads_unmapped++;
      }else if (num_cals > 0) {
	mapping_batch->targets[target_pos++] = mapping_batch->targets[i];
      }else if (!num_cals) {
	num_reads_unmapped++;
      }
    
    }
    
    mapping_batch->num_targets = target_pos;

    //sw_batch = sw_batch_new(reads_with_cals, allocate_cals, allocate_reads);
    //sw_batch_free(sw_batch_p);
    //sw_item = list_item_new(0, 0, sw_batch);
    
    if (time_on) { timing_stop(CAL_SEEKER, cal_id, timing_p); }      
    list_insert_item(item, input->sw_list); 
    //printf("Cal Seeker process batch finish!\n");
  }

  list_decr_writers(input->sw_list);
  LOG_DEBUG_F("cal_seeker_server(%lu reads unmapped): END\n", num_reads_unmapped); 
  // free memory for mapping list
   
}

//====================================================================================
// apply_caling
//====================================================================================
//int unmapped_by_max_cals_counter[100];
//int unmapped_by_zero_cals_counter[100];

void apply_caling(cal_seeker_input_t* input, mapping_batch_t *batch) {

  //  printf("START: apply_caling\n"); 
  int tid = omp_get_thread_num();

  array_list_t *list = NULL;
  size_t read_index, num_cals, min_seeds, max_seeds;
  int min_limit;

  cal_t *cal;
  array_list_t *cal_list;

  size_t num_targets = batch->num_targets;
  size_t *targets = batch->targets;
  size_t new_num_targets = 0;
  //  size_t *new_targets = (size_t *) calloc(num_targets, sizeof(size_t));
  
  // set to zero
  batch->num_to_do = 0;

  for (size_t i = 0; i < num_targets; i++) {

    read_index = targets[i];
    
    if (!list) {
      list = array_list_new(1000, 
			    1.25f, 
			    COLLECTION_MODE_ASYNCHRONIZED);
    }

    //    printf("cal_seeker.c: %s\n", &(batch->fq_batch->header[batch->fq_batch->header_indices[index]]));
    //    printf("\tcal_seeker.c: array_list_size = %d\n", array_list_size(list));
    /*
    {
      // for debugging
      region_t *region;
      int size = array_list_size(batch->mapping_lists[index]);
      for (int i = 0; i < size; i++) {
	region = array_list_get(i, batch->mapping_lists[index]);
	printf("region %i: strand %d chromosome %d [%d-%d], seq [%d-%d]\n", 
	       i, region->strand, region->chromosome_id, region->start, region->end,
	       region->seq_start, region->seq_end);
      }
    }
    */
    // optimized version
    num_cals = bwt_generate_cal_list_linkedlist(batch->mapping_lists[read_index], 
						input->cal_optarg,
						&min_seeds, &max_seeds,
						list);

    /*
    {
      // for debugging
      cal_t *cal;
      int size = array_list_size(list);
      for (int i = 0; i < size; i++) {
	cal = array_list_get(i, list);
	printf("\tcal %i: strand %d chromosome %d [%d-%d]\n", 
	       i, cal->strand, cal->chromosome_id, cal->start, cal->end);
      }
    }
    */

    //    printf("\tcal_seeker.c: num_cals = %d, array_list_size = %d, (MAX = %d), num seeds (min, max) = (%d, %d)\n", 
    //	   num_cals, array_list_size(list), MAX_CALS, min_seeds, max_seeds);

    // filter CALs by the number of seeds
    if (min_seeds == max_seeds) {
      cal_list = list;
      list = NULL;
    } else {
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (size_t j = 0; j < num_cals; j++) {
	cal = array_list_get(j, list);
	if (cal->num_seeds == max_seeds) {
	  array_list_insert(cal, cal_list);
	  array_list_set(j, NULL, list);
	}
      }
      array_list_clear(list, (void *) cal_free);
      num_cals = array_list_size(cal_list);
    }


    /*
    min_limit = max_seeds - 3;
    if (min_seeds == max_seeds || min_limit < min_seeds) {
      cal_list = list;
      list = NULL;
    } else {
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (size_t j = 0; j < num_cals; j++) {
	cal = array_list_get(j, list);
	if (cal->num_seeds > min_limit) {
	  array_list_insert(cal, cal_list);
	  array_list_set(j, NULL, list);
	}
      }
      array_list_clear(list, (void *)cal_free);
      num_cals = array_list_size(cal_list);
    }
    */
    
    if (num_cals > MAX_CALS) {
      for (size_t j = num_cals - 1; j >= MAX_CALS; j--) {
	cal = (cal_t *) array_list_remove_at(j, cal_list);
	cal_free(cal);
      }
      num_cals = array_list_size(cal_list);
    }

    //    printf("\tcal_seeker.c: after filter: num_cals = %d\n", num_cals);

    if (num_cals > 0 && num_cals <= MAX_CALS) {
      array_list_set_flag(2, cal_list);
      batch->num_to_do += num_cals;

      targets[new_num_targets++] = read_index;
      
      // we have to free the region list
      array_list_free(batch->mapping_lists[read_index], (void *) region_bwt_free);
      batch->mapping_lists[read_index] = cal_list;
    } else {
      array_list_set_flag(0, batch->mapping_lists[read_index]);
      /*
      if (num_cals > 0) {
	if (strncmp("@rand", &(batch->fq_batch->header[batch->fq_batch->header_indices[index]]), 5)) {
	  unmapped_by_max_cals_counter[tid]++;
	  //	  printf("--> %s\n", &(batch->fq_batch->header[batch->fq_batch->header_indices[index]]));
	}
      } else {
	if (strncmp("@rand", &(batch->fq_batch->header[batch->fq_batch->header_indices[index]]), 5)) {
	  unmapped_by_zero_cals_counter[tid]++;
	}
      }
      */
      // we have to free the region list
      array_list_clear(batch->mapping_lists[read_index], (void *) region_bwt_free);
      if (cal_list) array_list_free(cal_list, (void *) cal_free);
      if (list) array_list_clear(list, (void *) cal_free);
    }
  } // end for 0 ... num_seqs

  // update batch
  batch->num_targets = new_num_targets;
  //  batch->num_allocated_targets = num_targets;
  //  if (batch->targets) free(batch->targets);
  //  batch->targets = new_targets;
  
  // update counter
  //  thr_cal_items[tid] += batch->num_done;

  // free memory
  if (list) array_list_free(list, NULL);
  /*
  {
    // displaying for debugging
    cal_t *cal;
    for (size_t i = 0; i< num_outputs; i++) {
      printf("\tlist %d (read %d), size = %d\n", 
	     i, outputs[i], array_list_size(batch->mapping_lists[outputs[i]]));
      for (size_t j = 0; j < array_list_size(batch->mapping_lists[outputs[i]]) ; j++) {
	cal = array_list_get(j, batch->mapping_lists[outputs[i]]);
	printf("\t\tcal %d: strand = %d, start = %d, end = %d\n", j, cal->strand, cal->start, cal->end);
      }
    }
  }
  */  
  //  printf("END: apply_caling, (caling  %d reads)\n", num_outputs);
}

//------------------------------------------------------------------------------------

void cal_seeker_input_init(list_t *regions_list, cal_optarg_t *cal_optarg, 
			   list_t* write_list, unsigned int write_size, 
			   list_t *sw_list, list_t *pair_list, 
			   cal_seeker_input_t *input){
  input->regions_list = regions_list;
  input->cal_optarg = cal_optarg;
  input->batch_size = write_size;
  input->sw_list = sw_list;
  input->pair_list = pair_list;
  input->write_list = write_list;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
