#include "region_seeker.h"

void region_seeker_input_init(list_t *unmppaed_reads_list, cal_optarg_t *cal_optarg, 
			      bwt_optarg_t *bwt_optarg, bwt_index_t *bwt_index, 
			      list_t* region_list, region_seeker_input_t *input){
  input->unmppaed_reads_list_p = unmppaed_reads_list;
  input->cal_optarg_p = cal_optarg;
  input->bwt_optarg_p = bwt_optarg;
  input->bwt_index_p = bwt_index;
  input->region_list_p = region_list;
}

//--------------------------------------------------------------------------------------

void region_seeker_server(region_seeker_input_t *input){
  
  printf("region_seeker_server(%d): START\n", omp_get_thread_num());  
  list_item_t *item = NULL;
  list_item_t *region_item = NULL;
  fastq_batch_t *unmapped_batch;
  size_t num_reads;
  array_list_t **allocate_mapping;
  region_batch_t *region_batch;
  size_t num_mappings, total_mappings = 0, num_batches = 0;
  size_t num_threads = input->bwt_optarg_p->num_threads;
  size_t chunk = MAX(1, num_reads / (num_threads * 10));
  omp_set_num_threads(num_threads);

  while ( (item = list_remove_item(input->unmppaed_reads_list_p)) != NULL ) {
    //printf("Region Seeker Processing batch...\n");
    num_batches++;
    if (time_on) { timing_start(REGION_SEEKER, 0, timing_p); }
    
    unmapped_batch = (fastq_batch_t *)item->data_p;
    num_reads = unmapped_batch->num_reads;
    allocate_mapping = (array_list_t **) malloc(sizeof(array_list_t *) * num_reads);
    //printf("Region Seeker :: Process Batch with %d reads\n", num_reads); 
    #pragma omp parallel for private(num_mappings) reduction(+:total_mappings) 
    for (size_t i = 0; i < num_reads; i++) {
      allocate_mapping[i] = array_list_new(1000, 
					   1.25f, 
					   COLLECTION_MODE_ASYNCHRONIZED);
      
      num_mappings = bwt_map_exact_seeds_seq(&(unmapped_batch->seq[unmapped_batch->data_indices[i]]), 
					     input->cal_optarg_p->seed_size,
					     input->cal_optarg_p->min_seed_size,
					     input->bwt_optarg_p, input->bwt_index_p, 
					     allocate_mapping[i]);
      total_mappings += num_mappings;
      //printf("----------------->>>>>>>>>>>Regions found %d\n", num_mappings);
      
    }
    
    region_batch = (region_batch_t *) malloc(sizeof(region_batch_t));
    region_batch_init(allocate_mapping, unmapped_batch, region_batch);
          
    region_item = list_item_new(0, item->type, region_batch);

    list_item_free(item);

    if (time_on) { timing_stop(REGION_SEEKER, 0, timing_p); }
    list_insert_item(region_item, input->region_list_p);
    //printf("Region Seeker Processing batch finish!\n");
  }
  
  list_decr_writers(input->region_list_p);
 
  if (statistics_on) { 
    statistics_set(REGION_SEEKER_ST, 0, num_batches, statistics_p); 
    statistics_set(REGION_SEEKER_ST, 1, total_mappings, statistics_p); 
  }
 
  printf("region_seeker_server: END\n");
}

//====================================================================================
// apply_seeding
//====================================================================================

void apply_seeding(region_seeker_input_t* input, aligner_batch_t *batch) {

  //  printf("START: apply_seeding\n"); 

  char *seq;
  list_t *list = NULL;
  size_t index, num_mappings;
  fastq_batch_t *fq_batch = batch->fq_batch;
  size_t seed_size = input->cal_optarg_p->seed_size;
  size_t min_seed_size = input->cal_optarg_p->min_seed_size;
  size_t num_seqs = batch->num_targets;

  size_t num_outputs = 0;
  size_t *outputs = (size_t *) calloc(num_seqs, sizeof(size_t));

  // set to zero
  batch->num_done = 0;
  batch->num_to_do = 0;

  // omp parallel for !!
  for (size_t i = 0; i < num_seqs; i++) {
      
    index = batch->targets[i];
    list = batch->mapping_lists[index];

    seq = &(fq_batch->seq[fq_batch->data_indices[index]]);
    num_mappings = bwt_map_exact_seeds_seq(seq, seed_size,min_seed_size,
					   input->bwt_optarg_p, input->bwt_index_p, 
					   list);
    if (num_mappings > 0) {
      outputs[num_outputs++] = index;
      batch->num_to_do += num_mappings;
    }
    //printf("\tSEED  : read %d (%d items): %s\n", 
    //	   index, num_mappings, &(fq_batch->header[fq_batch->header_indices[index]]));

    batch->num_done += 2 * (strlen(seq) / seed_size);
    batch->num_to_do += num_mappings;
  }

  // update batch
  batch->num_allocated_targets = num_seqs;
  batch->num_targets = num_outputs;
  if (batch->targets != NULL) free(batch->targets);
  batch->targets = outputs;

  // update counter
  thr_seeding_items[omp_get_thread_num()] += batch->num_done;
  
  //  printf("END: apply_seeding, (seeding %d reads)\n", num_outputs);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
