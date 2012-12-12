#include "region_seeker.h"

#ifdef HPG_GPU
   void region_seeker_input_init(list_t *unmapped_read_list, cal_optarg_t *cal_optarg, 
			         bwt_optarg_t *bwt_optarg, bwt_index_t *bwt_index, 
			         list_t* region_list, unsigned int region_threads, 
			         unsigned int gpu_enable, gpu_context_t *gpu_context,
			         region_seeker_input_t *input_p) {

  input_p->unmapped_read_list_p = unmapped_read_list;
  input_p->cal_optarg_p = cal_optarg;
  input_p->bwt_optarg_p = bwt_optarg;
  input_p->bwt_index_p = bwt_index;
  input_p->region_list_p = region_list;
  input_p->region_threads = region_threads;
  input_p->gpu_enable = gpu_enable;
  input_p->gpu_context = gpu_context;
}
#else 
   void region_seeker_input_init(list_t *unmapped_read_list, cal_optarg_t *cal_optarg, 
			         bwt_optarg_t *bwt_optarg, bwt_index_t *bwt_index, 
			         list_t* region_list, unsigned int region_threads, 
			         unsigned int gpu_enable,
			         region_seeker_input_t *input_p) {

  input_p->unmapped_read_list_p = unmapped_read_list;
  input_p->cal_optarg_p = cal_optarg;
  input_p->bwt_optarg_p = bwt_optarg;
  input_p->bwt_index_p = bwt_index;
  input_p->region_list_p = region_list;
  input_p->region_threads = region_threads;
  input_p->gpu_enable = gpu_enable;
}
#endif
//--------------------------------------------------------------------------------------

void region_seeker_server(region_seeker_input_t *input_p){  
  printf("region_seeker_server(%d): START\n", omp_get_thread_num());  
  list_item_t *item = NULL;
  mapping_batch_t *mapping_batch;
  size_t num_reads;
  array_list_t **allocate_mapping_p;
  cal_batch_t *cal_batch_p;
  size_t num_mappings, total_mappings = 0, num_batches = 0;
  size_t num_threads = input_p->region_threads;
  size_t chunk;
  size_t total_reads = 0;
  size_t targets = 0, num_targets;
 
  omp_set_num_threads(num_threads);
  
  while ( (item = list_remove_item(input_p->unmapped_read_list_p)) != NULL ) {

    //printf("Region Seeker Processing batch...\n");
    num_batches++;
    if (time_on) { timing_start(REGION_SEEKER, 0, timing_p); }
    
    mapping_batch = (mapping_batch_t *)item->data_p;
    num_targets = mapping_batch->num_targets;
    total_reads += num_targets;
    //allocate_mapping_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
    
    //if (input_p->gpu_enable) {
    //******************************* GPU PROCESS *********************************//
    /*for (size_t i = 0; i < num_reads; i++) {
      allocate_mapping_p[i] = array_list_new(1000, 
      1.25f, 
      COLLECTION_MODE_ASYNCHRONIZED);
      }*/
    //#ifdef HPG_GPU
    /*num_mappings = bwt_map_exact_seed_batch_gpu(unmapped_batch_p,
      input_p->bwt_optarg_p, 
      input_p->cal_optarg_p,
      input_p->bwt_index_p,
      input_p->gpu_context,
      allocate_mapping_p);*/
    //#endif
    //****************************************************************************//
    //} else {
    //******************************* CPU PROCESS *********************************//
    //printf("Region Seeker :: Process Batch with %d reads\n", num_targets); 
    chunk = MAX(1, num_targets/(num_threads*10));
    
    #pragma omp parallel for private(num_mappings) reduction(+:total_mappings) schedule(dynamic, chunk)
    for (size_t i = 0; i < num_targets; i++) {
      //printf("Threads region zone: %d\n", omp_get_num_threads());
      fastq_read_t *read = array_list_get(mapping_batch->targets[i], mapping_batch->fq_batch);
      num_mappings = bwt_map_exact_seeds_seq(read->sequence, 
					     input_p->cal_optarg_p->seed_size,
					     input_p->cal_optarg_p->min_seed_size,
					     input_p->bwt_optarg_p, 
					     input_p->bwt_index_p, 
					     mapping_batch->mapping_lists[mapping_batch->targets[i]]);
      
      total_mappings += num_mappings;
    } 
      //****************************************************************************//
      //}
    
    if (time_on) { timing_stop(REGION_SEEKER, 0, timing_p); }
    list_insert_item(item, input_p->region_list_p);
    //printf("Region Seeker Processing batch finish!\n");

  } //End of while
  
  list_decr_writers(input_p->region_list_p);
  /*
  if (statistics_on) { 
    statistics_set(REGION_SEEKER_ST, 0, num_batches, statistics_p); 
    statistics_set(REGION_SEEKER_ST, 1, total_reads, statistics_p); 
  }
  */
  printf("region_seeker_server: END\n");
  
}

//====================================================================================
// apply_seeding
//====================================================================================

void apply_seeding(region_seeker_input_t* input, mapping_batch_t *batch) {

  //  printf("START: apply_seeding\n"); 

  array_list_t *list = NULL;
  size_t read_index, num_mappings;

  size_t min_num_seeds = input->cal_optarg_p->min_num_seeds;
  size_t max_num_seeds = input->cal_optarg_p->max_num_seeds;
  size_t seed_size = input->cal_optarg_p->seed_size;
  size_t min_seed_size = input->cal_optarg_p->min_seed_size;

  fastq_read_t *fq_read;
  array_list_t *fq_batch = batch->fq_batch;

  size_t num_targets = batch->num_targets;
  size_t *targets = batch->targets;
  size_t new_num_targets = 0;
  //  size_t *new_targets = (size_t *) calloc(num_targets, sizeof(size_t));

  // set to zero
  batch->num_to_do = 0;

  // omp parallel for !!
  for (size_t i = 0; i < num_targets; i++) {
      
    read_index = targets[i];
    list = batch->mapping_lists[read_index];
    fq_read = array_list_get(read_index, fq_batch);
    //printf("region_seeker.c: apply_seeding: list #%i size = %i\n", i, array_list_size(list));
    
    num_mappings = bwt_map_exact_seeds_seq_by_num(fq_read->sequence, min_num_seeds, max_num_seeds, 
						  seed_size, min_seed_size,
    						  input->bwt_optarg_p, input->bwt_index_p, 
    						  list);
    if (num_mappings > 0) {
      //      printf("\tregion_seeker.c: apply_seeding, setting flag to 2 for list %i\n", index);
      array_list_set_flag(2, list);
      targets[new_num_targets++] = read_index;

      batch->num_to_do += num_mappings;
    }
  }
    
  // update batch targets
  batch->num_targets = new_num_targets;
  //  batch->num_allocated_targets = num_targets;
  //  if (batch->targets) free(batch->targets);
  //  batch->targets = new_targets;
  
  // update counter
  //thr_seeding_items[omp_get_thread_num()] += batch->num_done;
  
  //  printf("region_seeker.c: apply_seeding: num_outputs = %i\n", num_outputs);
  //  printf("END: apply_seeding, (seeding %d reads)\n", num_outputs);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
