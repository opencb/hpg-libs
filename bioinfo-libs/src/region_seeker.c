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
  list_item_t *item_p = NULL;
  list_item_t *cal_item_p = NULL;
  fastq_batch_t *unmapped_batch_p;
  size_t num_reads;
  array_list_t **allocate_mapping_p;
  cal_batch_t *cal_batch_p;
  size_t num_mappings, total_mappings = 0, num_batches = 0;
  size_t num_threads = input_p->region_threads;
  size_t chunk;
  size_t total_reads = 0;

  omp_set_num_threads(num_threads);
  
  while ( (item_p = list_remove_item(input_p->unmapped_read_list_p)) != NULL ) {

    //printf("Region Seeker Processing batch...\n");
    num_batches++;
    if (time_on) { timing_start(REGION_SEEKER, 0, timing_p); }
    
    unmapped_batch_p = (fastq_batch_t *)item_p->data_p;
    num_reads = unmapped_batch_p->num_reads;
    total_reads += num_reads;
    allocate_mapping_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
    
    if (input_p->gpu_enable) {
      //******************************* GPU PROCESS *********************************//
      for (size_t i = 0; i < num_reads; i++) {
	allocate_mapping_p[i] = array_list_new(1000, 
					       1.25f, 
					       COLLECTION_MODE_ASYNCHRONIZED);
      }
      #ifdef HPG_GPU
      num_mappings = bwt_map_exact_seed_batch_gpu(unmapped_batch_p,
						  input_p->bwt_optarg_p, 
						  input_p->cal_optarg_p,
						  input_p->bwt_index_p,
						  input_p->gpu_context,
						  allocate_mapping_p);
      #endif
      //****************************************************************************//
    } else {

      //******************************* CPU PROCESS *********************************//
      //printf("Region Seeker :: Process Batch with %d reads\n", num_reads); 
      chunk = MAX(1, num_reads/(num_threads*10));
      
      //printf("Region Seeker :: Process Batch with %d reads\n", num_reads);
      #pragma omp parallel for private(num_mappings) reduction(+:total_mappings) schedule(dynamic, chunk)
      //#pragma omp parallel for private(num_mappings) reduction(+:total_mappings) schedule(static)
      for (size_t i = 0; i < num_reads; i++) {
	//printf("Threads region zone: %d\n", omp_get_num_threads());
	
	allocate_mapping_p[i] = array_list_new(1000, 
					       1.25f, 
					       COLLECTION_MODE_ASYNCHRONIZED);
	
	num_mappings = bwt_map_exact_seeds_seq(&(unmapped_batch_p->seq[unmapped_batch_p->data_indices[i]]), 
					       input_p->cal_optarg_p->seed_size,
					       input_p->cal_optarg_p->min_seed_size,
					       input_p->bwt_optarg_p, input_p->bwt_index_p, allocate_mapping_p[i]);
	
	total_mappings += num_mappings;
	//printf("----------------->>>>>>>>>>>Regions found %d\n", num_mappings);      
      }
      //****************************************************************************//
    
    }

    cal_batch_p = cal_batch_new(allocate_mapping_p, unmapped_batch_p);
      
    list_item_free(item_p);
    cal_item_p = list_item_new(0, 0, cal_batch_p);
    //region_batch_free(region_batch_p);    

    if (time_on) { timing_stop(REGION_SEEKER, 0, timing_p); }
    
    list_insert_item(cal_item_p, input_p->region_list_p);
    //printf("Region Seeker Processing batch finish!\n");

  } //End of while
  
  list_decr_writers(input_p->region_list_p);
 
  if (statistics_on) { 
    statistics_set(REGION_SEEKER_ST, 0, num_batches, statistics_p); 
    statistics_set(REGION_SEEKER_ST, 1, total_reads, statistics_p); 
  }
 
  printf("region_seeker_server: END\n");
  
}

//====================================================================================
// apply_seeding
//====================================================================================

void apply_seeding(region_seeker_input_t* input, aligner_batch_t *batch) {

  //  printf("START: apply_seeding\n"); 

  char *seq;
  array_list_t *list = NULL;
  size_t index, num_mappings;
  fastq_batch_t *fq_batch = batch->fq_batch;
  size_t num_seeds = input->cal_optarg_p->num_seeds;
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
    //printf("region_seeker.c: apply_seeding: list #%i size = %i\n", i, array_list_size(list));
    
    seq = &(fq_batch->seq[fq_batch->data_indices[index]]);
    /*
    num_mappings = bwt_map_exact_seeds_seq(seq, seed_size, min_seed_size,
					   input->bwt_optarg_p, input->bwt_index_p, 
    					   list);
    */
    num_mappings = bwt_map_exact_seeds_seq_by_num(seq, num_seeds, seed_size, min_seed_size,
    						  input->bwt_optarg_p, input->bwt_index_p, 
    						  list);
    if (num_mappings > 0) {
      //      printf("\tregion_seeker.c: apply_seeding, setting flag to 2 for list %i\n", index);
      array_list_set_flag(2, list);
      outputs[num_outputs++] = index;
      batch->num_to_do += num_mappings;
      //    } else {
      //	if (strncmp("@rand", &(fq_batch->header[fq_batch->header_indices[index]]), 5)) {
      //	  printf("\tno seeds for read #%d: %s\n", 
      //		 index, &(fq_batch->header[fq_batch->header_indices[index]]));
      //	} 
      //   } else {
      //      printf("\tregion_seeker.c: apply_seeding: %s: list #%i, size = %i\n", &(fq_batch->header[fq_batch->header_indices[index]]), i, array_list_size(list));
    }
    //    printf("\tSEED  : read %d (%d items): %s\n", 
    //    	   index, num_mappings, &(fq_batch->header[fq_batch->header_indices[index]]));

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
  
  //  printf("region_seeker.c: apply_seeding: num_outputs = %i\n", num_outputs);

  //  printf("END: apply_seeding, (seeding %d reads)\n", num_outputs);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
