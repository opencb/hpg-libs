#include "region_seeker.h"

void region_seeker_input_init(list_t *unmapped_read_list_p, cal_optarg_t *cal_optarg_p, bwt_optarg_t *bwt_optarg_p, bwt_index_t *bwt_index_p, list_t* region_list_p, unsigned int region_threads, region_seeker_input_t *input_p){
  input_p->unmapped_read_list_p = unmapped_read_list_p;
  input_p->cal_optarg_p = cal_optarg_p;
  input_p->bwt_optarg_p = bwt_optarg_p;
  input_p->bwt_index_p = bwt_index_p;
  input_p->region_list_p = region_list_p;
  input_p->region_threads = region_threads;
}

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
    //printf("Region Seeker :: Process Batch with %d reads\n", num_reads); 
    chunk = MAX(1, num_reads/(num_threads*10));

    //printf("Region Seeker :: Process Batch with %d reads\n", num_reads);
    #pragma omp parallel for private(num_mappings) reduction(+:total_mappings) schedule(dynamic, chunk)
    for (size_t i = 0; i < num_reads; i++) {
      //printf("Threads region zone: %d\n", omp_get_num_threads());

      allocate_mapping_p[i] = array_list_new(1000, 
					     1.25f, 
					     COLLECTION_MODE_SYNCHRONIZED);
      
      num_mappings = bwt_map_exact_seeds_seq(&(unmapped_batch_p->seq[unmapped_batch_p->data_indices[i]]), 
					     input_p->cal_optarg_p->seed_size,
					     input_p->cal_optarg_p->min_seed_size,
					     input_p->bwt_optarg_p, input_p->bwt_index_p, allocate_mapping_p[i]);
      total_mappings += num_mappings;
      //printf("----------------->>>>>>>>>>>Regions found %d\n", num_mappings);
      
    }
    
    cal_batch_p = cal_batch_new(allocate_mapping_p, unmapped_batch_p);
      
    
    list_item_free(item_p);
    cal_item_p = list_item_new(0, 0, cal_batch_p);
    //region_batch_free(region_batch_p);    

    if (time_on) { timing_stop(REGION_SEEKER, 0, timing_p); }
    
    list_insert_item(cal_item_p, input_p->region_list_p);
    //printf("Region Seeker Processing batch finish!\n");

  }
  
  list_decr_writers(input_p->region_list_p);
 
  if (statistics_on) { 
    statistics_set(REGION_SEEKER_ST, 0, num_batches, statistics_p); 
    statistics_set(REGION_SEEKER_ST, 1, total_reads, statistics_p); 
  }
 
  printf("region_seeker_server: END\n");
  
}
