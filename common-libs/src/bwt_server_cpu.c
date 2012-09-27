#include "bwt_server.h"


void bwt_server_cpu(bwt_server_input_t* input_p){
    
    printf("bwt_server_cpu(%d): START\n", omp_get_thread_num()); 
    list_item_t *item_p = NULL;
    list_item_t *write_item_p = NULL;
    list_item_t *unmapped_item_p = NULL;
    fastq_batch_t *fastq_batch_p;
    fastq_batch_t *unmapped_batch_p;
    array_list_t *mappings;
    size_t num_mappings;
    unsigned int write_size = input_p->batch_size;
    list_t *write_list_p = input_p->write_list_p;
    list_t *unmapped_read_list_p = input_p->unmapped_read_list_p;
    size_t num_mappings_tot = 0;
    size_t total_reads = 0;
    size_t reads_no_mapped = 0; 
    write_batch_t* write_batch_p = write_batch_new(write_size, MATCH_FLAG);
    size_t num_batches = 0;
    while ( (item_p = list_remove_item(input_p->read_list_p)) != NULL ) {
	num_batches++;
	if (time_on) { timing_start(BWT_SERVER, 0, timing_p); }
	
	fastq_batch_p = (fastq_batch_t *)item_p->data_p;
	//printf("BWT Processing batch...\n");
	unmapped_batch_p = fastq_batch_new(input_p->batch_size);

	//printf("Batch Init done. Processing batch ...\n");
	mappings = array_list_new(100000, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
	//start_timer(start_time);
	//printf("\tCall function process\n");
	num_mappings = bwt_map_inexact_batch(fastq_batch_p, input_p->bwt_optarg_p, input_p->bwt_index_p, unmapped_batch_p, mappings);
	//printf("\tEnd call\n");
	num_mappings_tot += num_mappings;
	total_reads += fastq_batch_p->num_reads;
	reads_no_mapped += unmapped_batch_p->num_reads;
	//Results
	//printf("Process Batch (bwt_server_cpu): (%d)Mappings - (%d)Unmappings\n", num_mappings, unmapped_batch_p->num_reads);
	for(int i = 0; i < num_mappings; i++){
	  
	  if ( write_batch_p->size >= write_batch_p->allocated_size - 1) {
	    write_item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
	    if (time_on) { timing_stop(BWT_SERVER, 0, timing_p); }
	    list_insert_item(write_item_p, write_list_p);
	    if (time_on) { timing_start(BWT_SERVER, 0, timing_p); }
	    
	    write_batch_p = write_batch_new(write_size, MATCH_FLAG);
	  }
	  ((alignment_t **)write_batch_p->buffer_p)[write_batch_p->size] = (alignment_t *)array_list_get(i, mappings);
	  write_batch_p->size++;
	}
	
	array_list_free(mappings, NULL);
	list_item_free(item_p);
	fastq_batch_free(fastq_batch_p);
	//	fastq_batch_free(unmapped_batch_p);

	unmapped_item_p = list_item_new(0, WRITE_ITEM, unmapped_batch_p); 
	if (time_on) { timing_stop(BWT_SERVER, 0, timing_p); }
	list_insert_item(unmapped_item_p, unmapped_read_list_p);
	
	//printf("BWT Processing batch finish! %d reads no mapped\n", unmapped_batch_p->num_reads);
	
    }
    
    if (write_batch_p != NULL) {
      if (write_batch_p->size > 0) {
	write_item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
	list_insert_item(write_item_p, write_list_p);
      } else {
	write_batch_free(write_batch_p);
      }
    }
    
    list_decr_writers(write_list_p);
    list_decr_writers(unmapped_read_list_p);
    
    if (statistics_on) { 
      statistics_set(BWT_SERVER_ST, 0, num_batches, statistics_p);
      statistics_set(BWT_SERVER_ST, 1, total_reads, statistics_p); 
      statistics_set(BWT_SERVER_ST, 2, total_reads - reads_no_mapped, statistics_p); 
      statistics_set(BWT_SERVER_ST, 3, reads_no_mapped, statistics_p); 
      statistics_set(BWT_SERVER_ST, 4, num_mappings_tot, statistics_p); 
      
      statistics_add(TOTAL_ST, 1, total_reads - reads_no_mapped, statistics_p); 
    }
    
    printf("bwt_server_cpu (Total reads process %i, Reads unmapped %i, Total mappings %i): END\n", total_reads, reads_no_mapped, num_mappings_tot); 
}

//====================================================================================
// bwt_server_input functions: init
//====================================================================================

void bwt_server_input_init(list_t* read_list_p, unsigned int batch_size, bwt_optarg_t *bwt_optarg_p, bwt_index_t *bwt_index_p, list_t* write_list_p, unsigned int write_size, list_t* unmapped_read_list_p, bwt_server_input_t* input_p) {
  input_p->read_list_p = read_list_p;
  input_p->batch_size = batch_size;
  input_p->bwt_optarg_p = bwt_optarg_p;
  input_p->write_list_p = write_list_p;
  input_p->write_size = write_size;
  input_p->bwt_index_p = bwt_index_p;
  input_p->unmapped_read_list_p = unmapped_read_list_p;
  //extern int num_of_chromosomes;

}

