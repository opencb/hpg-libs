#include "bwt_server.h"


void bwt_server_cpu(bwt_server_input_t* input){
    
    printf("bwt_server_cpu(%d): START\n", omp_get_thread_num()); 
    list_item_t *item = NULL;
    list_item_t *write_item = NULL;
    list_item_t *unmapped_item = NULL;
    fastq_batch_t *fastq_batch;
    fastq_batch_t *unmapped_batch;
    array_list_t *mappings;
    size_t num_mappings;
    size_t write_size = input->batch_size;
    list_t *list = NULL;
    list_t *pair_list = input->pair_list;
    list_t *write_list = input->write_list;
    list_t *unmapped_read_list = input->unmapped_read_list;
    size_t num_mappings_tot = 0;
    size_t total_reads = 0;
    size_t reads_no_mapped = 0; 
    write_batch_t* write_batch = write_batch_new(write_size, MATCH_FLAG);
    size_t num_batches = 0;

    while ( (item = list_remove_item(input->read_list)) != NULL ) {
      
        if (item->type & PAIR1_FLAG || item->type & PAIR2_FLAG) {
	  list = pair_list;
	  printf("bwt_server inserting into pair list\n");
	} else {
	  list = write_list;
	  printf("bwt_server inserting into write list\n");
	}

	num_batches++;
	if (time_on) { timing_start(BWT_SERVER, 0, timing_p); }
	
	fastq_batch = (fastq_batch_t *) item->data_p;
	//printf("BWT Processing batch...\n");
	unmapped_batch = fastq_batch_new(input->batch_size);

	//printf("Batch Init done. Processing batch ...\n");
	mappings = array_list_new(100000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	//start_timer(start_time);
	//printf("\tCall function process\n");
	num_mappings = bwt_map_inexact_batch(fastq_batch, input->bwt_optarg, 
					     input->bwt_index, unmapped_batch, mappings);
	//printf("\tEnd call\n");
	num_mappings_tot += num_mappings;
	total_reads += fastq_batch->num_reads;
	reads_no_mapped += unmapped_batch->num_reads;

	//Results
	//printf("Process Batch (bwt_server_cpu): (%d)Mappings - (%d)Unmappings\n", 
	//       num_mappings, unmapped_batch->num_reads);
	for(int i = 0; i < num_mappings; i++){
	  /*
	  printf("\t\t\tmapping found by bwt_server, name = %s, strand = %i, position = %i\n",
		 ((alignment_t *)array_list_get(i, mappings))->query_name,
		 ((alignment_t *)array_list_get(i, mappings))->seq_strand,
		 ((alignment_t *)array_list_get(i, mappings))->position);
	  */
	  
	  if ( write_batch->size >= write_batch->allocated_size - 1) {
	    write_item = list_item_new(0, item->type | WRITE_ITEM_FLAG, write_batch);
	    if (time_on) { timing_stop(BWT_SERVER, 0, timing_p); }
	    list_insert_item(write_item, list);
	    if (time_on) { timing_start(BWT_SERVER, 0, timing_p); }
	    
	    write_batch = write_batch_new(write_size, MATCH_FLAG);
	  }
	  ((alignment_t **) write_batch->buffer_p)[write_batch->size] = (alignment_t *)array_list_get(i, mappings);

	  write_batch->size++;
	}
	
	// insert mapped read batch into the list to write
	if (write_batch->size > 0) {
	  write_item = list_item_new(0, item->type | WRITE_ITEM_FLAG, write_batch);
	  if (time_on) { timing_stop(BWT_SERVER, 0, timing_p); }
	  list_insert_item(write_item, list);
	  if (time_on) { timing_start(BWT_SERVER, 0, timing_p); }

	  write_batch = write_batch_new(write_size, MATCH_FLAG);
	}

	// insert unmapped read batch into the list to write
	unmapped_item = list_item_new(0, item->type, unmapped_batch); 
	if (time_on) { timing_stop(BWT_SERVER, 0, timing_p); }
	list_insert_item(unmapped_item, unmapped_read_list);
	
	array_list_free(mappings, NULL);
	list_item_free(item);
	fastq_batch_free(fastq_batch);

	//printf("BWT Processing batch finish! %d reads no mapped\n", unmapped_batch->num_reads);
    }

    // free memory
    write_batch_free(write_batch);

    // decreasing writers
    if (write_list != NULL) list_decr_writers(write_list);
    if (unmapped_read_list != NULL) list_decr_writers(unmapped_read_list);
    if (pair_list != NULL) list_decr_writers(pair_list);

    if (statistics_on) { 
      statistics_set(BWT_SERVER_ST, 0, num_batches, statistics_p);
      statistics_set(BWT_SERVER_ST, 1, total_reads, statistics_p); 
      statistics_set(BWT_SERVER_ST, 2, total_reads - reads_no_mapped, statistics_p); 
      statistics_set(BWT_SERVER_ST, 3, reads_no_mapped, statistics_p); 
      statistics_set(BWT_SERVER_ST, 4, num_mappings_tot, statistics_p); 
      
      statistics_add(TOTAL_ST, 1, total_reads - reads_no_mapped, statistics_p); 
    }
    
    printf("bwt_server_cpu (Total reads process %i, Reads unmapped %i, Total mappings %i): END\n", 
	   total_reads, reads_no_mapped, num_mappings_tot); 
}

//====================================================================================
// bwt_server_input functions: init
//====================================================================================

void bwt_server_input_init(list_t *read_list, unsigned int batch_size, 
			   bwt_optarg_t *bwt_optarg, bwt_index_t *bwt_index, 
			   list_t *pair_list, list_t *write_list, size_t write_size,  
			   list_t* unmapped_read_list, bwt_server_input_t* input) {

  input->read_list = read_list;
  input->batch_size = batch_size;
  input->bwt_optarg = bwt_optarg;
  input->write_list = write_list;
  input->write_size = write_size;
  input->pair_list = pair_list;
  input->bwt_index = bwt_index;
  input->unmapped_read_list = unmapped_read_list;
}

//====================================================================================
// apply_bwt
//====================================================================================

void bwt_map_inexact_batch_by_filter(fastq_batch_t *batch,
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index,
				     size_t num_selected, size_t *selected,
				     size_t num_lists, array_list_t **lists,
				     size_t *num_mapped, size_t* mapped,
				     size_t *num_unmapped, size_t *unmapped);

//------------------------------------------------------------------------------------

void apply_bwt(bwt_server_input_t* input, aligner_batch_t *batch) {

  //    printf("START: apply_bwt\n"); 

    size_t num_mapped_reads = 0;
    size_t *mapped_reads = (size_t *) calloc(batch->num_mapping_lists, 
					     sizeof(size_t));

    // count the number of BWT to run
    batch->num_done = batch->num_mapping_lists;
    
    // run bwt
    bwt_map_inexact_batch_by_filter(batch->fq_batch, input->bwt_optarg,
				    input->bwt_index, 
				    batch->num_mapping_lists, NULL, 
				    batch->num_mapping_lists, batch->mapping_lists,
				    &num_mapped_reads, mapped_reads,
				    &batch->num_targets, batch->targets);

    batch->num_to_do = num_mapped_reads;

    // free memory
    free(mapped_reads);

    // update counter
    thr_bwt_items[omp_get_thread_num()] += batch->num_done;

    //    printf("END: apply_bwt, (mapped reads = %d, unmapped reads = %d)\n",
    //	   num_mapped_reads, batch->num_targets);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------



    /*
    // this actions must be performed by the scheduler
    if (p->num_targets == 0) {
      p->action = 
    }
    p->all_targets = (p->num_mapping_lists == p->num_targets ? 1 : 0);
    p->action = SEEDING_ACTION;
    */
