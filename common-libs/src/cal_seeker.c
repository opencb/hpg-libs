
#include "cal_seeker.h"
#include "libs/aligners/bwt/bwt.h"

//------------------------------------------------------------------------------------
// cal_seeker_input functions: init
//------------------------------------------------------------------------------------

void cal_seeker_server(cal_seeker_input_t* input_p) {
  
  extern short int cal_type;
  
  unsigned int cal_id = omp_get_thread_num();

  printf("cal_seeker_server(%d): START\n", cal_id); 
  
  list_item_t *item_p = NULL;
  list_item_t *write_item_p = NULL;
  list_item_t *sw_item_p = NULL;
  size_t num_reads;
  array_list_t **allocate_cals_p;
  cal_batch_t *cal_batch_p;
  size_t num_cals, total_cals = 0;
  fastq_read_t *read_p;

  char *seq;
  char *header;
  char *quality;
  unsigned int seq_len, header_len;

  char *cigar_p;
  char *header_id;
  unsigned int write_size = input_p->batch_size;
  write_batch_t* write_batch_p = write_batch_new(write_size, MATCH_FLAG);
  
  list_t *write_list_p = input_p->write_list_p;
  sw_batch_t *sw_batch_p;
  fastq_read_t **allocate_reads_p; 
  size_t reads_with_cals;
  alignment_t *alignment_p;
  size_t num_batches = 0, num_reads_unmapped = 0 ;
  size_t total_reads = 0;
  
  while ( (item_p = list_remove_item(input_p->regions_list_p)) != NULL ) {
    //printf("Cal Seeker %d processing batch...\n", omp_get_thread_num());
    num_batches++;
    if (time_on) { timing_start(CAL_SEEKER, cal_id, timing_p); }
    
    cal_batch_p = (cal_batch_t *)item_p->data_p;
    num_reads = cal_batch_p->unmapped_batch_p->num_reads;
    total_reads += num_reads;
    allocate_cals_p = (array_list_t **)calloc(num_reads, sizeof(array_list_t *));
    allocate_reads_p = (fastq_read_t **)calloc(num_reads, sizeof(fastq_read_t *));
    reads_with_cals = 0;
    for (size_t i = 0; i < num_reads; i++) {
      allocate_cals_p[reads_with_cals] = array_list_new(100, 
					                1.25f, 
							COLLECTION_MODE_ASYNCHRONIZED);
      /*      if (cal_type == 0) {
		num_cals = bwt_generate_cal_list(cal_batch_p->allocate_mapping_p[i], input_p->cal_optarg_p, allocate_cals_p[reads_with_cals]);
		} else {*/
      num_cals = bwt_generate_cal_list_linkedlist(cal_batch_p->allocate_mapping_p[i], input_p->cal_optarg_p, allocate_cals_p[reads_with_cals]);
		// }
      total_cals += num_cals;
      //printf("CAL SEEKER: %d cals for read %d\n", num_cals, i);
      if(num_cals <= 0 || num_cals > MAX_CALS){
	
	//printf("READ WITH cals <= 0 or cals %d > %d!!\n", array_list_size(allocate_cals_p[reads_with_cals]), MAX_CALS);
	
	num_reads_unmapped++;
	
	alignment_p = alignment_new();
	
	seq_len = cal_batch_p->unmapped_batch_p->data_indices[i + 1] - cal_batch_p->unmapped_batch_p->data_indices[i];

	header_len = cal_batch_p->unmapped_batch_p->header_indices[i + 1] - cal_batch_p->unmapped_batch_p->header_indices[i];

	header_id = (char *)malloc(sizeof(char)*header_len);
	header_len = get_to_first_blank(&(cal_batch_p->unmapped_batch_p->header[cal_batch_p->unmapped_batch_p->header_indices[i]]), header_len, header_id);
	
	seq = (char *)calloc(sizeof(char), seq_len);
	memcpy(seq, &(cal_batch_p->unmapped_batch_p->seq[cal_batch_p->unmapped_batch_p->data_indices[i]]), seq_len);
	
	quality = (char *)calloc(sizeof(char), seq_len);
	memcpy(quality, &(cal_batch_p->unmapped_batch_p->quality[cal_batch_p->unmapped_batch_p->data_indices[i]]), seq_len);
	
	cigar_p = (char *)malloc(sizeof(char)*10);
	sprintf(cigar_p, "%dX\0", seq_len - 1);
	
	alignment_init_single_end(header_id, seq, quality, 0, 
				  0, //index->karyotype.chromosome + (idx-1) * IDMAX,
				  0, 
				  cigar_p, 1, 255, 0, 0, alignment_p);

	//printf("List %d > %d?\n", write_batch_p->size, write_batch_p->allocated_size);		
	if ( write_batch_p->size >= write_batch_p->allocated_size - 1) {
	  write_item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
	  if (time_on) { timing_stop(CAL_SEEKER, cal_id, timing_p); }
	  list_insert_item(write_item_p, write_list_p);
	  if (time_on) { timing_start(CAL_SEEKER, cal_id, timing_p); }
	  
	  write_batch_p = write_batch_new(write_size, MATCH_FLAG);
	}
	
	((alignment_t **)write_batch_p->buffer_p)[write_batch_p->size] = alignment_p;
	write_batch_p->size++;
	
	array_list_free(allocate_cals_p[reads_with_cals], cal_free);
      }else{
	//array_list_free(allocate_cals_p[reads_with_cals], cal_free);
	
	read_p = fastq_read_new(&(cal_batch_p->unmapped_batch_p->header[cal_batch_p->unmapped_batch_p->header_indices[i]]), 
				&(cal_batch_p->unmapped_batch_p->seq[cal_batch_p->unmapped_batch_p->data_indices[i]]),
			        &(cal_batch_p->unmapped_batch_p->quality[cal_batch_p->unmapped_batch_p->data_indices[i]]));
	
	allocate_reads_p[reads_with_cals] = read_p;
	reads_with_cals++;
      }
       
      //array_list_free(regions_batch_p->allocate_mapping_p[i], NULL);
  }
    
    sw_batch_p = sw_batch_new(reads_with_cals, allocate_cals_p, allocate_reads_p);
    
    //sw_batch_free(sw_batch_p);
    cal_batch_free(cal_batch_p);
    free(item_p);
    sw_item_p = list_item_new(0, 0, sw_batch_p);
  
    if (time_on) { timing_stop(CAL_SEEKER, cal_id, timing_p); }
    
    list_insert_item(sw_item_p, input_p->sw_list_p); 
    //printf("Cal Seeker process batch finish!\n");
  }
  
  if (write_batch_p != NULL) {
    if (write_batch_p->size > 0) {
      write_item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
      list_insert_item(write_item_p, write_list_p);
    } else {
      write_batch_free(write_batch_p);
    }
  }

  
  list_decr_writers(input_p->sw_list_p);
  list_decr_writers(input_p->write_list_p);
  
  if (statistics_on) { 
   statistics_add(CAL_SEEKER_ST, 0, num_batches, statistics_p); 
   statistics_add(CAL_SEEKER_ST, 1, total_reads, statistics_p); 
   statistics_add(CAL_SEEKER_ST, 2, num_reads_unmapped, statistics_p); 
   
   statistics_add(TOTAL_ST, 2, num_reads_unmapped, statistics_p); 
  }
 
  printf("cal_seeker_server (%d reads unmapped): END\n", num_reads_unmapped); 
  // free memory for mapping list
  
}

//====================================================================================
// apply_caling
//====================================================================================
int unmapped_by_max_cals_counter[100];

void apply_caling(cal_seeker_input_t* input, aligner_batch_t *batch) {

  //  printf("START: apply_caling\n"); 
  int tid = omp_get_thread_num();

  array_list_t *list = NULL;
  size_t index, num_cals;
  size_t num_seqs = batch->num_targets;

  size_t num_outputs = 0;
  size_t *outputs = (size_t *) calloc(num_seqs, sizeof(size_t));

  // set to zero
  batch->num_done = batch->num_to_do;
  batch->num_to_do = 0;

  for (size_t i = 0; i < num_seqs; i++) {

    index = batch->targets[i];
    
    if (list == NULL) {
      list = array_list_new(1000, 
			    1.25f, 
			    COLLECTION_MODE_ASYNCHRONIZED);
    }


    // previous version
    //num_cals = bwt_generate_cal_list(batch->mapping_lists[index], 
    //input->cal_optarg,
    //				     list);

    // optimized version
    num_cals = bwt_generate_cal_list_linkedlist(batch->mapping_lists[index], 
						input->cal_optarg,
						list);
    /*
    printf("--------------------------------------------------------------\n");
    printf("cal_seeker.c: read %d (input %d seed-mappings -> output %d CALs): %s\n", 
	   index, array_list_size(batch->mapping_lists[index]), num_cals,
	   &(batch->fq_batch->header[batch->fq_batch->header_indices[index]]));
    printf("--------------------------------------------------------------\n");
    */
    //    printf("cal_seeker.c: num_cals = %d (MAX = %d)\n", num_cals, MAX_CALS);

    if (num_cals > 0 && num_cals <= MAX_CALS) {

      batch->num_to_do += num_cals;

      outputs[num_outputs++] = index;
      
      // we have to free the region list
      array_list_free(batch->mapping_lists[index], region_free);
      batch->mapping_lists[index] = list;
      list = NULL;
    } else {
      if (strncmp("@rand", &(batch->fq_batch->header[batch->fq_batch->header_indices[index]]), 5)) {
	unmapped_by_max_cals_counter[tid]++;
      }
      array_list_clear(batch->mapping_lists[index], region_free);
      array_list_clear(list, cal_free);
    }
  } // end for 0 ... num_seqs

  // update batch
  batch->num_allocated_targets = num_seqs;
  batch->num_targets = num_outputs;
  if (batch->targets != NULL) free(batch->targets);
  batch->targets = outputs;
  
  // update counter
  thr_cal_items[tid] += batch->num_done;

  // free memory
  if (list != NULL) array_list_free(list, NULL);
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

void cal_seeker_input_init(list_t *regions_list_p, cal_optarg_t *cal_optarg_p, 
			   list_t* write_list_p, unsigned int write_size, 
			   list_t *sw_list_p, list_t *pair_list_p, 
			   cal_seeker_input_t *input_p){
  input_p->regions_list_p = regions_list_p;
  input_p->cal_optarg_p = cal_optarg_p;
  input_p->batch_size = write_size;
  input_p->sw_list_p = sw_list_p;
  input_p->pair_list_p = pair_list_p;
  input_p->write_list_p = write_list_p;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
