#include "pair_server.h"

//------------------------------------------------------------------------------------
// this functions writes reads to disk, reads came from
// a list
//------------------------------------------------------------------------------------

//unsigned long alignment_hash_code(void *p);
//int alignment_compare(void *p1, void *p2);

void pair_server(pair_server_input_t* input) {

  cp_hashtable *hashtable = cp_hashtable_create_by_mode(COLLECTION_MODE_NOSYNC, 1000, 
							cp_hash_istring, 
							(cp_compare_fn) strcasecmp);
  
  
  printf("pair_server (%i): START\n", omp_get_thread_num());

  size_t total_pairs = 0, num_alignments, num_reads, num_cals;
  alignment_t **alignments, *alig;
  cal_t *cal;
		
  list_t* pair_list = input->pair_list;
  list_t* sw_list = input->sw_list;
  list_t* write_list = input->write_list;

  list_item_t *pair_item = NULL;
  write_batch_t* write_batch = NULL;
  sw_batch_t* sw_batch = NULL;


  // main loop
  while ( (pair_item = list_remove_item(pair_list)) != NULL ) {

    if (pair_item->type & WRITE_ITEM_FLAG) {

      // batch for writer
      printf("\t\tpair_server, batch to batch writer\n");

      write_batch = (write_batch_t*) pair_item->data_p;

      if (write_batch->flag == MATCH_FLAG) {
	printf("\t\t\t\tmatch flag !!!\n");      
      } else if (write_batch->flag == MISMATCH_FLAG) {
	printf("\t\t\t\tmismatch flag !!!\n");      
      } else {
	printf("\t\t\t\tunknown flag !!!\n");      
      }
      
      num_alignments = write_batch->size;
      alignments = (alignment_t **) write_batch->buffer_p;
      for (size_t i = 0; i < num_alignments; i++) {
	alig = alignments[i];
	printf("\t\t\t\t\t%s (chr = %i, strand = %d, pos = %d, cigar = %s)\n", 
	       alig->query_name, alig->chromosome, alig->seq_strand, 
	       alig->position, alig->cigar);
      }

    } else if (pair_item->type & SW_ITEM_FLAG) {

      // batch for Smith-Waterman server
      printf("\t\tpair_server, batch to sw server\n");      

      sw_batch = (sw_batch_t*) pair_item->data_p;

      num_reads = sw_batch->num_reads;
      for (size_t i = 0; i < num_reads; i++) {
	num_cals = array_list_size(sw_batch->allocate_cals_p[i]);
	printf("\t\t\t\t\t%s (num_cals = %d)\n", 
	       sw_batch->allocate_reads_p[i]->id, num_cals);
	for (size_t j = 0; j < num_cals; j++) {
	  cal = array_list_get(j, sw_batch->allocate_cals_p[i]);

	  printf("\t\t\t\t\t\tcal %d: chr = %d, strand = %d, start = %d, end = %d)\n", 
		 j, cal->chromosome_id, cal->strand, cal->start, cal->end);

	}
      }


    } else {
      printf("\t\tpair_server, batch to unknown !!!\n");      
    }
    /*    

    pairs = 0;
    if ( (item->type & PAIR1_FLAG) || (item->type & PAIR2_FLAG) ) {
      pairs = 1;
    }
    
    if (time_on) { timing_start(BATCH_WRITER, 0, timing_p); }


    if (batch->flag == SPLICE_FLAG) { 
      fwrite((char *) batch->buffer_p, batch->size, 1, splice_fd);
    } else if (batch->flag == MATCH_FLAG || batch->flag == MISMATCH_FLAG) {

      alignments = (alignment_t **) batch->buffer_p;
      num_alignments = batch->size;


      if ( pairs ) {
	process_pair((item->type & PAIR1_FLAG ? 1 : 2), 
		     alignments, num_alignments, 
		     bam_file, hashtable);
      } else {

	total_mappings += num_alignments;
	for (size_t i = 0; i < num_alignments; i++) {
	  //alignment_print(buffer_p[i]);
	  //sprintf(alignments[i]->cigar, "100=");
	  //printf("+++++ cigar = %s\n", alignments[i]->cigar);
	  
	  bam1 = convert_to_bam(alignments[i], 33);
	  bam_fwrite(bam1, bam_file);
	  bam_destroy1(bam1);
	  alignment_free(alignments[i]);
	}
	}
    }
    */

    if (time_on) { timing_stop(BATCH_WRITER, 0, timing_p); }
  } // end of batch loop
  
    // free memory
  //if (write_batch != NULL) write_batch_free(write_batch);
  //if (sw_batch != NULL) write_batch_free(sw_batch);
  if (pair_item != NULL) list_item_free(pair_item);
  
  // decreasing writers
  if (sw_list != NULL) list_decr_writers(sw_list);
  if (write_list != NULL) list_decr_writers(write_list);
  
  printf("pair_server (Total pairs %d): END\n", total_pairs);
}

//------------------------------------------------------------------------------------
/*
unsigned long alignment_hash_code(void *p) {
  alignment_t *alignment = (alignment_t *) p;

  char name[strlen(alignment->query_name)];
  strcpy(name, alignment->query_name);
  char *s = strrchr(name, '/');
  if (s != NULL) {
    *s = '\0';
  }

  printf("name = %s\n", name);

  return cp_hash_istring(name);
}    

//------------------------------------------------------------------------------------

int alignment_compare(void *p1, void *p2) {

  alignment_t *alignment1 = (alignment_t *) p1;
  alignment_t *alignment2 = (alignment_t *) p2;

  char name1[strlen(alignment1->query_name)];
  strcpy(name1, alignment1->query_name);
  char *s = strrchr(name1, '/');
  if (s != NULL) {
    *s = '\0';
  }

  char name2[strlen(alignment2->query_name)];
  strcpy(name2, alignment2->query_name);
  *s = strrchr(name2, '/');
  if (s != NULL) {
    *s = '\0';
  }

  return cp_hash_compare_istring(name1, name2);
}
*/
//------------------------------------------------------------------------------------
/*
void process_pair(int pair_id, alignment_t **alignments, size_t num_alignments,
		  bam_file_t *bam_file, cp_hashtable *hashtable) {

  bam1_t* bam1;

  printf("******************* arriving data from pair %i\n", pair_id);

  alignment_t *value;
  char *s;
  for (size_t i = 0; i < num_alignments; i++) {
    char name[strlen(alignments[i]->query_name)];
    strcpy(name, alignments[i]->query_name);
    if ( (s = strrchr(name, '/')) != NULL) {
      *s = '\0';
    }
    printf("\t\t %s\n", name);
    
    value = (alignment_t *) cp_hashtable_get(hashtable, (void *) name);
    if (value == NULL) {
      printf("\t\t\tNOT FOUND IN HASHTABLE, then insert it...\n");
      cp_hashtable_put(hashtable, (void *) name, (void *) alignments[i]);
    } else {
      printf("\t\t\tFOUND IN HASHTABLE !!!\n");

      if (pair_id == 1) {
	alignment_update_paired_end(alignments[i], value);
      } else {
	alignment_update_paired_end(value, alignments[i]);
      }

      printf("\t\t\t\t%s   :::    %s\n", value->query_name, alignments[i]->query_name);

      bam1 = convert_to_bam(alignments[i], 33);
      bam_fwrite(bam1, bam_file);
      bam_destroy1(bam1);
      alignment_free(alignments[i]);

      bam1 = convert_to_bam(value, 33);
      bam_fwrite(bam1, bam_file);
      bam_destroy1(bam1);
      alignment_free(value);

      cp_hashtable_remove(hashtable, name);
    }
  }
}

//------------------------------------------------------------------------------------
*/

//------------------------------------------------------------------------------------

inline void remove_items(size_t *valid_items, array_list_t *list) {
  void *item;
  size_t num_items = array_list_size(list);
  int flag = array_list_get_flag(list);

  printf("list of %i items\n", num_items);
  for (int k = num_items - 1; k >= 0; k--) {
    if (valid_items[k] == 0) {
      printf("%i, ", k);
      item = (void *) array_list_remove_at(k, list);
      if (flag == 1) {
	alignment_free((alignment_t *) item);
      } else {
	cal_free((cal_t *) item);
      }
    }
  }
  printf("\n");
}

//====================================================================================
// apply_pair
//====================================================================================

void apply_pair(pair_server_input_t* input, aligner_batch_t *batch) {

  //  printf("START: apply_pair\n"); 

  char *seq;
  list_t *list = NULL;
  size_t index, num_mappings;
  fastq_batch_t *fq_batch = batch->fq_batch;

  int pair_mode = input->pair_mng->pair_mode;
  size_t min_distance = input->pair_mng->min_distance;
  size_t max_distance = input->pair_mng->max_distance;

  size_t num_reads = batch->num_mapping_lists; // it must be equal to fq_batch->num_reads
  size_t num_targets = batch->num_targets;
  printf("total seqs = %d, num_seqs to process = %d\n", num_reads, batch->num_targets);

  int flag1, flag2;
  array_list_t *list1, *list2;
  size_t num_items1 = 0, num_items2 = 0;
  size_t num_allocated_items1 = 0, num_allocated_items2 = 0;
  alignment_t *alig;
  cal_t *cal;
   // size_t num_outputs = 0;
  //  size_t *outputs = (size_t *) calloc(num_seqs, sizeof(size_t));

  size_t end1, start2, distance;
  int pair_found, chr1, chr2, strand1, strand2;
  
  size_t *mapped1 = NULL, *mapped2 = NULL;
  int mapped1_counter, mapped2_counter;

  printf("pair_server.c:apply_pair: pair_mode = %i, min_distance = %lu, max_distance = %lu\n",
	 pair_mode, min_distance, max_distance);

  for (size_t i = 0; i < num_reads; i += 2) {

    list1 = batch->mapping_lists[i];
    list2 = batch->mapping_lists[i + 1];

    flag1 = array_list_get_flag(list1);
    flag2 = array_list_get_flag(list2);

    num_items1 = array_list_size(list1);
    num_items2 = array_list_size(list2);

    if (num_items1 > 1 && num_items2 > 1) {

      // allocated memory for items from list #1
      mapped1_counter = 0;
      if (num_allocated_items1 < num_items1) {
	free(mapped1);
	num_allocated_items1 = num_items1;
	mapped1 = (size_t *) calloc(num_items1, sizeof(size_t));
      } else {
	memset(mapped1, 0, num_items1 * sizeof(size_t));
      }

      // allocated memory for items from list #2
      mapped2_counter = 0;
      if (num_allocated_items2 < num_items2) {
	free(mapped2);
	num_allocated_items2 = num_items2;
	mapped2 = (size_t *) calloc(num_items2, sizeof(size_t));
      } else {
	memset(mapped2, 0, num_items2 * sizeof(size_t));
      }

      pair_found = 0;

      for (size_t j1 = 0; j1 < num_items1; j1++) {

	if (flag1 == 1) {
	  alig = (alignment_t *) array_list_get(j1, list1);
	  chr1 = alig->chromosome;
	  strand1 = alig->seq_strand;
	  end1 = alig->position + strlen(alig->sequence);
	} else if (flag1 == 2) {
	  cal = (cal_t *) array_list_get(j1, list1);
	  chr1 = cal->chromosome_id - 1;
	  strand1 = cal->strand;
	  end1 = cal->end;
	} else {
	  end1 = -1;
	}
	
	//for (size_t j2 = num_items2 - 1; j2 > 0; j2--) {
	for (size_t j2 = 0; j2 < num_items2; j2++) {
	  if (mapped2[j2] == 1) continue;

	  if (flag2 == 1) {
	    alig = (alignment_t *) array_list_get(j2, list2);
	    chr2 = alig->chromosome;
	    strand2 = alig->seq_strand;
	    start2 = alig->position;
	  } else if (flag2 == 2) {
	    cal = (cal_t *) array_list_get(j2, list2);
	    chr2 = cal->chromosome_id - 1;
	    strand2 = cal->strand;
	    start2 = cal->start;
	  } else {
	    start2 = -1;
	  }

	  distance = (start2 > end1 ? start2 - end1 : end1 - start2); // abs
	  if ( (chr1 == chr2) &&
	       (distance >= min_distance) && (distance <= max_distance) &&
	       ((strand1 != strand2 && pair_mode == PAIRED_END_MODE) ||
		(strand1 == strand2 && pair_mode == MATE_PAIR_MODE )   ) ) {

	    
	    mapped1[j1] = 1;
	    mapped2[j2] = 1;
	    mapped1_counter++;
	    mapped2_counter++;

	    printf("***** reads %i, %i : pair1 (mapping #%i of %i), pair2 (mapping #%i of %i) : abs(end1 - start2) = (%lu - %lu) = %lu (%lu - %lu)\n", 
		   i, i+1, j1, num_items1, j2, num_items2, end1, start2, distance, min_distance, max_distance);

	    pair_found = 1;
	    break;
	  }
	}
      }

      if (pair_found) {
	void *item;
	printf("before removing: counters (mapped1, mapped2) = (%i of %i, %i of %i)\n", 
	       mapped1_counter, num_items1, mapped2_counter, num_items2);

	// removing no valid items
	//	if (mapped1_counter != num_items1) remove_items(mapped1, list1);
	//	if (mapped2_counter != num_items2) remove_items(mapped2, list2);


	printf("after removing: counters (mapped1, mapped2) = (%i of %i, %i of %i)\n", 
	       mapped1_counter, array_list_size(list1), mapped2_counter, array_list_size(list2));
	printf("\n");
      }      
    }
    //    index = batch->targets[i];
    //    printf("read to process %d\n", index);
  }

  if (mapped1 != NULL) free(mapped1);
  if (mapped2 != NULL) free(mapped2);
}

//------------------------------------------------------------------------------------

void prepare_alignments(pair_server_input_t* input, aligner_batch_t *batch) {
}

//------------------------------------------------------------------------------------

void pair_server_input_init(pair_mng_t *pair_mng, list_t* pair_list, list_t *sw_list,
			    list_t *write_list, pair_server_input_t* input) {

  input->pair_mng = pair_mng;
  input->pair_list = pair_list;
  input->sw_list = sw_list;
  input->write_list = write_list;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
