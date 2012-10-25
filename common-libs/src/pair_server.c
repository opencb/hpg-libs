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
*/

static void prepare_single_alignments(pair_server_input_t *input, aligner_batch_t *batch) {

  static char aux[4096];

  fastq_batch_t *fq_batch = batch->fq_batch;
  size_t index, num_items, num_targets = batch->num_targets;

  int secondary_alignment;
  size_t read_len, header_len, mapped_len, pos;
  size_t deletion_n, num_cigar_ops;
  char *cigar, *header_match, *read_match, *quality_match;

  alignment_t *alignment;

  sw_output_t *sw_output;
  array_list_t *sw_list, *alignment_list;

  // convert the SW output to alignments
  for (size_t i = 0; i < num_targets; i++) {
    //    printf("pair_server.c, prepare_single_alignments: target = #%i of %i\n", i, num_targets);
    index = batch->targets[i];
    sw_list = batch->mapping_lists[index];
    
    header_len = fq_batch->header_indices[index + 1] - fq_batch->header_indices[index] - 1;   
    num_items = array_list_size(sw_list);

    alignment_list = array_list_new(1000, 
				    1.25f, 
				    COLLECTION_MODE_ASYNCHRONIZED);

    //    printf("pair_server.c, prepare_single_alignments: process read #%i with %i mappings\n", 
    //	   index, num_items);

    for (size_t j = 0; j < num_items; j++) {
      sw_output = (sw_output_t *) array_list_get(j, sw_list);
      
      mapped_len = sw_output->mref_len;
      
      // gaps counting
      deletion_n = 0;
      for (int ii = 0; ii < mapped_len; ii++){
	if (sw_output->mquery[ii] == '-') {
	  deletion_n++;
	}
      }
      
      secondary_alignment = (j == 0 ? 1 : 0);
	    
      header_len = get_to_first_blank(&(fq_batch->header[fq_batch->header_indices[index]]), header_len, aux);
      
      header_match = (char *) malloc(sizeof(char) * (header_len + 1));
      memcpy(header_match, aux, header_len);
      header_match[header_len] = '\0';
      
      read_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
      memcpy(read_match, 
	     &fq_batch->seq[fq_batch->data_indices[index]] + sw_output->mquery_start, 
	     mapped_len - deletion_n);
      read_match[mapped_len - deletion_n ] = '\0';

      quality_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
      memcpy(quality_match, 
	     &fq_batch->quality[fq_batch->data_indices[index]] + sw_output->mquery_start, 
	     mapped_len - deletion_n);
      quality_match[mapped_len - deletion_n] = '\0';
      
      cigar =  generate_cigar_str(sw_output->mquery, 
				  sw_output->mref, 
				  sw_output->mquery_start, 
				  read_len, 
				  sw_output->mref_len, 
				  &num_cigar_ops);
      
      pos = sw_output->ref_start + sw_output->mref_start;
      /*
      if (sw_output->strand == 0) {
	pos = sw_output->ref_start + sw_output->mref_start - sw_output->strand;
      } else {
	pos = sw_output->ref_start + (sw_output->ref_len - 
				      sw_output->mref_len - 
				      sw_output->mref_start);
      }
      */
      // create the alignment and insert into the list
      alignment = alignment_new();
      alignment_init_single_end(header_match, read_match, quality_match, 
				sw_output->strand, 
				sw_output->chromosome - 1, 
				pos - 1,
				cigar, num_cigar_ops, 255, secondary_alignment, 1, alignment);
      array_list_insert(alignment, alignment_list);

      // free memory (sw output)
      sw_output_free(sw_output);
    } // end for sw items

    // free the sw list, and update the mapping list with the alignment list
    array_list_free(sw_list, NULL);
    batch->mapping_lists[index] = alignment_list;
  } // end for targets
  //  printf("pair_server.c, prepare_single_alignments: Done\n");
}

//------------------------------------------------------------------------------------
inline void update_mispaired_pair(int pair_num, size_t num_items, array_list_t *list) {
  alignment_t *alig;

  for (size_t i = 0; i < num_items; i++) {
    alig = (alignment_t *) array_list_get(i, list);

    // set pair fields
    alig->mate_position = 0;
    alig->mate_chromosome = 0;
    alig->template_length = 0;
     
    alig->is_paired_end = 1;
    alig->is_paired_end_mapped = 0;
    alig->is_mate_mapped = 0;
    alig->mate_strand = 0;
    alig->pair_num = pair_num;
  }
}

//------------------------------------------------------------------------------------

inline void update_mispaired_pairs(size_t num_items1, size_t num_items2,
				   array_list_t *list1, array_list_t *list2) {
  alignment_t *alig;
  alignment_t *first1 = array_list_get(0, list1);
  alignment_t *first2 = array_list_get(0, list2);

  for (size_t i = 0; i < num_items1; i++) {
    alig = (alignment_t *) array_list_get(i, list1);

    // set pair1 fields
    alig->mate_position = first2->position;
    alig->mate_chromosome = first2->chromosome;
    alig->template_length = first2->position - alig->position;;
     
    alig->is_paired_end = 1;
    alig->is_paired_end_mapped = 0;
    alig->is_mate_mapped = 1;
    alig->mate_strand = first2->seq_strand;
    alig->pair_num = 1;
  }

  for (size_t i = 0; i < num_items2; i++) {
    alig = (alignment_t *) array_list_get(i, list2);

    // set pair1 fields
    alig->mate_position = first1->position;
    alig->mate_chromosome = first1->chromosome;
    alig->template_length = first1->position - alig->position;;
     
    alig->is_paired_end = 1;
    alig->is_paired_end_mapped = 0;
    alig->is_mate_mapped = 1;
    alig->mate_strand = first1->seq_strand;
    alig->pair_num = 2;
  }
}

//------------------------------------------------------------------------------------

void prepare_paired_alignments(pair_server_input_t *input, aligner_batch_t *batch) {

  fastq_batch_t *fq_batch = batch->fq_batch;
  size_t num_items1, num_items2, num_reads = fq_batch->num_reads;

  int distance;
  int min_distance = input->pair_mng->min_distance;
  int max_distance = input->pair_mng->max_distance;
  int pair_mode = input->pair_mng->pair_mode;

  array_list_t *list1, *list2;

  alignment_t *alig1, *alig2;
  size_t mapped1_counter = 0, mapped2_counter = 0;
  size_t allocated_mapped1 = 100, allocated_mapped2 = 100;
  size_t *mapped1 = (size_t *) malloc(allocated_mapped1 * sizeof(size_t));
  size_t *mapped2 = (size_t *) malloc(allocated_mapped2 * sizeof(size_t));

  short int chr1, chr2, strand1, strand2;
  size_t end1, start2;

  int pair_found;

  for (int i = 0; i < num_reads; i += 2) {

    list1 = batch->mapping_lists[i];
    list2 = batch->mapping_lists[i+1];

    num_items1 = 0;
    if (list1 != NULL)  num_items1 = array_list_size(list1);
    num_items2 = 0;
    if (list2 != NULL) num_items2 = array_list_size(list2);

    printf("prepare_paired_alignments:  reads %i, %i : pair1 (%i mappings, allocated = %i), pair2 (%i mappins, allocated = %i)\n",
	   i, i+1, num_items1, allocated_mapped1, num_items2, allocated_mapped2);
    
    if (num_items1 > 0 && num_items2 > 0) {

      // initalizes memory and counters
      if (allocated_mapped1 < num_items1) {
	free(mapped1);
	mapped1 = (size_t *) malloc(num_items1 * sizeof(size_t));
	allocated_mapped1 = num_items1;
      }
      memset(mapped1, 0, num_items1 * sizeof(size_t));
      mapped1_counter = 0;

      if (allocated_mapped2 < num_items2) {
	free(mapped2);
	mapped2 = (size_t *) malloc(num_items2 * sizeof(size_t));
	allocated_mapped2 = num_items1;
      }
      memset(mapped2, 0, num_items2 * sizeof(size_t));
      mapped2_counter = 0;

      pair_found = 0;

      // search for pairs properly aligned
      for (size_t j1 = 0; j1 < num_items1; j1++) {
	alig1 = (alignment_t *) array_list_get(j1, list1);
	chr1 = alig1->chromosome;
	strand1 = alig1->seq_strand;
	end1 = alig1->position;
	
	for (size_t j2 = 0; j2 < num_items2; j2++) {
	  if (mapped2[j2] == 1) continue;
	  alig2 = (alignment_t *) array_list_get(j2, list2);
	  chr2 = alig2->chromosome;
	  strand2 = alig2->seq_strand;
	  start2 = alig2->position;
	  
	  // computes distance between alignments,
	  // is a valid distance ?
	  distance = (start2 > end1 ? start2 - end1 : end1 - start2); // abs                                      
	  if ( (chr1 == chr2) &&
	       (distance >= min_distance) && (distance <= max_distance) &&
	       ((strand1 != strand2 && pair_mode == PAIRED_END_MODE) ||
		(strand1 == strand2 && pair_mode == MATE_PAIR_MODE )   ) ) {
	    
            mapped1[j1] = 1;
	    mapped2[j2] = 1;
	    
	    mapped1_counter++;
	    mapped2_counter++;
	    
	    // set pair1 fields
	    alig1->mate_position = alig2->position;
	    alig1->mate_chromosome = alig2->chromosome;
	    alig1->template_length = alig2->position - alig1->position;
     
	    alig1->is_paired_end = 1;
	    alig1->is_paired_end_mapped = 1;
	    alig1->is_mate_mapped = 1;
	    alig1->mate_strand = alig2->seq_strand;
	    alig1->pair_num = 1;

	    // set pair2 fields
	    alig2->mate_position = alig1->position;
	    alig2->mate_chromosome = alig1->chromosome;
	    alig2->template_length = alig1->position - alig2->position;
	    
	    alig2->is_paired_end = 1;
	    alig2->is_paired_end_mapped = 1;
	    alig2->is_mate_mapped = 1;
	    alig2->mate_strand = alig1->seq_strand;
	    alig2->pair_num = 2;
	    
	    printf("***** reads %i, %i : pair1 (mapping #%i of %i), pair2 (mapping #%i of %i) : abs(end1 - start2)\n",
                   i, i+1, j1, num_items1, j2, num_items2, end1, start2, distance, min_distance, max_distance);
	    
	    pair_found = 1;
            break;
	  }
	}
      }

      // check if there are unproperly aligned pairs
      if (pair_found) {
	// remove unpaired alignments
	if (mapped1_counter != num_items1) remove_alignments(mapped1, list1); 
	if (mapped2_counter != num_items2) remove_alignments(mapped2, list2); 
      } else {
	// all aligments are unpaired
	update_mispaired_pairs(num_items1, num_items2, list1, list2);
      }
    } else {
      // pairs are not properly aligned, only one is mapped
      update_mispaired_pair(1, num_items1, list1);
      update_mispaired_pair(2, num_items2, list2);
    }
  } // end for num_reads

  // free memory
  free(mapped1);
  free(mapped2);
}

//------------------------------------------------------------------------------------

inline int remove_alignments(size_t *valid_items, array_list_t *list) {
  alignment_t *alig;
  int num = 0;

  size_t num_items = array_list_size(list);

  for (int k = num_items - 1; k >= 0; k--) {
    if (valid_items[k] == 0) {
      alig = (alignment_t *) array_list_remove_at(k, list);
      alignment_free(alig);
      ++num;
    }
  }
  return num;
}

//------------------------------------------------------------------------------------

inline int remove_items(size_t *valid_items, array_list_t *list) {
  void *item;
  int num = 0;

  size_t num_items = array_list_size(list);
  int flag = array_list_get_flag(list);

  printf("list of %i items\n", num_items);
  for (int k = num_items - 1; k >= 0; k--) {
    if (valid_items[k] == 0) {
      printf("%i, ", k);
      item = (void *) array_list_remove_at(k, list);
      ++num;
      if (flag == 1) {
	alignment_free((alignment_t *) item);
      } else {
	cal_free((cal_t *) item);
      }
    }
  }
  printf("\n");
  return num;
}

//====================================================================================
// main functions: apply pair and prperare alignments
//====================================================================================

void apply_pair(pair_server_input_t* input, aligner_batch_t *batch) {

  {
    size_t index, num_seqs = batch->num_targets;
    for (size_t i = 0; i < num_seqs; i++) {
      index = batch->targets[i];
      printf("apply_pair: read #%i of %i: with %i cals\n", index, num_seqs, 
	     array_list_size(batch->mapping_lists[index]));
    }
  }
    
  //  printf("START: apply_pair\n"); 
  char *seq;
  list_t *list = NULL;
  fastq_batch_t *fq_batch = batch->fq_batch;

  int pair_mode = input->pair_mng->pair_mode;
  size_t min_distance = input->pair_mng->min_distance;
  size_t max_distance = input->pair_mng->max_distance;
  int distance;

  size_t num_targets = batch->num_targets;
  size_t num_items1, num_items2, num_reads = fq_batch->num_reads;


  int flag1, flag2;
  array_list_t *list1, *list2;

  size_t end1, start2;
  short int chr1, chr2, strand1, strand2;

  size_t mapped1_counter = 0, mapped2_counter = 0;
  size_t allocated_mapped1 = 100, allocated_mapped2 = 100;
  size_t *mapped1 = (size_t *) malloc(allocated_mapped1 * sizeof(size_t));
  size_t *mapped2 = (size_t *) malloc(allocated_mapped2 * sizeof(size_t));

  int pair_found;

  alignment_t *alig;
  cal_t *cal;
  
  int total_removed = 0;
  printf("pair_server.c:apply_pair: pair_mode = %i, min_distance = %lu, max_distance = %lu\n",
	 pair_mode, min_distance, max_distance);

  for (size_t i = 0; i < num_reads; i += 2) {

    list1 = batch->mapping_lists[i];
    list2 = batch->mapping_lists[i + 1];

    flag1 = array_list_get_flag(list1);
    flag2 = array_list_get_flag(list2);

    num_items1 = 0;
    if (list1 != NULL)  num_items1 = array_list_size(list1);
    num_items2 = 0;
    if (list2 != NULL) num_items2 = array_list_size(list2);

    if (num_items1 > 1 && num_items2 > 1) {

      // initalizes memory and counters
      mapped1_counter = 0;
      if (allocated_mapped1 < num_items1) {
	free(mapped1);
	mapped1 = (size_t *) malloc(num_items1 * sizeof(size_t));
	allocated_mapped1 = num_items1;
      }
      memset(mapped1, 0, num_items1 * sizeof(size_t));
      mapped1_counter = 0;

      mapped2_counter = 0;
      if (allocated_mapped2 < num_items2) {
	free(mapped2);
	mapped2 = (size_t *) malloc(num_items2 * sizeof(size_t));
	allocated_mapped2 = num_items1;
      }
      memset(mapped2, 0, num_items2 * sizeof(size_t));
      mapped2_counter = 0;

      pair_found = 0;

      // search for pairs properly aligned
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
	  printf("Error in pair_server.c, apply_sw function\n");
	  abort();
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

	  // computes distance between alignments,
	  // is a valid distance ?
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
	} // end for j2..num_items2
      } // end for j1..num_item1


      if (pair_found) {
	printf("before removing: counters (mapped1, mapped2) = (%i of %i, %i of %i)\n", 
	       mapped1_counter, num_items1, mapped2_counter, num_items2);

	// removing no valid items
	if (mapped1_counter != num_items1) {
	  remove_items(mapped1, list1);
	  if (flag1 == 2) total_removed += (num_items1 - mapped1_counter);
	}
	if (mapped2_counter != num_items2) {
	  remove_items(mapped2, list2);
	  if (flag1 == 2) total_removed += (num_items2 - mapped2_counter);
	}


	printf("after removing: counters (mapped1, mapped2) = (%i of %i, %i of %i)\n", 
	       mapped1_counter, array_list_size(list1), mapped2_counter, array_list_size(list2));
	printf("\n");
      }      
    }
  }

  printf("sw to do = %i, total removed = %i\n", batch->num_to_do, total_removed);
  batch->num_to_do -= total_removed;

  // free memory
  free(mapped1);
  free(mapped2);

  {
    size_t index, num_seqs = batch->num_targets;
    for (size_t i = 0; i < num_seqs; i++) {
      index = batch->targets[i];
      printf("apply_pair: read #%i of %i: with %i cals\n", index, num_seqs, 
	     array_list_size(batch->mapping_lists[index]));
    }
  }

}

//------------------------------------------------------------------------------------

void prepare_alignments(pair_server_input_t *input, aligner_batch_t *batch) {
  prepare_single_alignments(input, batch);
  //printf("pair_server.c: 0: after prepare_single_alignments\n");
  /*
  if (input->pair_mng->pair_mode != SINGLE_END_MODE) {
    //prepare_paired_alignments(input, batch);
  }
  */
  //printf("pair_server.c: prepare_alignments done (pair mode = %i)\n", input->pair_mng->pair_mode);
  //  printf("pair_server.c: 1: after prepare_single_alignments\n");
}

//------------------------------------------------------------------------------------
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
