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

void pair_server_input_init(pair_mng_t *pair_mng, list_t* pair_list, list_t *sw_list,
			    list_t *write_list, pair_server_input_t* input) {

  input->pair_mng = pair_mng;
  input->pair_list = pair_list;
  input->sw_list = sw_list;
  input->write_list = write_list;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
