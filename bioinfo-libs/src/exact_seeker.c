
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "commons.h"
#include "list.h"
#include "system_utils.h"

#include "BW_io.h"

#include "timing.h"

#include "fastq_file.h"
#include "fastq_ex_batch.h"
#include "exact_seeker.h"
#include "string_utils.h"

//------------------------------------------------------------------------------------

seed_batch_t* seed_read(char* read_p, char* quality_p, unsigned int read_len, char* header_p, unsigned int header_len, list_t* list_p, seed_batch_t* current_batch_p);

//------------------------------------------------------------------------------------
// this functions search for read matching from 
// k and l vector store in a list
//------------------------------------------------------------------------------------

void exact_seeker(exact_seeker_input_t* input_p) {

  char header[1024], search[4096], quality[4096], header_id[1024];
  char *header_id_match_p, *search_match_p, *quality_match_p, *cigar_p;
  char plusminus[] = "+-";
  unsigned int num_gpus = 2;
  unsigned int k_aux, l_aux;
  unsigned int i, j, k, z, bytes, num_reads, found = 0;
  unsigned int len, len_header, index, key, total_invalids = 0, invalid = 0;
  unsigned int len_id_header;
  short int gpu;
  short int chromosome;
  unsigned short primary_alignment;
  printf("exact_seeker (%i): START\n", omp_get_thread_num());

  struct timespec ts;
  ts.tv_sec = 1;
  ts.tv_nsec = 0;

  unsigned int seed_batch_size = 1000000;
  unsigned int seed_size = input_p->seed_size;
  unsigned int num_seeds_per_read = input_p->num_seeds;

  // lists and items
  //
  list_t* kl_list_p = input_p->kl_list_p;
  list_item_t *kl_item_p = NULL;
  kl_t* kl_p = NULL;

  list_item_t* item_p = NULL;

  list_t* seed_list_p = input_p->seed_list_p;
  fastq_ex_batch_t* fastq_ex_batch_p = NULL;

  list_t* write_list_p = input_p->write_list_p;
  unsigned int write_size = input_p->write_size;
  write_batch_t* write_batch_p = write_batch_new(write_size, 1);
  
  seed_batch_t *seed_p, *seed_batch_p = seed_batch_new(seed_batch_size, seed_size, num_seeds_per_read);

  bwt_index_t* index_p = input_p->index_p;
  genome_t* genome_p = input_p->genome_p;

  int no_exacts = 0, total = 0;
  alignment_t *alignment_p;
  alignment_t **buffer_p;
  
  while ( (kl_item_p = list_remove_item(kl_list_p)) != NULL ) {
    //printf("EXACT SEEKER EXTRACT ONE ELEMENT kl_list_p\n");
    if (time_on) { timing_start(EXACT_SEEKER_INDEX, 0, timing_p); }

    //    if (time_on) { start_timer(t1_exact_seeker); }
    kl_p = (kl_t*) kl_item_p->data_p;
    fastq_ex_batch_p = (fastq_ex_batch_t*) kl_p->batch_p;
    num_reads = fastq_ex_batch_p->num_reads;
    total += num_reads;

    for (i=0; i < num_reads; i++) {
      invalid = 0;
      found=0;
      
      len = fastq_ex_batch_p->data_indices[i+1] - fastq_ex_batch_p->data_indices[i] - 1;
      decodeBases(search, &(fastq_ex_batch_p->seq[fastq_ex_batch_p->data_indices[i]]), len);
      memcpy(quality, &(fastq_ex_batch_p->quality[fastq_ex_batch_p->data_indices[i]]), len);
      
      len_header = fastq_ex_batch_p->header_indices[i+1] - fastq_ex_batch_p->header_indices[i] - 1;
      memcpy(header, &(fastq_ex_batch_p->header[fastq_ex_batch_p->header_indices[i]]), len_header);
      
      for(gpu=0 ; gpu<num_gpus ; gpu++) {
	
	k_aux = kl_p->k_p[i + (num_reads * gpu)];
	l_aux = kl_p->l_p[i + (num_reads * gpu)];

	if (l_aux - k_aux + 1 < 100) {
	  
	  //header_match_p = (char *)malloc(sizeof(char)*len_header);
	  len_id_header = get_to_first_blank(&header, len_header, &header_id);
	  
	  //printf("ID READ::%s\n", header_match_p);
	  primary_alignment = 0;
	  for (j=k_aux; j<=l_aux; j++) {
	    key = index_p->S[j];
	    index = binsearch(index_p->offset_p, index_p->size, key);
	    //printf("exact_seeker.c: binarysearch: index = %i, k=%d, strand=%i\n", index, j, gpu);
	    
	    if(key + len <= index_p->offset_p[index]) {
	      found=1;
	      // found-reads handling,
	      //
	      //printf("found: k=%d,l=%d, %s\t%c\t%s %d %s\n", k_aux, l_aux, header, plusminus[gpu], index_p->chromosome_p + (index-1)*index_p->max_id, index_p->start_p[index-1] + (key - index_p->offset_p[index-1]), search);
	      
	      if ( write_batch_p->size > write_batch_p->allocated_size ) {
		item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
		if (time_on) { timing_stop(EXACT_SEEKER_INDEX, 0, timing_p); }
		list_insert_item(item_p, write_list_p);
		if (time_on) { timing_start(EXACT_SEEKER_INDEX, 0, timing_p); }
		
		write_batch_p = write_batch_new(write_size, MATCH_FLAG);
		
	      }
	      
	      search_match_p = (char *)malloc(sizeof(char)*(len + 1));
	      memcpy(search_match_p, &search, len);
	      search_match_p[len] = '\0';
	      
	      quality_match_p = (char *)malloc(sizeof(char)*(len + 1));
	      memcpy(quality_match_p, &quality, len);
	      quality_match_p[len] = '\0';
	  
	      header_id_match_p = (char *)malloc(sizeof(char)*len_id_header);
	      memcpy(header_id_match_p, &header_id, len_id_header);
	      
	      alignment_p = alignment_new();
	      cigar_p = (char *)malloc(sizeof(char)*10);
	      sprintf(cigar_p, "%d%c\0", len, '=');
	      
	      chromosome = atoi(genome_p->chr_name[index - 1]);
	     
	      alignment_init_single_end(header_id_match_p, search_match_p, quality_match_p, 0, chromosome - 1, index_p->start_p[index-1] + (key - index_p->offset_p[index-1]), cigar_p, 1, 255, 0, 1, alignment_p);
	      //printf("seq: %s\n", alignment_p->sequence);
	      ((alignment_t **)write_batch_p->buffer_p)[write_batch_p->size] = alignment_p;
	      write_batch_p->size++;
	      
	      //break; // uncomment this if you only want the first alignment
	    }
	  } // end of for k_aux..l_aux
	} else {
	  invalid = 1;
	}// end of if mappings < 100
      } // end of for 0..num_gpus
      
      // not-found-reads handling, i.e., seed the read and insert in the corresponding list
      // to be processed by gpu_server
      //
      if (!found) {
	if (invalid) total_invalids++;
	no_exacts++;
	seed_p = seed_read(&(fastq_ex_batch_p->seq[fastq_ex_batch_p->data_indices[i]]), 
			     &(fastq_ex_batch_p->quality[fastq_ex_batch_p->data_indices[i]]), 
			     len,
			     header, len_header,
			     seed_list_p, seed_batch_p);
	
	if (seed_p == NULL) {
	  // seed not possible, write to disk
	  //
	  /*if ( (2*(len_header + len) + notfound_write_p->size) > write_size) {
	    item_p = list_item_new(0, WRITE_ITEM, notfound_write_p);
	    if (time_on) { timing_stop(EXACT_SEEKER_INDEX, 0, timing_p); }
	    list_insert_item(item_p, write_list_p);
	    if (time_on) { timing_start(EXACT_SEEKER_INDEX, 0, timing_p); }
	    
	    notfound_write_p = write_batch_new(write_size, MISMATCH_FLAG);
	  }
	  memcpy(quality, &(fastq_ex_batch_p->quality[fastq_ex_batch_p->data_indices[i]]), len);

	  bytes = packFastQ(header, search, quality, len_header, len, &(((char *)notfound_write_p->buffer_p)[notfound_write_p->size]));
	  notfound_write_p->size += bytes;
	  if (notfound_write_p->size > notfound_write_p->allocated_size) {
	    printf("ERROR: exceed, %i over %i\n", notfound_write_p->size, notfound_write_p->allocated_size);
	  }*/
	  
	  if ( write_batch_p->size > write_batch_p->allocated_size ) {
	    item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
	    if (time_on) { timing_stop(EXACT_SEEKER_INDEX, 0, timing_p); }
	    list_insert_item(item_p, write_list_p);
	    if (time_on) { timing_start(EXACT_SEEKER_INDEX, 0, timing_p); }
	    
	    write_batch_p = write_batch_new(write_size, MATCH_FLAG);
	  }
	  
	  search_match_p = (char *)malloc(sizeof(char)*(len + 1));
	  memcpy(search_match_p, &search, len);
	  search_match_p[len] = '\0';
	  
	  quality_match_p = (char *)malloc(sizeof(char)*(len + 1));
	  memcpy(quality_match_p, &quality, len);
	  quality_match_p[len] = '\0';
      
	  len_id_header = get_to_first_blank(&header, len_header, &header_id);
	  
	  header_id_match_p = (char *)malloc(sizeof(char)*len_id_header);
	  memcpy(header_id_match_p, &header_id, len_id_header);
	  
	  alignment_p = alignment_new();
	  cigar_p = (char *)malloc(sizeof(char)*10);
	  sprintf(cigar_p, "%d%c\0", len, 'X');
	  //TODO:chromosome 0??
	  chromosome = 0;
	  
	  alignment_init_single_end(header_id_match_p, search_match_p, quality_match_p, gpu, chromosome - 1, 0, cigar_p, 1, 255, primary_alignment, 0, alignment_p);
	  primary_alignment = 1;
	  //printf("seq: %s\n", alignment_p->sequence);
	  ((alignment_t **)write_batch_p->buffer_p)[write_batch_p->size] = alignment_p;
	  write_batch_p->size++;
	} else {
	  seed_batch_p = seed_p;
	}
      } // end of !found
    } // end of for num_reads
      
    //printf("num of bytes in not found = %i, allocated = %i\n", notfound_write_p->size, notfound_write_p->allocated_size);
    
    fastq_ex_batch_free(fastq_ex_batch_p);
    kl_free(kl_p);

    //    if (time_on) { stop_timer(t1_exact_seeker, t2_exact_seeker, exact_seeker_time); }
    //printf("exact_seeker: processing kl item (%x) done !!\n", kl_item_p);

    list_item_free(kl_item_p);
    
    //printf(" EXACT SEEKER FINISH Process ONE Batch\n");
    if (time_on) { timing_stop(EXACT_SEEKER_INDEX, 0, timing_p); }
  } // end of while 

  // insert or free memory
  //
  if (seed_batch_p != NULL) {
    if (seed_batch_p->num_reads > 0) {
      item_p = list_item_new(0, SEED_ITEM, seed_batch_p);
      list_insert_item(item_p, seed_list_p);
    } else {
      seed_batch_free(seed_batch_p);
    }
  }  

  if (write_batch_p != NULL) {
    if (write_batch_p->size > 0) {
      item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
      list_insert_item(item_p, write_list_p);
    } else {
      write_batch_free(write_batch_p);
    }
  }

 /* if (notfound_write_p != NULL) {
    if (notfound_write_p->size > 0) {
      item_p = list_item_new(0, WRITE_ITEM, notfound_write_p);
      list_insert_item(item_p, write_list_p);
    } else {
      write_batch_free(notfound_write_p);
    }
  }*/

  list_decr_writers(seed_list_p);
  list_decr_writers(write_list_p);

  printf("exact_seeker: END (total reads = %i, exacts = %i, no-exacts = %i, invalids = %i)\n", total, total - no_exacts, no_exacts, total_invalids);
}

//------------------------------------------------------------------------------------

seed_batch_t* seed_read(char* read_p, char* quality_p, unsigned int read_len, char* header_p, unsigned int header_len, list_t* list_p, seed_batch_t* current_batch_p) {
  seed_batch_t* batch_p = current_batch_p;

  unsigned int num_seeds;


  if ( (read_len <= batch_p->seed_size) || ( (num_seeds = read_len - batch_p->seed_size) < batch_p->num_seeds_per_read) ) {
   printf("ERROR: impossible to seed (%i, min. %i) read\n", num_seeds, batch_p->num_seeds_per_read);
   return NULL;
  }

  if (num_seeds > batch_p->num_seeds_per_read) {
    num_seeds = batch_p->num_seeds_per_read;
  }

  if ( (batch_p->num_reads + 1 >= batch_p->allocated_read_indices) ||
       (batch_p->num_seeds + num_seeds > batch_p->allocated_seed_indices) ||
       (batch_p->header_size + header_len > batch_p->allocated_header_size) ||
       (batch_p->read_size + read_len > batch_p->allocated_read_size) ) {
       
    //printf("--------- not space enough !!, allocate new batch\n");

    list_item_t* item_p = list_item_new(0, SEED_ITEM, batch_p);
    if (time_on) { timing_stop(EXACT_SEEKER_INDEX, 0, timing_p); }
    list_insert_item(item_p, list_p);
    if (time_on) { timing_start(EXACT_SEEKER_INDEX, 0, timing_p); }
    
    batch_p = seed_batch_new(batch_p->max_size, batch_p->seed_size, batch_p->num_seeds_per_read);
  }
  
  batch_p->read_indices_p[batch_p->num_reads] = batch_p->read_size;
  memcpy(&batch_p->read_p[batch_p->read_size], read_p, read_len);
  memcpy(&batch_p->quality_p[batch_p->read_size], quality_p, read_len);

  batch_p->header_indices_p[batch_p->num_reads] = batch_p->header_size;
  memcpy(&batch_p->header_p[batch_p->header_size], header_p, header_len);

  int i;
  double shift = (read_len - batch_p->seed_size) / num_seeds;
  for(i=0 ; i<num_seeds ; i++) {
    //printf("mini read %i index: %i\n", batch_p->num_seeds + i, batch_p->read_size + (unsigned int) (shift*i));    
    batch_p->seed_indices_p[batch_p->num_seeds + i] = batch_p->read_size + (unsigned int) (shift*i);    
  }

  batch_p->num_reads++;
  batch_p->num_seeds += num_seeds;
  batch_p->header_size += header_len;
  batch_p->read_size += read_len;

  batch_p->read_indices_p[batch_p->num_reads] = batch_p->read_size;
  batch_p->header_indices_p[batch_p->num_reads] = batch_p->header_size;

  return batch_p;
}

//------------------------------------------------------------------------------------
// exact_seeker_input functions: init
//------------------------------------------------------------------------------------

void exact_seeker_input_init(list_t* kl_list_p, list_t* seed_list_p, list_t* write_list_p, unsigned int num_seeds, unsigned int seed_size, unsigned int write_size, bwt_index_t* bwt_index_p, genome_t* genome_p, exact_seeker_input_t* input_p) {
  input_p->kl_list_p = kl_list_p;
  input_p->seed_list_p = seed_list_p;
  input_p->write_list_p = write_list_p;

  input_p->num_seeds = num_seeds;
  input_p->seed_size = seed_size;
  input_p->write_size = write_size;

  input_p->index_p = bwt_index_p;
  input_p->genome_p = genome_p;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
