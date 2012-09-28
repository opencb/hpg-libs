#include "batch_writer.h"

//------------------------------------------------------------------------------------
// this functions writes reads to disk, reads came from
// a list
//------------------------------------------------------------------------------------

void batch_writer(batch_writer_input_t* input_p) {

  struct timespec ts;
  ts.tv_sec = 1;
  ts.tv_nsec = 0;

  alignment_t **buffer_p;
  bam1_t* bam1_p;
  bam_header_t* bam_header_p;
  bam_file_t* bam_file_p;
  
  char* match_filename = input_p->match_filename;
  //char* mismatch_filename = input_p->mismatch_filename;

  char* splice_exact_filename = input_p->splice_exact_filename;
  char* splice_extend_filename = input_p->splice_extend_filename;

  list_t* list_p = input_p->list_p;

  printf("batch_writer (%i): START\n", omp_get_thread_num());
		
  list_item_t *item_p = NULL;
  write_batch_t* batch_p;

  FILE* fd;
  FILE* splice_exact_fd  = fopen(splice_exact_filename, "w");
  FILE* splice_extend_fd = fopen(splice_extend_filename, "w");

  bam_header_p = bam_header_new(HUMAN, NCBI37);
  //bam_file_p = bam_fopen(match_filename);
  bam_file_p = bam_fopen_mode(match_filename, bam_header_p, "w");
  bam_fwrite_header(bam_header_p, bam_file_p);

  // main loop
  while ( (item_p = list_remove_item(list_p)) != NULL ) {
    
    if (time_on) { timing_start(BATCH_WRITER, 0, timing_p); }

    batch_p = (write_batch_t*) item_p->data_p;
    //printf("*********************************Extract one item*********************************\n");
    if (batch_p->flag == MATCH_FLAG || batch_p->flag == MISMATCH_FLAG) { //fd = match_fd; 
      //printf("start write alignment. Total %d\n", batch_p->size);
      buffer_p = (alignment_t **)batch_p->buffer_p;
      for(int i = 0; i < batch_p->size; i++){
	//alignment_print(buffer_p[i]);
	bam1_p = convert_to_bam(buffer_p[i], 33);
	bam_fwrite(bam1_p, bam_file_p);
	bam_destroy1(bam1_p);
	alignment_free(buffer_p[i]);
      }
    }else{
      if (batch_p->flag == SPLICE_EXACT_FLAG){ fd = splice_exact_fd; }
      else if (batch_p->flag == SPLICE_EXTEND_FLAG){ fd = splice_extend_fd; }
      else { fd = NULL; }
      
      if (fd != NULL){
	//printf("start write batch, %i bytes...\n", batch_p->size);
	fwrite((char *)batch_p->buffer_p, batch_p->size, 1, fd);
	//printf("write done !!\n");
	//if (time_on) { stop_timer(t1_write, t2_write, write_time); }
      }
    }
    //printf("Free batch\n");
    write_batch_free(batch_p);
    list_item_free(item_p);

    if (time_on) { timing_stop(BATCH_WRITER, 0, timing_p); }
  } // end of batch loop
  
  //fclose(match_fd);
  //fclose(mismatch_fd);
  fclose(splice_exact_fd);
  fclose(splice_extend_fd);

  bam_fclose(bam_file_p);
  //bam_header_free(bam_header_p);
  printf("batch_writer: END\n");
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
//------------------------------------------------------------------------------------

void batch_writer2(batch_writer_input_t* input) {

  //  printf("START: batch_writer (%i)\n", omp_get_thread_num());

  bam1_t *bam1;
  bam_header_t *bam_header;
  bam_file_t *bam_file;
  alignment_t *alig, *alignment = alignment_new();

  char* match_filename = input->match_filename;
  //  char* splice_filename = input->splice_filename;
  
  list_t *write_list = input->list_p;
  array_list_t *array_list;
  
  list_item_t *item = NULL;
  aligner_batch_t *batch = NULL;
  fastq_batch_t *fq_batch = NULL;

  FILE* fd;
  //  FILE* splice_fd = fopen(splice_filename, "w");

  char aux[2048], *seq;
  
  sw_output_t *sw_output;
  size_t read_len, header_len, mapped_len, pos;
  size_t deletion_n, primary_alignment, num_cigar_ops;
  char *cigar, *header_match, *read_match, *quality_match;

  bam_header = bam_header_new(HUMAN, NCBI37);
  bam_file = bam_fopen_mode(match_filename, bam_header, "w");
  bam_fwrite_header(bam_header, bam_file);

  int flag = 0;
  size_t num_reads = 0, num_items = 0, total_mappings = 0;

  // main loop
  while ( (item = list_remove_item(write_list)) != NULL ) {

    batch = (aligner_batch_t *) item->data_p;
    fq_batch = batch->fq_batch;
    num_reads = batch->num_mapping_lists;

    // single mode (single-end)
    for (size_t i = 0; i < num_reads; i++) {

      array_list = batch->mapping_lists[i];
      //if (array_list == NULL) printf("READ %d, writer, list is NULL\n", i);

      //      printf("----> list == NULL ? %d\n", (array_list == NULL));
      num_items = (array_list == NULL ? 0 : array_list_size(array_list));
      //printf("----> number of items = %d, num_items <= 0 ? %d\n", num_items, num_items <= 0);

      read_len = fq_batch->data_indices[i + 1] - fq_batch->data_indices[i] - 1;
      if (num_items == 0) {

	// unmapped read
	//printf("\tWRITE : read %d (%d items): unmapped!!\n", i, num_items);

	// calculating cigar
	sprintf(aux, "%dX", read_len);
	
	alignment_init_single_end(&(fq_batch->header[fq_batch->header_indices[i]])+1,
				  &(fq_batch->seq[fq_batch->data_indices[i]]),
				  &(fq_batch->quality[fq_batch->data_indices[i]]),
				  0, 
				  0,
				  0, 
				  aux, 1, 255, 0, 0, alignment);

	bam1 = convert_to_bam(alignment, 33);
	bam_fwrite(bam1, bam_file);
	bam_destroy1(bam1);

      } else {

	// mapped read by BWT or SW ??
	flag = array_list_get_flag(array_list);

	if (flag == 1) {
	  // mapped by BWT
	  // here we have alignment_t objects 
	  //printf("WRITE : read %d (%d items): mapped by BWT\n", i, num_items);

	  for (size_t j = 0; j < num_items; j++) {
	    alig = (alignment_t *) array_list_get(j, array_list);
	    if (alig != NULL) {

	      bam1 = convert_to_bam(alig, 33);
	      bam_fwrite(bam1, bam_file);
	      bam_destroy1(bam1);

	      alignment_free(alig);
	    }
	  }
	} else {
	  // mapped by SW
	  // here we have sw_output_t objects
	  //printf("WRITE : read %d (%d items): mapped by SW\n", i, num_items);

	  header_len = fq_batch->header_indices[i + 1] - fq_batch->header_indices[i] - 1;

	  for (size_t j = 0; j < num_items; j++) {
	    sw_output = (sw_output_t *) array_list_get(j, array_list);

	    mapped_len = sw_output->mref_len;

	    deletion_n = 0;
	    for (int ii = 0; ii < mapped_len; ii++){
	      if (sw_output->mquery[ii] == '-') {
		deletion_n++;
	      }
	    }
	    
	    if (j == 0) {
	      primary_alignment = 1;
	    } else {
	      primary_alignment = 0;
	    }
	    
	    header_len = get_to_first_blank(&(fq_batch->header[fq_batch->header_indices[i]]), header_len, aux);
	    //	    aux = &(fq_batch->header_indices[i]);

	    header_match = (char *) malloc(sizeof(char) * (header_len + 1));
	    memcpy(header_match, aux, header_len);
	    header_match[header_len] = '\0';

	    read_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
	    memcpy(read_match, 
		   &fq_batch->seq[fq_batch->data_indices[i]] + sw_output->mquery_start, 
		   mapped_len - deletion_n);
	    read_match[mapped_len - deletion_n ] = '\0';

	    quality_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
	    memcpy(quality_match, 
		   &fq_batch->quality[fq_batch->data_indices[i]] + sw_output->mquery_start, 
		   mapped_len - deletion_n);
	    quality_match[mapped_len - deletion_n] = '\0';

	    cigar =  generate_cigar_str(sw_output->mquery, 
					sw_output->mref, 
					sw_output->mquery_start, 
					read_len, 
					sw_output->mref_len, 
					&num_cigar_ops);


	    if (sw_output->strand == 0) {
	      pos = sw_output->ref_start + sw_output->mref_start - sw_output->strand;
	    } else {
	      pos = sw_output->ref_start + (sw_output->ref_len - 
					   sw_output->mref_len - 
					   sw_output->mref_start);
	    }
	    
	    alignment_init_single_end(header_match, read_match, quality_match, 
				      sw_output->strand, 
				      sw_output->chromosome - 1, 
				      pos - 1,
				      cigar, num_cigar_ops, 255, primary_alignment, 1, alignment);

	    bam1 = convert_to_bam(alignment, 33);
	    bam_fwrite(bam1, bam_file);
	    bam_destroy1(bam1);

	    // free memory
	    free(header_match);
	    free(read_match);
	    free(quality_match);
	    free(cigar);

	    sw_output_free(sw_output);
	  }
	}
      }
      if (array_list != NULL) array_list_free(array_list, NULL);
    }

    if (batch != NULL) aligner_batch_free(batch);
    if (item != NULL) list_item_free(item);

    if (time_on) { timing_stop(BATCH_WRITER, 0, timing_p); }
  } // end of batch loop                                                                                                                                                            

  // some cosmetic things before freeing the alignment,
  // to be sure to free twice
  alignment->query_name = NULL;
  alignment->sequence = NULL;
  alignment->quality = NULL;
  alignment->cigar = NULL;
  alignment_free(alignment);
  
  //  fclose(splice_fd);
  
  bam_fclose(bam_file);
  //  printf("END: batch_writer (total mappings %d)\n", total_mappings);
}

//------------------------------------------------------------------------------------

void batch_writer_input_init(char* match_filename, char* splice_exact_filename, 
			     char* splice_extend_filename, list_t* list_p, batch_writer_input_t* input_p) {
  input_p->match_filename = match_filename;
  input_p->splice_exact_filename = splice_exact_filename;
  input_p->splice_extend_filename = splice_extend_filename;
  input_p->list_p = list_p;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
