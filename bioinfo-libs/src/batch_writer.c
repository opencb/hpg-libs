#include "batch_writer.h"

//------------------------------------------------------------------------------------
// this functions writes reads to disk, reads came from
// a list
//------------------------------------------------------------------------------------
/*
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
  
  //printf("HEADER FROM WRITE: %s\n", input_p->header_filename);
  //  bam_header_p = bam_header_new(HUMAN, NCBI37, input_p->header_filename);
  bam_header_p = create_bam_header_by_genome(input_p->genome);

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
      for(size_t i = 0; i < batch_p->size; i++){
	//alignment_print(buffer_p[i]);
	//if (!buffer_p[i]) {exit(-1);}
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
*/
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

  printf("START: batch_writer (%i): START, for file %s\n", 
	 omp_get_thread_num(), input->match_filename);

  bam1_t *bam1;
  bam_header_t *bam_header;
  bam_file_t *bam_file;
  alignment_t *alig;

  char* match_filename = input->match_filename;
  //  char* splice_filename = input->splice_filename;
  
  list_t *write_list = input->list_p;
  array_list_t *array_list;
  
  list_item_t *item = NULL;
  mapping_batch_t *batch = NULL;
  //  fastq_batch_t *fq_batch = NULL;

  FILE* fd;

  fastq_read_t *fq_read;

  char *id, *sequence, *quality;
  static char aux[8096];
  
  size_t read_index, read_len, header_len;

  //  bam_header = bam_header_new(HUMAN, NCBI37, input_p->header_filename);
  bam_header = create_bam_header_by_genome(input->genome);
  bam_file = bam_fopen_mode(match_filename, bam_header, "w");

  bam_fwrite_header(bam_header, bam_file);

  size_t num_reads = 0, num_mapped_reads = 0, num_items = 0;
  size_t total_reads = 0, total_mappings = 0, total_batches = 0;
  size_t reads_mapped = 0;
  size_t limit_print = 500000;
  // main loop
  while ( (item = list_remove_item(write_list)) != NULL ) {
    //    if (array_list == NULL) printf("batch_writer.c...\n");
    //printf("######Extract item Write\n");
    total_batches++;
    batch = (mapping_batch_t *) item->data_p;

    num_reads = array_list_size(batch->fq_batch);
    total_reads += num_reads;
    
    for (size_t i = 0; i < num_reads; i++) {
      array_list = batch->mapping_lists[i];
      num_items = array_list_size(array_list);
      total_mappings += num_items;
      // mapped or not mapped ?
      if (num_items == 0) {
	total_mappings++;

	fq_read = (fastq_read_t *) array_list_get(i, batch->fq_batch);

	// calculating cigar
	sprintf(aux, "%luX", fq_read->length);
	
	alig = alignment_new();

	header_len = strlen(fq_read->id);
        id = (char *)malloc(sizeof(char)*(header_len + 1));
	get_to_first_blank(fq_read->id, header_len, id);
	//free(fq_read->id);
	alignment_init_single_end(id, fq_read->sequence, fq_read->quality, 
				  0, -1, -1, aux, 1, 0, 0, 0, 0, NULL, alig);

	bam1 = convert_to_bam(alig, 33);
	bam_fwrite(bam1, bam_file);
	bam_destroy1(bam1);

	// some cosmetic stuff before freeing the alignment,
	// (in order to not free twice some fields)
	//alig->query_name = NULL;
	alig->sequence = NULL;
	alig->quality = NULL;
	alig->cigar = NULL;
	alignment_free(alig);

	//	printf("\tWRITE : read %i (%d items): unmapped...done !!\n", i, num_items);

      } else {
	num_mapped_reads++;
	//	printf("\tWRITE : read %d (%d items): mapped...\n", i, num_items);
	for (size_t j = 0; j < num_items; j++) {
	  alig = (alignment_t *) array_list_get(j, array_list);
	  //printf("\t%s\n", alig->cigar);
	  if (alig != NULL) {
	    
	    bam1 = convert_to_bam(alig, 33);
	    bam_fwrite(bam1, bam_file);
	    bam_destroy1(bam1);
	  
	    alignment_free(alig);
	  }
	}
	//	printf("\tWRITE : read %d (%d items): mapped...done !!\n", i, num_items);
      }
      if (array_list) array_list_free(array_list, NULL);
    }
    if (total_reads >= limit_print) {
      printf("TOTAL READS PROCESS: %lu\n", total_reads);
      printf("\tTotal Reads Mapped: %lu(%.2f%)\n", num_mapped_reads, (float)(num_mapped_reads*100)/(float)(total_reads));
      limit_print += 500000;
    }
    //printf("Batch Write OK!\n");
    
    if (batch != NULL) mapping_batch_free(batch);
    if (item != NULL) list_item_free(item);
    
    //if (time_on) { timing_stop(BATCH_WRITER, 0, timing_p); }
  } // end of batch loop                                                                                           
  bam_fclose(bam_file);

  basic_statistics_init(total_reads, num_mapped_reads, total_mappings, &basic_st);
}

//------------------------------------------------------------------------------------

bam_header_t *create_bam_header_by_genome(genome_t *genome) {

  bam_header_t *bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));

  int num_targets = genome->num_chromosomes;

  bam_header->n_targets = num_targets;
  bam_header->target_name = (char **) calloc(num_targets, sizeof(char *));
  bam_header->target_len = (uint32_t*) calloc(num_targets, sizeof(uint32_t));
  for (int i = 0; i < num_targets; i++) {
    bam_header->target_name[i] = strdup(genome->chr_name[i]);
    bam_header->target_len[i] = genome->chr_size[i] + 1;
  }
  bam_header->text = strdup("@PG\tID:HPG-Aligner\tVN:1.0\n");
  bam_header->l_text = strlen(bam_header->text);

  return bam_header;
}

//------------------------------------------------------------------------------------

void batch_writer_input_init(char* match_filename, char* splice_exact_filename, 
			     char* splice_extend_filename, 
			     list_t* list_p, genome_t* genome, 
			     batch_writer_input_t* input_p) {

  input_p->match_filename = match_filename;
  input_p->splice_exact_filename = splice_exact_filename;
  input_p->splice_extend_filename = splice_extend_filename;
  input_p->list_p = list_p;

  input_p->genome = genome;

}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
