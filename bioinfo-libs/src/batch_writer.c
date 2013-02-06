#include "batch_writer.h"

//------------------------------------------------------------------------------------
// this functions writes reads to disk, reads came from
// a list
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void batch_writer2(batch_writer_input_t* input) {

  LOG_DEBUG_F("START: batch_writer (%i): START, for file %s\n", 
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
  size_t limit_print = 1000000;
  // main loop
  while ( (item = list_remove_item(write_list)) != NULL ) {
    //    if (array_list == NULL) printf("batch_writer.c...\n");
    //printf("######Extract item Write\n");
    total_batches++;
    batch = (mapping_batch_t *) item->data_p;

    num_reads = array_list_size(batch->fq_batch);
    total_reads += num_reads;
    //printf("Read %i reads\n", num_reads);
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
	//printf("\tWRITE : read %d (%d items): mapped...\n", i, num_items);
	if (num_items > 1000){alig = (alignment_t *) array_list_get(0, array_list); printf("%s\n", alig->sequence); }
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
      LOG_DEBUG_F("TOTAL READS PROCESS: %lu\n", total_reads);
      LOG_DEBUG_F("\tTotal Reads Mapped: %lu(%.2f%)\n", num_mapped_reads, (float)(num_mapped_reads*100)/(float)(total_reads));
      limit_print += 1000000;
    }
    
    //printf("Batch Write OK!\n");
    
    if (batch != NULL) mapping_batch_free(batch);
    if (item != NULL) list_item_free(item);
    
    //if (time_on) { timing_stop(BATCH_WRITER, 0, timing_p); }
  } // end of batch loop                                                                                           
  
  bam_fclose(bam_file);

  basic_statistics_init(total_reads, num_mapped_reads, total_mappings, &basic_st);
  LOG_DEBUG_F("Batch_writer  (Total batches %lu): END\n", total_batches);
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
