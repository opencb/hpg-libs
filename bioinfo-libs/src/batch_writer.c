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

void batch_writer_input_init(char* match_filename, char* splice_exact_filename, char* splice_extend_filename, list_t* list_p, batch_writer_input_t* input_p) {
  input_p->match_filename = match_filename;
  input_p->splice_exact_filename = splice_exact_filename;
  input_p->splice_extend_filename = splice_extend_filename;
  input_p->list_p = list_p;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
