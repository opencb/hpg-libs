
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <stdint.h>

#include "bam_writer.h"

//=====================================================
// functions to manage bam write server
//
//=====================================================

//----------------------------------------------
// bam writer new
//--------------------------------------------

bam_writer_t* bam_writer_new(char* filename, alignments_list_t* alignments_list_p, bam_header_t* bam_header_p, int mode, int chromosome) {  
  bam_writer_t* writer_p = (bam_writer_t*) calloc(1, sizeof(bam_writer_t));

  // open the output file....
  //
  writer_p->bam_file_p = bam_fopen_mode(filename, bam_header_p, "w");
  writer_p->mode = mode;
  writer_p->chromosome = chromosome;

  writer_p->alive = 1;
  //writer_p->alive_lock = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_init(&(writer_p->alive_lock), NULL);

  writer_p->alignments_list_p = alignments_list_p;
  writer_p->bam_write_alignment_count = (int*) calloc(NUM_OF_CHROMOSOMES, sizeof(int));

  bam_fwrite_header(bam_header_p, writer_p->bam_file_p);

  return writer_p;
}

//----------------------------------------------
// bam writer free
//--------------------------------------------

void bam_writer_free(bam_writer_t* writer_p) {

  // close the output file and exiting....
  //
  bam_fclose(writer_p->bam_file_p);
  free(writer_p->alignments_list_p);
}

//----------------------------------------------
// bam reader start
//--------------------------------------------

void bam_writer_start(bam_writer_t* writer_p) {

  // create and launch pthread....
  //
  LOG_DEBUG("Launching bam write server thread....\n");
  
  if (writer_p->mode == CHROMOSOME_MODE) {
      pthread_create(&(writer_p->thread), NULL, bam_writer_chromosome_thread_function, (void*) writer_p);
  } else if (writer_p->mode == SEQUENTIAL_MODE) {
      pthread_create(&(writer_p->thread), NULL, bam_writer_sequential_thread_function, (void*) writer_p);
  }
}

//----------------------------------------------

unsigned int bam_writer_join(bam_writer_t* writer_p) {
  void* r;
  pthread_join(writer_p->thread, &r);
  return (uintptr_t) r;
}

//----------------------------------------------
// 'secure' set and get for alive
//--------------------------------------------

void bam_writer_set_alive(bam_writer_t* writer_p, int alive) {
  pthread_mutex_lock(&(writer_p->alive_lock));
  writer_p->alive = alive;
  pthread_mutex_unlock(&(writer_p->alive_lock));
}

int bam_writer_get_alive(bam_writer_t* writer_p) {
  int alive;
  pthread_mutex_lock(&(writer_p->alive_lock));
  alive = writer_p->alive;
  pthread_mutex_unlock(&(writer_p->alive_lock));
  return alive;
}

//-----------------------------------------------------
// bam writer thread function
//
// this thread writes aligments to a BAM file in disk,
// 
//-----------------------------------------------------

//----------------------------------------------
// one thread for all chromosomes
// each chromosome is read sequentially
//----------------------------------------------
/*
void* bam_writer_sequential_thread_function(void* param_p) {

  struct timespec ts;
  ts.tv_sec = 0;
  ts.tv_nsec = 1000000;

  bam_writer_t* writer_p = (bam_writer_t*) param_p;
  bam1_t* bam_alignment_p;
    
  bam_writer_set_alive(writer_p, 1);
  
  printf("Thread-WRITE: START, for file '%s' and chromosome %i\n", writer_p->bam_file_p->filename, writer_p->chromosome);
  
  int actual_chromosome = 0;
  int write_bytes = 0;
  long batch_write_bytes = 0, total_write_bytes = 0;
   
  chrom_alignments_t chrom_alignments = (chrom_alignments_t) writer_p->bam_alignments_p[actual_chromosome];
// sleep(10); 
// printf("after sleep...\n");
  while (true) {

    // allocationg memory for the current fastq batch
    //
    //printf("Thread-WRITE: Starting while loop...\n");
    if (time_flag) { start_timer(t1_write); }

    bam_file_write_batch_by_chromosome(writer_p->bam_file_p, &(writer_p->bam_alignments_p[actual_chromosome]), &(writer_p->bam_write_alignment_count[actual_chromosome]), &(writer_p->bam_read_alignment_count[actual_chromosome]), actual_chromosome);

    if (time_flag) { stop_timer(t1_write, t2_write, write_time); }
//printf("bam_read_alignment_count[%i]: %i...\n", actual_chromosome, writer_p->bam_read_alignment_count[actual_chromosome]);
//printf("bam_write_alignment_count[%i]: %i...\n", actual_chromosome, writer_p->bam_write_alignment_count[actual_chromosome]);

    if ((writer_p->bam_read_alignment_count[actual_chromosome] == writer_p->bam_write_alignment_count[actual_chromosome])) {
     
      int chromosome_position_incr = 0;
      for (int i=actual_chromosome; i<NUM_OF_CHROMOSOMES; i++) {
	  if (writer_p->bam_read_alignment_count[i] != 0) {
	      chromosome_position_incr = i - actual_chromosome ;
	      break;
	  }
      }
     
     //printf("bam_reader_alive: %i\n", bam_reader_alive);
     
      if ((chromosome_position_incr == 0) && (!bam_reader_alive)) {
	  //printf("exit from writer file, bam_reader_alive: %i\n", bam_reader_alive);
	  break;
      } else if (chromosome_position_incr == 0) {
	  //printf("Thread-WRITE: go to sleep for a while...\n");
	  sched_yield();
	  usleep(10000);
      } else {
	  actual_chromosome += chromosome_position_incr;
	  for (int i=(actual_chromosome-chromosome_position_incr); i<actual_chromosome; i++) {
	      chrom_alignments_free(writer_p->bam_alignments_p[i], writer_p->max_estimated_alignments);
	  }
	  chrom_alignments = (chrom_alignments_t) writer_p->bam_alignments_p[actual_chromosome];
	  //bam_alignment_chromosome = (bam1_t**) writer_p->bam_alignment_matrix_p[actual_chromosome];
      }
      
    }
  } 
      
  bam_file_close(writer_p->bam_file_p);
  
  bam_writer_set_alive(writer_p, 0);
  
  printf("Thread-WRITE: END \n");
  
  pthread_exit(0);
}*/

//----------------------------------------------
// writer thread by splitted chromosome
//----------------------------------------------

void* bam_writer_chromosome_thread_function(void* param_p) {

  bam_writer_t* writer_p = (bam_writer_t*) param_p;
  bam1_t* alignment_p;
    
  bam_writer_set_alive(writer_p, 1);
  
  char log_message[200];
  sprintf(log_message, "Thread-WRITE: START, for file '%.150s' and chromosome %i\n", writer_p->bam_file_p->filename, writer_p->chromosome);
  LOG_DEBUG(log_message);  
  
  int current_chromosome = writer_p->chromosome;
  int num_alignments = 0;
  int write_bytes = 0, total_write_bytes = 0;

  chrom_alignments_t* chrom_alignments_p;

  while ( (chrom_alignments_p = alignments_list_get_chrom_alignment(current_chromosome, writer_p->alignments_list_p)) == NULL) {
    sched_yield();
    usleep(10000);      
  }

  while (1) {

    //printf("Thread-WRITE: Starting while loop...\n");
    if ((chrom_alignments_is_complete(chrom_alignments_p)) && (alignment_p = chrom_alignments_get_alignment(chrom_alignments_p, num_alignments)) != NULL) { 
      //printf("while... alignment_p->core.tid: %i, alignment_p->core.pos: %i\n", alignment_p->core.tid, alignment_p->core.pos);
      if (time_flag) { start_timer(t1_write); }
      write_bytes = bam_fwrite(alignment_p, writer_p->bam_file_p);
      total_write_bytes += write_bytes;
      num_alignments++;     
      if (time_flag) { stop_timer(t1_write, t2_write, write_time); }
            
      //printf("chrom_alignments_is_complete(chrom_alignments_p): %i, num_alignments: %i, chrom_alignments_p->alignment_count: %i\n", chrom_alignments_is_complete(chrom_alignments_p), num_alignments, chrom_alignments_p->alignment_count);
      
      if ((chrom_alignments_is_complete(chrom_alignments_p)) && (num_alignments == chrom_alignments_p->alignment_count)) {
	break;
      } else if (!chrom_alignments_is_complete(chrom_alignments_p)) {
	//printf("Thread-WRITE: go to sleep for a while...\n");
	sched_yield();
	usleep(10000);
      }
      
    } else {
      //printf("waiting for alignments..., num_alignments: %i\n", num_alignments);
      if ((chrom_alignments_is_complete(chrom_alignments_p)) && (num_alignments == chrom_alignments_p->alignment_count)) {
	break;
      } else {
	//printf("Thread-WRITE: go to sleep for a while...\n");
	sched_yield();
	usleep(10000);
      }    
    }

  } 
      
  bam_fclose(writer_p->bam_file_p);
  
  bam_writer_set_alive(writer_p, 0);
  
  sprintf(log_message, "Thread-WRITE: END for chromosome %i \n", current_chromosome);
  LOG_DEBUG(log_message);    
  
  //printf("Thread-WRITE: total alignments for chromosome %i: %i \n", total_alignments, total_alignments);

  pthread_exit((void*)num_alignments);
}

//----------------------------------------------
// writer thread, one file for all chromosomes
//----------------------------------------------

void* bam_writer_sequential_thread_function(void* param_p) {
  
  struct timespec ts;
  ts.tv_sec = 0;
  ts.tv_nsec = 1000000;

  bam_writer_t* writer_p = (bam_writer_t*) param_p;
  bam1_t* alignment_p;
    
  bam_writer_set_alive(writer_p, 1);
  
  char log_message[200];
  sprintf(log_message, "Thread-WRITE: START, for file '%.150s' and chromosome %i\n", writer_p->bam_file_p->filename, writer_p->chromosome);
  LOG_DEBUG(log_message);  
  
  int current_chromosome = 0;
  int num_alignments = 0;
  int write_bytes = 0, total_write_bytes = 0;

  chrom_alignments_t* chrom_alignments_p;

  while ( (chrom_alignments_p = alignments_list_get_chrom_alignment(current_chromosome, writer_p->alignments_list_p)) == NULL) {
    sched_yield();
    usleep(10000);   
  }

  while (1) {

    //printf("Thread-WRITE: Starting while loop... for chromosome: %i\n", current_chromosome);
    //if ((current_chromosome < NUM_OF_CHROMOSOMES) && (chrom_alignments_is_complete(chrom_alignments_p)) && (alignment_p = chrom_alignments_get_alignment(chrom_alignments_p, num_alignments)) != NULL) { 
    if ((current_chromosome < NUM_OF_CHROMOSOMES) && (chrom_alignments_is_complete(chrom_alignments_p))) { 
      
      if (time_flag) { start_timer(t1_write); }
      
      total_write_bytes = bam_fwrite_sorted_array(chrom_alignments_p->bam_alignments_p, chrom_alignments_p->indices_p, chrom_alignments_p->alignment_count, writer_p->bam_file_p);
      num_alignments += chrom_alignments_p->alignment_count;
      
      if (time_flag) { stop_timer(t1_write, t2_write, write_time); }

      current_chromosome++;
      
      while ( (chrom_alignments_p = alignments_list_get_chrom_alignment(current_chromosome, writer_p->alignments_list_p)) == NULL) {
	sched_yield();
	usleep(10000);      
      }
    } else if (current_chromosome == NUM_OF_CHROMOSOMES) {
      break;
    } else {
      sched_yield();
      usleep(10000);  
    }


  } 
      
  bam_fclose(writer_p->bam_file_p);
  
  bam_writer_set_alive(writer_p, 0);
  
  sprintf(log_message, "Thread-WRITE: END with total alignments: %i \n", num_alignments);
  LOG_DEBUG(log_message);  

  pthread_exit((void*)num_alignments);
}