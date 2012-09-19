#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

//#include "list.h"
//

//#include "fastq_file.h"
//#include "fastq_ex_batch.h"
#include "fastq_batch_reader.h"

//------------------------------------------------------------------------------------
// init structure for the batch reader
//------------------------------------------------------------------------------------

void fastq_batch_reader_input_init(char *filename1, char *filename2,
				   int flags, int batch_size, list_t *list, 
				   fastq_batch_reader_input_t *input) {
  input->filename1 = filename1;
  input->filename2 = filename2;
  input->flags = flags;
  input->batch_size = batch_size;
  input->list = list;
}

//------------------------------------------------------------------------------------
// this functions reads read from disk, and save them
// into a list to be processed by the next thread
//------------------------------------------------------------------------------------

void fastq_batch_reader_single(fastq_batch_reader_input_t* input);
void fastq_batch_reader_pair(fastq_batch_reader_input_t* input);

//------------------------------------------------------------------------------------

void fastq_batch_reader(fastq_batch_reader_input_t* input) {
  
  if (input->filename1 != NULL && 
      input->filename2 != NULL &&
      (input->flags & PAIRED_END_FLAG ||
       input->flags & MATE_PAIR_FLAG)) {
    
    fastq_batch_reader_pair(input);
    
  } else if (input->filename1 != NULL ) {
    
    fastq_batch_reader_single(input);
    
  } else {
    printf("ERROR: FastQ batch reader\n");
    exit(-1);
  }
}

//------------------------------------------------------------------------------------

void fastq_batch_reader_single(fastq_batch_reader_input_t* input) {
  unsigned int total_reads = 0;

  /*struct timespec ts;
  ts.tv_sec = 1;
  ts.tv_nsec = 0;*/

  char *filename = input->filename1;
  int flags = input->flags & 3;
  int batch_size = input->batch_size;
  list_t *list = input->list;

  printf("fastq_batch_reader (%i): START, for file %s\n", 
	 omp_get_thread_num(), filename);
		
  int num_reads = 0, num_batches = 0;
  fastq_batch_t *batch;
  list_item_t *item = NULL;

  fastq_file_t *file = fastq_fopen(filename);
  
  while (1) {
    //    if (time_on) { timing_start(FASTQ_READER, 0, timing_p); }
    // allocationg memory for the current fastq batch
    //
    batch = fastq_batch_new(batch_size);
    //batch_p->source_id = reader_p->source_id;
      
    // read reads from file
    //
    num_reads = fastq_fread_batch_max_size(batch, batch_size, file);
    total_reads	+= num_reads;
    
    // if no reads, free memory and go out....
    //
    if (num_reads==0) {
      fastq_batch_free(batch);
      break;
    }

    // otherwise, create a new batch object..
    // and insert this batch to the corresponding list
    //
    item = list_item_new(num_batches, flags, batch);
    //    if (time_on) { timing_stop(FASTQ_READER, 0, timing_p); }

    list_insert_item(item, list);
  
    //printf("fastq_ex_batch_reader: reading batch %i done !!\n", num_batches);
    num_batches++;
    
    //if (num_batches==2) break;
  } // end of batch loop
  
  list_decr_writers(list);
  fastq_fclose(file);
  
  //  if (statistics_on) { 
  //    statistics_set(FASTQ_READER_ST, 0, num_batches, statistics_p); 
  //    statistics_set(FASTQ_READER_ST, 1, total_reads, statistics_p); 
  //    statistics_set(TOTAL_ST, 0, total_reads, statistics_p); 
  //  }
  
  printf("fastq_batch_reader: END, %i total reads (%i batches), for file %s\n", 
	 total_reads, num_batches, filename);
}

//------------------------------------------------------------------------------------

void fastq_batch_reader_pair(fastq_batch_reader_input_t* input) {
  unsigned int total_reads = 0;
  unsigned int total_batches = 0;

  /*struct timespec ts;
  ts.tv_sec = 1;
  ts.tv_nsec = 0;*/

  char *filename1 = input->filename1;
  char *filename2 = input->filename2;
  int flags = input->flags & 3;
  int batch_size = input->batch_size;
  list_t *list = input->list;

  printf("fastq_batch_reader (%i): START, for files %s, %s\n", 
	 omp_get_thread_num(), filename1, filename2);
		
  int num_reads1 = -1, num_batches1 = 0;
  int num_reads2 = -1, num_batches2 = 0;
  fastq_batch_t *batch;
  list_item_t *item = NULL;

  fastq_file_t *file1 = fastq_fopen(filename1);
  fastq_file_t *file2 = fastq_fopen(filename2);
  
  while (1) {
    //    if (time_on) { timing_start(FASTQ_READER, 0, timing_p); }

    // allocationg memory for the current fastq batch
    batch = fastq_batch_new(batch_size);
      
    // read reads from file #1
    if (num_reads1 != 0) {
      num_reads1 = fastq_fread_batch_max_size(batch, batch_size, file1);
    }
    
    if (num_reads1 > 0) {
      num_batches1++;
      total_reads += num_reads1;
      item = list_item_new(num_batches1, flags | PAIR1_FLAG, batch);

      //      if (time_on) { timing_stop(FASTQ_READER, 0, timing_p);  }
      list_insert_item(item, list);
      //      if (time_on) { timing_start(FASTQ_READER, 0, timing_p); }

      // allocationg memory for the current fastq batch
      batch = fastq_batch_new(batch_size);
    }

    // read reads from file #2
    if (num_reads2 != 0) {
      num_reads2 = fastq_fread_batch_max_size(batch, batch_size, file2);
    }
    
    if (num_reads2 > 0) {
      num_batches2++;
      total_reads += num_reads2;
      item = list_item_new(num_batches2, flags | PAIR2_FLAG, batch);

      //      if (time_on) { timing_stop(FASTQ_READER, 0, timing_p);  }
      list_insert_item(item, list);
      //      if (time_on) { timing_start(FASTQ_READER, 0, timing_p); }
    }

    // if no reads, free memory and go out....
    //
    if (num_reads1 == 0 && num_reads2 == 0) {
      fastq_batch_free(batch);
      break;
    }
    //if (num_batches==2) break;
  } // end of batch loop
  
  list_decr_writers(list);
  fastq_fclose(file1);
  fastq_fclose(file2);
  
  //  if (statistics_on) { 
  //    statistics_set(FASTQ_READER_ST, 0, num_batches1 + num_batches2, statistics_p); 
  //    statistics_set(FASTQ_READER_ST, 1, total_reads, statistics_p); 
  //    statistics_set(TOTAL_ST, 0, total_reads, statistics_p); 
  //  }
  
  printf("fastq_batch_reader: END, %i total reads (%i batches), for files %s, %s\n", 
	 total_reads, num_batches1 + num_batches2, filename1, filename2);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
