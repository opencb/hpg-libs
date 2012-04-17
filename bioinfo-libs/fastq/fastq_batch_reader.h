
#ifndef FASTQ_BATCH_READER_H
#define FASTQ_BATCH_READER_H

#include <stdio.h>
#include <stdio.h>
#include <pthread.h>

#include "fastq_file.h"
#include "fastq_read.h"
#include "fastq_batch_list.h"
#include "list.h"

//=====================================================
// structure and functions to manage fastq read batch server
//
//=====================================================

typedef struct fastq_batch_reader {
  int alive;
  pthread_mutex_t alive_lock;

  int eof;
  pthread_mutex_t eof_lock;

  size_t batch_size; // in bytes
  int batch_list_max_length;

  pthread_t thread;

  int source_id;

  fastq_file_t* fastq_file_p;
  fastq_batch_list_t* batch_list_p;
  list_t* qc_batch_list_p;
} fastq_batch_reader_t;


fastq_batch_reader_t* fastq_batch_reader_new(char* filename, int source_id, fastq_batch_list_t* batch_list_p, size_t batch_size, list_t* qc_batch_list_p, int batch_list_max_length);
void fastq_batch_reader_free(fastq_batch_reader_t* reader_p);

void fastq_batch_reader_start(fastq_batch_reader_t* reader_p);
unsigned int fastq_batch_reader_join(fastq_batch_reader_t* reader_p);

fastq_batch_list_item_t* fastq_batch_reader_next_batch(fastq_batch_reader_t* reader_p);

int fastq_batch_reader_get_eof(fastq_batch_reader_t* reader_p);
void fastq_batch_reader_set_eof(fastq_batch_reader_t* reader_p, int eof);

int fastq_batch_reader_get_alive(fastq_batch_reader_t* reader_p);
void fastq_batch_reader_set_alive(fastq_batch_reader_t* reader_p, int alive);


// thread function
void* fastq_batch_reader_thread_function(void* param_p);


#endif	/*  FASTQ_BATCH_READER_H  */
