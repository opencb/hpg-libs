
#ifndef FASTQ_BATCH_READER_OMP_H
#define FASTQ_BATCH_READER_OMP_H

#include <stdio.h>
#include <stdio.h>
#include <pthread.h>

#include "fastq_file.h"
#include "fastq_read.h"
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
  list_t* batch_list_p;
} fastq_batch_reader_t;


fastq_batch_reader_t* fastq_batch_reader_new(char* filename, int source_id, list_t* batch_list_p,  size_t batch_size, int batch_list_max_length);
void fastq_batch_reader_free(fastq_batch_reader_t* reader_p);

// thread function
void* fastq_batch_reader_thread_function(void* param_p);

#endif	/*  FASTQ_BATCH_READER_OMP_H  */
