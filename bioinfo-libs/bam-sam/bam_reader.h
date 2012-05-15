
#ifndef BAM_READER_H
#define BAM_READER_H

#include <stdio.h>
#include <pthread.h>

#include "bam_data_batch_list.h"
#include "bam_file.h"

#include "chrom_alignments.h"
#include "commons.h"
#include "list.h"
#include "system_utils.h"


//===========================================================
// structure and functions to manage bam read server
//
//===========================================================

#define SEQUENTIAL_MODE		0
#define CHROMOSOME_MODE		1
#define LIST_INSERT_MODE	2

#define NO_SORT 0
#define SORT_BY_POSITION 	1
#define SORT_BY_ID 		2

#define MIN_FREE_MEMORY_ALLOWED_FOR_BAM_READ_KB		200

typedef struct bam_reader {
  int alive;
  int mode;
  int sort;
  int chromosome;
  int max_estimated_alignments;
  size_t batch_size;
  int base_quality;
  pthread_mutex_t alive_lock;
  pthread_t thread;
  
  bam_file_t* bam_file_p;
  alignments_list_t* alignments_list_p;
  list_t* bam_batch_list_p;
  list_t* bam_data_batch_list_p;  
} bam_reader_t;

//bam_reader_t* bam_reader_new(char* filename, int max_estimated_alignments, chrom_alignments_t* bam_alignments_p, int* bam_read_alignment_count_p, int split, int chromosome);
bam_reader_t* bam_reader_new(char* filename, size_t batch_size, int base_quality, alignments_list_t* alignments_list_p, int mode, int sort, int chromosome);
bam_reader_t* bam_reader_by_batch_new(char * filename, size_t batch_size, int base_quality, list_t* bam_data_batch_list_p, int mode);
void bam_reader_free(bam_reader_t* reader_p);

void bam_reader_start(bam_reader_t* reader_p);
unsigned int bam_reader_join(bam_reader_t* reader_p);

int bam_reader_get_alive(bam_reader_t* reader_p);
void bam_reader_set_alive(bam_reader_t* reader_p, int alive);


// thread function
//
void* bam_reader_sequential_thread_function(void* param_p);
void* bam_reader_chromosome_thread_function(void* param_p);
void* bam_reader_list_insert_thread_function(void* param_p);
void* bam_reader_list_insert_by_chromosome_thread_function(void* param_p);
//void* bam_reader_thread_function(void* param_p);

#endif
