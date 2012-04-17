
#ifndef BAM_WRITER_H
#define BAM_WRITER_H

#include <stdio.h>
#include <pthread.h>

#include "bam.h"
#include "bam_commons.h"
#include "bam_file.h"
#include "chrom_alignments.h"
#include "commons.h"
#include "system_utils.h"


//===========================================================
// structure and functions to manage bam write server
//
//===========================================================

#define SEQUENTIAL_MODE	0
#define CHROMOSOME_MODE	1


typedef struct bam_writer {
  int alive;
  int mode;
  int chromosome;
  pthread_mutex_t alive_lock;
  pthread_t thread;
  
  bam_file_t* bam_file_p;
  alignments_list_t* alignments_list_p;
  int* bam_write_alignment_count;
} bam_writer_t;

bam_writer_t* bam_writer_new(char* filename, alignments_list_t* alignments_list_p, bam_header_t* bam_header_p, int mode, int chromosome);
void bam_writer_free(bam_writer_t* writer_p);

void bam_writer_start(bam_writer_t* writer_p);
unsigned int bam_writer_join(bam_writer_t* writer_p);

int bam_writer_get_alive(bam_writer_t* writer_p);
void bam_writer_set_alive(bam_writer_t* writer_p, int alive);


// thread function
//
void* bam_writer_sequential_thread_function(void* param_p);
void* bam_writer_chromosome_thread_function(void* param_p);



#endif
