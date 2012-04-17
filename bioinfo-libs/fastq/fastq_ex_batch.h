
#ifndef FASTQ_EX_BATCH_H
#define FASTQ_EX_BATCH_H

#include <stdio.h>

#include "fastq_file.h"

typedef struct fastq_ex_batch {
  int source_id; 
  unsigned long num_reads;
  unsigned long data_indices_size; // data_indices_size = seq_indices_size = quality_indices_size
  unsigned long data_size;         // data_size = seq_size = quality_size
  int *data_indices;               // data_indices = seq_indices = quality_indices
  int *header_indices;
  char *header;
  char *seq;
  char *quality;
} fastq_ex_batch_t;


fastq_ex_batch_t* fastq_ex_batch_new(unsigned long size);
void fastq_ex_batch_free(fastq_ex_batch_t* fastq_ex_batch_p);

int fastq_ex_batch_read(fastq_ex_batch_t* batch_p, int encode, fastq_file_t* file_p);
int fastq_ex_batch_write(fastq_ex_batch_t* batch_p, int decode, fastq_file_t* file_p);

/*
void fastq_ex_batch_print(fastq_ex_batch_t* fastq_ex_batch_p, FILE* fd);

void fprintf_read(FILE* fd, fastq_ex_batch_t* batch_p, int index);

void fprintf_rtrim_read(FILE* fd, fastq_ex_batch_t* batch_p, int index, int rtrim_length);

void fprintf_ltrim_read(FILE* fd, fastq_ex_batch_t* batch_p, int index, int ltrim_length);

void fprintf_trim_read(FILE* fd, fastq_ex_batch_t* batch_p, int index, int rtrim_length, int ltrim_length);
*/
#endif  /* FASTQ_EX_BATCH_H */ 
