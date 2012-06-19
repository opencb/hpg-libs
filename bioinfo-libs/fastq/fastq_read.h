
#ifndef _FASTQ_READ_H
#define _FASTQ_READ_H

#define ID_MAX_BYTES		1*1024*1024
#define INDEX_MAX		100000
#define SEQUENCE_MAX_BYTES	2*1024*1024

//--------------------------------------------
//  structures
//--------------------------------------------

typedef struct fastq_read {
  char *id;
  char *sequence;
  char *quality;
} fastq_read_t;

//--------------------------------------------

fastq_read_t *fastq_read_new(char *id, char *sequence, char *quality);
void fastq_read_free(fastq_read_t *fq_read);

//--------------------------------------------
//--------------------------------------------


typedef struct fastq_read_stats {
  char *id;
  short int seq_length;
  float quality_average;
} fastq_read_stats_t;

//--------------------------------------------
//  functions
//--------------------------------------------

float fastq_quality_average(fastq_read_t fq_read_t);

void fastq_nt_quality_average(fastq_read_t fq_read_t);

void fastq_nt_count(fastq_read_t fq_read_t);

void fastq_length(fastq_read_t fq_read_t);

void fastq_remove_barcode(fastq_read_t fq_read_t, char *barcode);

void fastq_stats(fastq_read_t fq_read_t, fastq_read_stats_t *read_stats_t);


#endif
