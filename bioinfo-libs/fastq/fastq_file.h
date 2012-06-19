
#ifndef FASTQ_FILE_H
#define FASTQ_FILE_H

#include <stdio.h>

//#include "qc_batch.h"
#include "file_utils.h"

#include "fastq_batch.h"
#include "fastq_read.h"

//====================================================================================
//  fastq_t.h
//
//  fastq_t structures and prototypes
//====================================================================================

#define MAX_FASTQ_FILENAME_LENGTH		64		// Maximum filenname length
#define MAX_READ_ID_LENGTH			256		// Maximum read ID length
#define MAX_READ_SEQUENCE_LENGTH		2048	// Maximum read sequence length

#define MAX_NUM_PRODUCERS			10

#define FQ_SEEK_BEGIN 	0
#define FQ_SEEK_CURR  	1
#define FQ_SEEK_RND   	2


typedef struct fastq_file {
	char* filename;
	//int source_id;
	char* mode;
	FILE* fd;
	unsigned long num_reads;
	unsigned long num_lines;
	char* quality_encoding;	// illumina v1.5, solid, ...
} fastq_file_t;


typedef struct source {
	int id;
	char filename[MAX_FASTQ_FILENAME_LENGTH];
} source_t;


fastq_file_t *fastq_fopen(char *filename);

fastq_file_t *fastq_fopen_mode(char *filename, char *mode);

int fastq_fread(fastq_read_t *read, fastq_file_t *fq_file);

int fastq_fread_num_reads(fastq_read_t* buffer_reads, int num_reads, fastq_file_t *fq_file);

int fastq_fread_max_size(fastq_read_t *buffer_fq_reads, unsigned long max_size, fastq_file_t *fq_file);

int fastq_fread_batch_max_size(fastq_batch_t *buffer_fq_read_batch, unsigned long max_size, fastq_file_t *fq_file);

int fastq_fread_index_positions(fastq_read_t* buffer_reads, int *index_positions, fastq_file_t *fq_file);

int fastq_fwrite(fastq_read_t* buffer_reads, int num_writes, fastq_file_t *fq_file);

unsigned int fastq_fcount(fastq_file_t *fq_file);

//void fastq_remove_Ns(fastq_read_t* buffer_reads, qc_read_t* qc_read, int max_N_per_read);

void fastq_fclose(fastq_file_t *fq_file);

#endif    // FASTQ_FILE_H
