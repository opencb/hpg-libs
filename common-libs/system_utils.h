
#ifndef SYSTEM_UTILS_H
#define SYSTEM_UTILS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "commons.h"
#include "cuda_commons.h"
#include "log.h"

#define FASTQ_QC	1
#define FASTQ_PREPRO	2
#define BAM_QC		3

#define QC_MEMORY_USAGE_FACTOR 		4
#define PREPRO_MEMORY_USAGE_FACTOR	6

#define BAM_BATCH_SIZE_TO_FREE_MEMORY_RATIO	20

#define BAM_SAM_COMPRESSION_RATIO	4
#define MEAN_COMPRESSED_ALIGNMENT_SIZE 	100

unsigned long int get_free_memory();
unsigned long int get_estimated_memory_needed(int process, int batch_size, int max_list_length);
int get_max_estimated_alignments_by_chromosome(char* input_filename);
int get_optimal_cpu_num_threads();
int get_optimal_gpu_num_threads();
int get_optimal_batch_size(int process, int max_list_length);

#endif	/*  SYSTEM_UTILS_H   */
