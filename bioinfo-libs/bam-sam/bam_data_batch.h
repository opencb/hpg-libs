/*
 * bam_data_batch.h
 *
 *  Created on: Aug 5, 2011
 *      Author: victor
 */

#ifndef BAM_DATA_BATCH_H
#define BAM_DATA_BATCH_H

#include <pthread.h>
#include <stdint.h>

#include "bam_file.h"
#include "commons.h"

#define MAX_CIGAR_POSITIONS	10


//====================================================================================
//  bam_data_batch.h
//
//  structures and methods for inserting, removing and manage bam_data_batch objects
//====================================================================================

typedef struct bam_data_core {
  int map_quality;
  short int strand;  
  short int chromosome;
  int start_coordinate;
  int alignment_length;
  //int pair_distance;
  int id_seq_index;
  int cigar_index;
  short int paired_end;
} bam_data_core_t;

typedef struct bam_data_batch {
  int num_alignments;
  int num_cigar_operations;
  int id_seq_length;
  int last_alignments_position;
  short int num_chromosomes_in_batch;
  short int* chromosomes;
  int* start_positions;
  int* end_positions;
  bam_data_core_t* core_data_p;
  char* id_seq_data_p;
  uint32_t* cigar_data_p;
} bam_data_batch_t;


// ------------------------------------------------
// bam_data_batch functions
// ------------------------------------------------

bam_data_batch_t* bam_data_batch_new(int num_alignments);
bam_data_batch_t* bam_data_batch_init(bam_data_batch_t* bam_data_batch_p, bam_batch_t* bam_batch_p);
bam_data_batch_t* bam_data_batch_by_chromosome_init(bam_data_batch_t* bam_data_batch_p, bam_batch_t* bam_batch_p, int start_alignment);
void bam_data_batch_free(bam_data_batch_t* bam_data_batch_p);

bam_data_core_t* bam_data_batch_get_core(bam_data_batch_t* bam_data_batch_p);
uint32_t* bam_data_batch_get_cigar_data(bam_data_batch_t* bam_data_batch_p);
char* bam_data_batch_get_id_seq_data(bam_data_batch_t* bam_data_batch_p);

#endif /* BAM_DATA_BATCH_H_ */

























