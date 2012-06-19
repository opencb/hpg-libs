#ifndef BURROWS_WHEELER_H
#define BURROWS_WHEELER_H

#include <stdio.h>
#include <stdlib.h>

#include "fastq_read.h"
#include "alignment.h"
#include "BW_io.h"

typedef fastq_read_t read_t;
typedef alignment_t mapping_t;

//-----------------------------------------------------------------------------
// Read structure: id, sequence and quality
//-----------------------------------------------------------------------------

typedef struct qread {
  char *id;
  char *sequence;
  char *quality;
} qread_t;

qread_t *qread_new(const char *id, 
		   const char *sequence, 
		   const char *quality);
void qread_free(qread_t *read);

//-----------------------------------------------------------------------------
// Seeding method
//-----------------------------------------------------------------------------

typedef struct seeding {
  unsigned int seed_size;
} seeding_t;

seeding_t *seeding_new(const unsigned int seed_size);
void seeding_free(seeding_t *seeding);

//-----------------------------------------------------------------------------
// CAL policy
//-----------------------------------------------------------------------------

typedef struct cal_optarg {
  unsigned int min_num_seeds;
  unsigned int max_seed_distance;
  unsigned int left_flank_length;
  unsigned int right_flank_length;
} cal_optarg_t;

cal_optarg_t *cal_optarg_new(const unsigned int num_seeds, 
			     const unsigned int max_seed_distance,
			     const unsigned int left_flank_length, 
			     const unsigned int right_flank_length);
void cal_optarg_free(cal_optarg_t *optarg);

//-----------------------------------------------------------------------------

typedef struct cal {
  unsigned int chromosome_id;
  unsigned int strand;
  size_t start;
  size_t end;
} cal_t;

cal_t *cal_new(const unsigned int chromosome_id, 
	       const unsigned int strand,
	       const size_t start, 
	       const size_t end);
void cal_free(cal_t *cal);

//-----------------------------------------------------------------------------

typedef struct qread_cals {
  qread_t *qread;
  array_list_t *cal_list; // array list of cal_t structures
} qread_cals_t;

qread_cals_t *qread_cals_new(const qread_t *qread);
void qread_cals_free(qread_cals_t *qread_cals);


//-----------------------------------------------------------------------------
// Burrows-Wheeler Transform
//-----------------------------------------------------------------------------

typedef struct bwt_optarg {
  unsigned int num_errors;
  unsigned int num_threads;
  unsigned int max_alignments_per_read;
} bwt_optarg_t;

bwt_optarg_t *bwt_optarg_new(const unsigned int num_errors,
			     const unsigned int num_threads,
			     const unsigned int max_alginments_per_read);
void bwt_optarg_free(bwt_optarg_t *optarg);

//-----------------------------------------------------------------------------

typedef struct bwt_index {
  comp_matrix h_O, h_rO, h_Oi, h_rOi;
  vector h_C, h_rC, h_C1, h_rC1;
  comp_vector S, Si;
  exome karyotype;
  char *dirname;
} bwt_index_t;

bwt_index_t *bwt_index_new(const char *dirname);
void bwt_index_free(bwt_index_t *index);

//-----------------------------------------------------------------------------

typedef struct mapping {
  cal_t *cal;
  char *cigar;
} mapping_t;

//-----------------------------------------------------------------------------

typedef struct qread_mappings {
  qread_t *qread;
  array_list_t *mapping_list; // array list of mapping_t structures
} qread_mappings_t;

qread_cals_t *qread_cals_new(const qread_t *qread);
void qread_cals_free(qread_cals_t *qread_cals);


//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------

unsigned int bwt_map_seq_cpu(char *seq, 
			     bwt_optarg_t *bwt_optarg, 
			     bwt_index_t *index, 
			     array_list_t *mapping_list);

unsigned int bwt_map_read_cpu(read_t *read, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list);

unsigned int bwt_map_seqs_cpu(char **seqs, 
			      unsigned int num_reads,
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list);


unsigned int bwt_map_batch_cpu(read_batch_t *batch,
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       read_batch_t *unmapped_batch,
			       array_list_t *mapping_list);

//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_seq_cpu(char *seq, 
				   bwt_optarg_t *bwt_optarg, 
				   bwt_index_t *index, 
				   array_list_t *mapping_list);

unsigned int bwt_map_exact_read_cpu(read_t *read, 
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index, 
				    array_list_t *mapping_list);

unsigned int bwt_map_exact_seqs_cpu(char **seqs, 
				    unsigned int num_reads,
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index, 
				    array_list_t *mapping_list);


unsigned int bwt_map_exact_batch_cpu(read_batch_t *batch,
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index, 
				     read_batch_t *unmapped_batch,
				     array_list_t *mapping_list);

//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seq_cpu(char *seq, 
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index, 
				     array_list_t *mapping_list);

unsigned int bwt_map_inexact_read_cpu(read_t *read, 
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list);

unsigned int bwt_map_inexact_seqs_cpu(char **seqs, 
				      unsigned int num_reads,
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list);


unsigned int bwt_map_inexact_batch_cpu(read_batch_t *batch,
				       bwt_optarg_t *bwt_optarg, 
				       bwt_index_t *index, 
				       read_batch_t *unmapped_batch,
				       array_list_t *mapping_list);

//-----------------------------------------------------------------------------

unsigned int bwt_find_cals_from_seq_cpu(char *seq, 
					bwt_optarg_t *bwt_optarg, 
					bwt_index_t *index, 
					cal_optarg_t *cal_optarg, 
					array_list_t *cal_list);

unsigned int bwt_find_cals_from_read_cpu(read_t *read, 
					 bwt_optarg_t *bwt_optarg, 
					 bwt_index_t *index, 
					 cal_optarg_t *cal_optarg, 
					 array_list_t *cal_list);

unsigned int bwt_find_cals_from_seqs_cpu(char **seqs, 
					 unsigned int num_reads,
					 bwt_optarg_t *bwt_optarg, 
					 bwt_index_t *index, 
					 cal_optarg_t *cal_optarg, 
					 array_list_t *cal_list);


unsigned int bwt_find_cals_from_batch_cpu(read_batch_t *batch,
					  bwt_optarg_t *bwt_optarg, 
					  bwt_index_t *index, 
					  read_batch_t *unmapped_batch,
					  cal_optarg_t *cal_optarg, 
					  array_list_t *cal_list);


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#endif // BURROWS_WHEELER_H
