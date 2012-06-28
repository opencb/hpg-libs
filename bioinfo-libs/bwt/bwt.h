#ifndef BURROWS_WHEELER_H
#define BURROWS_WHEELER_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

#include "array_list.h"
#include "fastq_read.h"
#include "fastq_batch.h"
#include "alignment.h"
#include "BW_io.h"
#include "BW_search.h"

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

double global_parallel, global_sequential;

//-----------------------------------------------------------------------------
// Parameters for seed
//-----------------------------------------------------------------------------
typedef struct seed {
  unsigned int starts;
  unsigned int end;
}seed_t;

void seed_init(unsigned int start, unsigned int end, seed_t *seed_p);

//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Paratemers for the candidate alignment localizations (CALs)
//-----------------------------------------------------------------------------


typedef struct cal_optarg {
  unsigned int min_cal_size;
  unsigned int max_cal_distance;
  unsigned int seed_size;
  unsigned int min_seed_size;
} cal_optarg_t;

cal_optarg_t *cal_optarg_new(const unsigned int min_cal_size, 
			     const unsigned int max_cal_distance, 
			     const unsigned int seed_size,
			     const unsigned int min_seed_size);

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

typedef cal_t region_t;

region_t *region_new(const unsigned int chromosome_id, 
	             const unsigned int strand,
	             const size_t start, 
	             const size_t end);

void region_free(region_t *cal);

//-----------------------------------------------------------------------------

typedef struct read_cals {
  fastq_read_t *read;
  array_list_t *cal_list; // array list of cal_t structures
} read_cals_t;

read_cals_t *read_cals_new(const fastq_read_t *read);
void read_cals_free(read_cals_t *read_cals);


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

typedef struct read_mappings {
  fastq_read_t *read;
  array_list_t *mapping_list; // array list of mapping_t structures
} read_mappings_t;

read_cals_t *read_cals_new(const fastq_read_t *read);
void read_cals_free(read_cals_t *read_cals);

//-------------------------------------------------------------------

//-----------------------------------------------------------------------------
// general functions
//-----------------------------------------------------------------------------

char * bwt_error_type(char error_kind);

unsigned int bwt_map_seq_cpu(char *seq, 
			     bwt_optarg_t *bwt_optarg, 
			     bwt_index_t *index, 
			     array_list_t *mapping_list);

unsigned int bwt_map_read_cpu(fastq_read_t *read, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list);

unsigned int bwt_map_seqs_cpu(char **seqs, 
			      unsigned int num_reads,
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list);

unsigned int bwt_map_batch_cpu(fastq_batch_t *batch,
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       fastq_batch_t *unmapped_batch,
			       array_list_t *mapping_list);

//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_seq_cpu(char *seq, 
				   bwt_optarg_t *bwt_optarg, 
				   bwt_index_t *index, 
				   array_list_t *mapping_list);

unsigned int bwt_map_exact_read_cpu(fastq_read_t *read, 
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index, 
				    array_list_t *mapping_list);

unsigned int bwt_map_exact_seqs_cpu(char **seqs, 
				    unsigned int num_reads,
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index, 
				    array_list_t *mapping_list);

unsigned int bwt_map_exact_batch_cpu(fastq_batch_t *batch,
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index, 
				     fastq_batch_t *unmapped_batch,
				     array_list_t *mapping_list);

//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seeds_seq_cpu(char *seq, seed_t *seeds, unsigned int num_seeds,
					bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
					array_list_t *mapping_list);

unsigned int bwt_map_inexact_seed_cpu(char *seq,
                                     bwt_optarg_t *bwt_optarg,
                                     bwt_index_t *index,
                                     array_list_t *mapping_list);

unsigned int bwt_map_inexact_seq_cpu(char *seq, 
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index, 
				     array_list_t *mapping_list);

unsigned int bwt_map_inexact_read_cpu(fastq_read_t *read, 
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list);

unsigned int bwt_map_inexact_seqs_cpu(char **seqs, 
				      unsigned int num_reads,
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list);

unsigned int bwt_map_inexact_batch_cpu(fastq_batch_t *batch,
				       bwt_optarg_t *bwt_optarg, 
				       bwt_index_t *index, 
				       fastq_batch_t *unmapped_batch,
				       array_list_t *mapping_list);

//-----------------------------------------------------------------------------
// cal functions
//-----------------------------------------------------------------------------

unsigned int bwt_find_cals_from_seq_cpu(char *seq, 
					bwt_optarg_t *bwt_optarg, 
					bwt_index_t *index, 
					cal_optarg_t *cal_optarg, 
					array_list_t *cal_list);

unsigned int bwt_find_cals_from_read_cpu(fastq_read_t *read, 
					 bwt_optarg_t *bwt_optarg, 
					 bwt_index_t *index, 
					 cal_optarg_t *cal_optarg, 
					 array_list_t *cal_list);

unsigned int bwt_find_cals_from_seqs_cpu(char **seqs, 
					 unsigned int num_reads,
					 bwt_optarg_t *bwt_optarg, 
					 bwt_index_t *index, 
					 unsigned int *num_unmapped,
					 int *unmapped_array,
					 cal_optarg_t *cal_optarg, 
					 array_list_t *cal_list);

unsigned int bwt_find_cals_from_batch_cpu(fastq_batch_t *batch,
					  bwt_optarg_t *bwt_optarg, 
					  bwt_index_t *index, 
					  fastq_batch_t *unmapped_batch,
					  cal_optarg_t *cal_optarg, 
					  array_list_t *cal_list);

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#endif // BURROWS_WHEELER_H
