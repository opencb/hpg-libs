#ifndef BWT_H
#define BWT_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <pthread.h> 

#include "commons/string_utils.h"
#include "containers/array_list.h"
#include "containers/linked_list.h"
#include "bioformats/fastq/fastq_read.h"
#include "bioformats/fastq/fastq_batch.h"
#include "bioformats/bam/alignment.h"

#include "BW_io.h"
#include "BW_search.h"
#include "BW_preprocess.h"

#define NONE_HARD_CLIPPING 0
#define START_HARD_CLIPPING 1
#define END_HARD_CLIPPING 2

#define NO_CALS 1
#define EXTRA_CALS 2


#ifndef MAX
  #define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

double global_parallel, global_sequential;

double time_bwt, time_search, time_bwt_seed, time_search_seed;
//-----------------------------------------------------------------------------
// Paratemers for the candidate alignment localizations (CALs)
//-----------------------------------------------------------------------------


typedef struct cal_optarg {
  size_t min_cal_size;
  size_t max_cal_distance;
  size_t min_num_seeds;
  size_t max_num_seeds;
  size_t seed_size;
  size_t min_seed_size;
  size_t num_errors;
} cal_optarg_t;

cal_optarg_t *cal_optarg_new(const size_t min_cal_size, 
			     const size_t max_cal_distance, 
			     const size_t min_num_seeds,
			     const size_t max_num_seeds,
			     const size_t seed_size,
			     const size_t min_seed_size,
			     const size_t num_errors);

void cal_optarg_free(cal_optarg_t *optarg);

//-----------------------------------------------------------------------------

typedef struct cal {
  size_t chromosome_id;
  short int strand;
  size_t start;
  size_t end;
  size_t num_seeds;
  size_t flank_start;
  size_t flank_end;
} cal_t;

cal_t *cal_new(const size_t chromosome_id, 
	       const short int strand,
	       const size_t start, 
	       const size_t end,
	       const size_t num_seeds,
	       const size_t f_start,
	       const size_t f_end);

void cal_free(cal_t *cal);

//-----------------------------------------------------------------------------

typedef struct region {
  size_t chromosome_id;
  short int strand;
  size_t start;
  size_t end;

  size_t seq_start;
  size_t seq_end;
  size_t seq_len;
} region_t;

region_t *region_bwt_new(const size_t chromosome_id, 
	             const short int strand,
	             const size_t start, 
	             const size_t end,
		     const size_t seq_start,
		     const size_t seq_end,
		     const size_t seq_len);

void region_bwt_free(region_t *region);

//-----------------------------------------------------------------------------

typedef struct short_cal {
  size_t start;
  size_t end;

  size_t seq_start;
  size_t seq_end;
  size_t seq_len;

  size_t num_seeds;
} short_cal_t;

short_cal_t *short_cal_new(const size_t start, 
	                   const size_t end,
			   const size_t seq_start,
			   const size_t seq_end,
			   const size_t seq_len,
			   const size_t num_seeds);

void short_cal_free(short_cal_t *short_cal_p);

//-----------------------------------------------------------------------------

typedef struct read_cals {
  fastq_read_t *read;
  array_list_t *cal_list; // array list of cal_t structures
} read_cals_t;

read_cals_t *read_cals_new(fastq_read_t *read);
void read_cals_free(read_cals_t *read_cals);


//-----------------------------------------------------------------------------
// Burrows-Wheeler Transform
//-----------------------------------------------------------------------------

typedef struct bwt_optarg {
  size_t num_errors;
  size_t num_threads;
  size_t max_alignments_per_read;
  size_t report_best;
  size_t report_n_hits;
  char report_all;
} bwt_optarg_t;

bwt_optarg_t *bwt_optarg_new(const size_t num_errors,
			     const size_t num_threads,
			     const size_t max_alginments_per_read,
			     const size_t report_best,
			     const size_t report_n_hits,
			     const size_t report_all
			     );

void bwt_optarg_free(bwt_optarg_t *optarg);

//-----------------------------------------------------------------------------

typedef struct bwt_index {
  comp_matrix h_O, h_rO, h_Oi, h_rOi;
  vector h_C, h_rC, h_C1, h_rC1;
  byte_vector B;
  comp_vector S, Si;
  exome karyotype;
  char *dirname;
} bwt_index_t;

bwt_index_t *bwt_index_new(const char *dirname);
void bwt_index_free(bwt_index_t *index);


void bwt_generate_index_files(char *ref_file, char *output_dir, 
			      unsigned int s_ratio);


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

read_cals_t *read_cals_new(fastq_read_t *read);
void read_cals_free(read_cals_t *read_cals);

//-----------------------------------------------------------------------------
// general functions
//-----------------------------------------------------------------------------

/**
 * @brief  Makes the reverse and complementary from the input sequence.
 * @param  seq input sequence
 * @param  len sequence length input
 * 
 * Makes the reverse and complementary read from input sequence. For it
 * the read input is walked from end to start and all nucleotides are 
 * changed for their complementaries. 
 */
void seq_reverse_complementary(char *seq, unsigned int len);

/**
 */
char* reverse_str(char *src, char *dsp, size_t length);


size_t bwt_map_seq(char *seq, 
		   bwt_optarg_t *bwt_optarg, 
		   bwt_index_t *index, 
		   array_list_t *mapping_list);

size_t bwt_map_read(fastq_read_t *read, 
		    bwt_optarg_t *bwt_optarg, 
		    bwt_index_t *index, 
		    array_list_t *mapping_list);

size_t bwt_map_seqs(char **seqs, 
		    size_t num_reads,
		    bwt_optarg_t *bwt_optarg, 
		    bwt_index_t *index, 
		    char *out_status,
		    array_list_t *mapping_list);

size_t bwt_map_reads(fastq_read_t **reads, 
		     bwt_optarg_t *bwt_optarg, 
		     bwt_index_t *index, 
		     char *out_status,
		     array_list_t *mapping_list);

size_t bwt_map_batch(fastq_batch_t *batch,
		     bwt_optarg_t *bwt_optarg, 
		     bwt_index_t *index, 
		     fastq_batch_t *unmapped_batch,
		     array_list_t *mapping_list);

size_t alignments_filter(char report_all, 
			 size_t report_best, 
			 size_t report_n_hits,
			 array_list_t *mapping_list);

size_t bwt_map_inexact_batch(fastq_batch_t *batch,
			     bwt_optarg_t *bwt_optarg, 
			     bwt_index_t *index, 
			     fastq_batch_t *unmapped_batch,
			     array_list_t *mapping_list);


//-----------------------------------------------------------------------------
// seed functions
//-----------------------------------------------------------------------------

size_t bwt_map_exact_seeds_seq(int padding_left,
			       int padding_right,
			       char *seq, 
			       size_t seed_size, size_t min_seed_size,
			       bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
			       array_list_t *mapping_list, unsigned char step_id);

size_t bwt_map_exact_seeds_seq_by_num(char *seq, 
				      size_t min_num_seeds, size_t max_num_seeds, 
				      size_t seed_size, size_t min_seed_size,
				      bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				      array_list_t *mapping_list);

//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seeds_seq(char *seq, size_t seed_size, size_t min_seed_size,
				 bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				 array_list_t *mapping_list);

//-----------------------------------------------------------------------------
// cal functions
//-----------------------------------------------------------------------------

size_t bwt_find_cals_from_seq(char *seq, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      cal_optarg_t *cal_optarg, 
			      array_list_t *cal_list);

size_t bwt_find_cals_from_read(fastq_read_t *read, 
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       cal_optarg_t *cal_optarg, 
			       array_list_t *cal_list);

size_t bwt_find_cals_from_seqs(char **seqs, 
			       size_t num_reads,
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       cal_optarg_t *cal_optarg, 
			       char *out_status,
			       array_list_t *cal_list);

size_t bwt_find_cals_from_reads(fastq_read_t **reads, 
				bwt_optarg_t *bwt_optarg, 
				bwt_index_t *index, 
				cal_optarg_t *cal_optarg, 
				char *out_status,
				array_list_t *cal_list);

size_t bwt_find_cals_from_batch(fastq_batch_t *batch,
				bwt_optarg_t *bwt_optarg, 
				bwt_index_t *index, 
				fastq_batch_t *unmapped_batch,
				cal_optarg_t *cal_optarg, 
				array_list_t *cal_list);


size_t bwt_generate_cal_list_rna_linkedlist(array_list_t *mapping_list,
					    cal_optarg_t *cal_optarg,
					    array_list_t *cal_list,
					    size_t read_length, size_t nchromosomes);


size_t bwt_generate_cal_list_linkedlist(array_list_t *mapping_list,
					cal_optarg_t *cal_optarg,
					size_t *min_seeds, size_t *max_seeds,
					size_t nchromosomes,
					array_list_t *cal_list);


size_t bwt_generate_cal_list_linked_list_rna(array_list_t *mapping_list,
					     cal_optarg_t *cal_optarg,
					     size_t *min_seeds, size_t *max_seeds,
					     array_list_t *cal_list);


size_t bwt_map_inexact_array_list(array_list_t *reads,
				  bwt_optarg_t *bwt_optarg, 
				  bwt_index_t *index,
				  array_list_t **lists,
				  size_t *num_unmapped, 
				  size_t *unmapped_indices);

void bwt_map_inexact_array_list_by_filter(array_list_t *reads,
					  bwt_optarg_t *bwt_optarg, 
					  bwt_index_t *index,
					  array_list_t **lists,
					  size_t *num_unmapped, 
					  size_t *unmapped_indices);

size_t bwt_generate_cal_list_rna_linked_list(array_list_t *mapping_list,
					     cal_optarg_t *cal_optarg,
					     array_list_t *cal_list,
					     size_t read_length,
					     size_t nchromosomes);

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#endif // BWT_H
