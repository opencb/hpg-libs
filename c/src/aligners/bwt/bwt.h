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

#include "search/search.h"
#include "search/preprocess.h"
#include "search/runtime.h"
#include "bwt_commons.h"

#define NONE_HARD_CLIPPING 0
#define START_HARD_CLIPPING 1
#define END_HARD_CLIPPING 2

#define NO_CALS 1
#define EXTRA_CALS 2

#define BACKWARD_ANCHOR 0
#define FORWARD_ANCHOR  1

#define EXTRA_SEED_NONE 0
#define EXTRA_SEED_START 1
#define EXTRA_SEED_END 2

#ifndef MAX
  #define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#define SEED_MISMATCH    0
#define SEED_INSERTION   1
#define SEED_DELETION    2

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

double global_parallel, global_sequential;

double time_bwt, time_search, time_bwt_seed, time_search_seed;
//-----------------------------------------------------------------------------
// Paratemers for the candidate alignment localizations (CALs)
//-----------------------------------------------------------------------------

typedef struct cal_optarg {
  size_t min_cal_size;
  size_t max_cal_distance;
  size_t num_seeds;
  size_t min_num_seeds_in_cal;
  size_t seed_size;
  size_t min_seed_size;
  size_t num_errors;
  size_t max_intron_size;
  size_t min_intron_size;
} cal_optarg_t;

cal_optarg_t *cal_optarg_new(const size_t min_cal_size, 
			     const size_t max_cal_distance, 
			     const size_t num_seeds,
			     const size_t min_num_seeds_in_cal,
			     const size_t seed_size,
			     const size_t min_seed_size,
			     const size_t num_errors, 
			     const size_t max_intron_size,
			     const size_t min_intron_size);

void cal_optarg_free(cal_optarg_t *optarg);

//-----------------------------------------------------------------------------

typedef struct bwt_err {
  int pos;
  char name;
} bwt_err_t;

bwt_err_t *bwt_err_new(int pos, char name);

void bwt_err_free(bwt_err_t *p);

//-----------------------------------------------------------------------------

typedef struct seed_region {
  int read_start;
  int read_end;
  size_t genome_start;
  size_t genome_end;
  int id;
  int fusion_left;
  int fusion_right;

  int strand;
  int chromosome_id;
  int num_mismatches;

  int pos_err;
  int type_err;
  array_list_t *errors_list;

  void *info;
} seed_region_t;

seed_region_t *seed_region_new(int read_start, int read_end, 
			       size_t genome_start, size_t genome_end, int id,
			       int pos_err, int type_err);

void seed_region_free();

void seed_region_simple_free(seed_region_t *seed_region);

//-----------------------------------------------------------------------------

typedef struct cal {
  size_t chromosome_id;
  short int strand;
  size_t start;
  size_t end;
  size_t num_seeds;
  int read_area;
  int l_flank;
  int r_flank;
  int fill_gaps;
  int num_targets;
  int type_seeds;
  int num_mismatches;

  linked_list_t *sr_list;
  linked_list_t *sr_duplicate_list;
  array_list_t *candidates_seeds_start;
  array_list_t *candidates_seeds_end;
  void *info;
} cal_t;

cal_t *cal_new(const size_t chromosome_id,
               const short int strand,
               const size_t start,
               const size_t end,
               const size_t num_seeds,
               linked_list_t *sr_list,
               linked_list_t *sr_duplicate_list);

cal_t *cal_simple_new(const size_t chromosome_id, 
		      const short int strand,
		      const size_t start, 
		      const size_t end);

void cal_free(cal_t *cal);
void cal_simple_free(cal_t *cal);
void cal_print(cal_t *cal);

//-----------------------------------------------------------------------------

typedef struct simple_seed {
  size_t start;
  size_t end;
} simple_seed_t;

simple_seed_t *simple_seed_new(size_t read_start, size_t read_end);
void simple_seed_free(simple_seed_t *simple_seed);

//-----------------------------------------------------------------------------

typedef struct bwt_anchor {
  int strand;
  int chromosome;
  size_t start;
  size_t end;
  int type;
} bwt_anchor_t;

bwt_anchor_t *bwt_anchor_new(int strand, int chromosome, size_t start, size_t end, int type);
void bwt_anchor_free(bwt_anchor_t *bwt_anchor);

//-----------------------------------------------------------------------------

typedef struct region {
  size_t chromosome_id;
  short int strand;
  size_t start;
  size_t end;
  size_t seq_start;
  size_t seq_end;
  size_t seq_len;
  int id;

  int pos_err;
  int type_err;
} region_t;

region_t *region_bwt_new(const size_t chromosome_id, 
			 const short int strand,
			 const size_t start, 
			 const size_t end,
			 const size_t seq_start,
			 const size_t seq_end,
			 const size_t seq_len,
			 const int id);

void region_bwt_free(region_t *region);

//-----------------------------------------------------------------------------

typedef struct short_cal {
  size_t start;
  size_t end;
  size_t seq_len;
  size_t num_seeds;
  size_t seq_start;
  size_t seq_end;

  int pos_err;
  int type_err;
  
  linked_list_t *sr_list;
  linked_list_t *sr_duplicate_list;
  unsigned char *seeds_ids_array;
} short_cal_t;

short_cal_t *short_cal_new(const size_t start, 
	                   const size_t end,
			   const size_t seq_start,
			   const size_t seq_end,
			   const size_t seq_len,
			   const int max_seeds,
			   const int id,
			   const int pos_err,
			   const int pos_type);

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
  int filter_read_mappings;
  int filter_seed_mappings;
} bwt_optarg_t;

bwt_optarg_t *bwt_optarg_new(const size_t num_errors,
			     const size_t num_threads,
			     const int filter_read_mappings, 
			     const int filter_seed_mappings);

void bwt_optarg_free(bwt_optarg_t *optarg);


//-----------------------------------------------------------------------------

typedef struct {
  exome *karyotype;
  bwt_config_t bwt_config;
  bwt_index *backward;
  bwt_index *forward;
  bwt_index *backward_rev;
  bwt_index *forward_rev;
  char *dirname;
} bwt_index_t;

bwt_index_t *bwt_index_new(const char *dirname, bool inverse_sa);

void bwt_index_free(bwt_index_t *index);

void bwt_generate_index_files(char *ref_file, char *output_dir, 
			      unsigned int s_ratio, bool duplicate_strand,
			      char *bases);

//void bwt_generate_index_files_bs(char *ref_file, char *output_dir, 
//				 unsigned int s_ratio, char *bases);

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

alignment_t* add_optional_fields(alignment_t *alignment, size_t n_mappings, int read_length);

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

size_t bwt_map_inexact_batch(fastq_batch_t *batch,
			     bwt_optarg_t *bwt_optarg, 
			     bwt_index_t *index, 
			     fastq_batch_t *unmapped_batch,
			     array_list_t *mapping_list);


size_t bwt_map_inexact_read(fastq_read_t *read, 
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    array_list_t *mapping_list);

size_t bwt_map_inexact_read_bs(fastq_read_t *read, 
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       array_list_t *mapping_list, 
			       int type);

//-----------------------------------------------------------------------------
// seed functions
//-----------------------------------------------------------------------------

size_t bwt_map_seeds_IA(int padding_left,
			int padding_right,
			char *seq, 
			size_t seed_size, size_t min_seed_size,
			bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
			array_list_t *mapping_list, unsigned char step_id);

size_t bwt_map_exact_seeds_seq(int padding_left,
			       int padding_right,
			       char *seq, 
			       size_t seed_size, size_t min_seed_size,
			       bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
			       array_list_t *mapping_list, unsigned char step_id);

size_t bwt_map_exact_seeds_seq_by_num(char *seq, size_t num_seeds,
				      size_t seed_size, size_t min_seed_size,
				      bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				      array_list_t *mapping_list);

size_t bwt_map_exact_seeds_between_coords(int start_position, int end_position, 
					  char *seq, int seed_size, int min_seed_size,
					  bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
					  array_list_t *mapping_list, int extra_seed,
					  int *last_seed_id);

size_t bwt_map_exact_seeds_seq_by_num_bs(char *seq, size_t num_seeds,
					 size_t seed_size, size_t min_seed_size,
					 bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
					 array_list_t *mapping_list);

//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seeds_seq(char *seq, size_t seed_size, size_t min_seed_size,
				 bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				 array_list_t *mapping_list);


size_t bwt_map_inexact_seeds_by_region(int start_position, int end_position, 
				       int strand, 
				       int chromosome_target, int start_target,
				       int end_target,
				       char *seq, size_t seed_size,
				       size_t min_seed_size,
				       bwt_optarg_t *bwt_optarg,
				       bwt_index_t *index, 
				       array_list_t *mapping_list);

size_t bwt_map_exact_seeds_by_region(int start_position, int end_position, 
				     char *seq, int seed_size, int min_seed_size,
				     bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				     array_list_t *mapping_list);

void insert_seeds_and_merge(array_list_t *mapping_list, linked_list_t ***cals_list,  size_t max_cal_distance);

//-----------------------------------------------------------------------------
// cal functions
//-----------------------------------------------------------------------------

size_t bwt_generate_cals(char *seq, size_t seed_size, 
			 bwt_optarg_t *bwt_optarg, 
			 cal_optarg_t *cal_optarg,			
			 bwt_index_t *index, 
			 array_list_t *cal_list, 
			 unsigned int nchromosomes);


size_t bwt_generate_cals_bs(char *seq, char *seq2, size_t seed_size, bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, bwt_index_t *index2, array_list_t *cal_list);

size_t bwt_generate_cals_between_coords(int strand_target, int chromosome_target,
					size_t start_target, size_t end_target, 
					int start_position, int end_posistion, 
					char *seq, int seed_size, int min_seed_size,
					bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
					array_list_t *init_list, array_list_t *cal_list);

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


/*size_t bwt_generate_cal_list_rna_linkedlist(array_list_t *mapping_list,
					    cal_optarg_t *cal_optarg,
					    array_list_t *cal_list,
					    size_t read_length, size_t nchromosomes);

*/

size_t bwt_generate_cal_list_linked_list(array_list_t *mapping_list,
					 cal_optarg_t *cal_optarg,
					 int *min_seeds, int *max_seeds,
					 size_t nchromosomes,
					 array_list_t *cal_list,
					 size_t read_length,
					 size_t min_cal_size,
					 int type_seeds);


/*size_t bwt_generate_cal_list_linked_list_rna(array_list_t *mapping_list,
					     cal_optarg_t *cal_optarg,
					     size_t *min_seeds, size_t *max_seeds,
					     array_list_t *cal_list);

*/

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

void bwt_map_inexact_array_list_by_filter_bs(array_list_t *reads,
					     bwt_optarg_t *bwt_optarg, 
					     bwt_index_t *index,
					     array_list_t **lists,
					     size_t *num_unmapped, 
					     size_t *unmapped_indices);

size_t bwt_map_forward_inexact_seq(char *seq, 
				   bwt_optarg_t *bwt_optarg, 
				   bwt_index_t *index, 
				   array_list_t *mapping_list);
  
size_t bwt_generate_cal_list_rna_linked_list(array_list_t *mapping_list,
					     cal_optarg_t *cal_optarg,
					     array_list_t *cal_list,
					     size_t read_length,
					     size_t nchromosomes);


size_t bwt_generate_cal_rna_list_linked_list(array_list_t *mapping_list,
                                             cal_optarg_t *cal_optarg,
                                             size_t *min_seeds, int *max_seeds,
                                             size_t nchromosomes,
                                             array_list_t *cal_list,
                                             size_t read_length);

size_t bwt_map_inexact_read_2(fastq_read_t *read, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list);

/*void append_seed_region_linked_list(linked_list_t* sr_list,
				    size_t read_start, size_t read_end, 
				    size_t genome_start, size_t genome_end, 
				    int seed_id);
*/

//-----------------------------------------------------------------------------
void bwt_cigar_cpy(alignment_t *mapping, char *quality);

/*void initReplaceTable_bs(const char *str);

char * readNucleotide(const char *directory, const char *name);

void saveNucleotide(char *nucleotide, const char *directory, const char *name);
*/
//-----------------------------------------------------------------------------

#endif // BWT_H
