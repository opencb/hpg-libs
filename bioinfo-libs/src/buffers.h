#ifndef BUFFERS_H
#define BUFFERS_H

#include "containers/array_list.h"
#include "containers/list.h"

#include "bioformats/bam-sam/alignment.h"
#include "aligners/bwt/bwt.h"

#include "genome.h"
#include "timing.h"
#include "statistics.h"
//#include "bam_commons.h"

//====================================================================================
//  globals
//====================================================================================

#define DNA_FLAG           1
#define RNA_FLAG           2
#define SINGLE_END_FLAG    4
#define PAIRED_END_FLAG    8
#define MATE_PAIR_FLAG    16
#define PAIR1_FLAG        32
#define PAIR2_FLAG        64
#define WRITE_ITEM_FLAG  128
#define SW_ITEM_FLAG     256

//------------------------------------------------------------------------------------

#define SINGLE_END 0
#define PAIRED_END 1
#define MATE_PAIR  2

//------------------------------------------------------------------------------------

#define UNKNOWN_ITEM    0
#define READ_ITEM       1
#define KL_ITEM         2
#define SEED_ITEM       3
#define SPLIT_KL_ITEM   4
#define SW_ITEM         5
#define WRITE_ITEM      6

//------------------------------------------------------------------------------------

#define MISMATCH_FLAG 0
#define MATCH_FLAG    1
#define SPLICE_FLAG   2
//------------------------------------------------------------------------------------

#define NORMAL_MODE 0
#define SEED_MODE   1

//------------------------------------------------------------------------------------

#define INIT_BWT_INDEX         	0
#define INIT_GENOME_INDEX       1
#define FASTQ_READER	 	2
#define BWT_SERVER	   	3
#define REGION_SEEKER		4
#define CAL_SEEKER	   	5
#define SW_SERVER	    	6
#define BATCH_WRITER	 	7
#define FREE_MAIN         	8
#define MAIN_INDEX         	9

//------------------------------------------------------------------------------------

#define FASTQ_READER_ST	 	0
#define BWT_SERVER_ST	   	1
#define REGION_SEEKER_ST	2
#define CAL_SEEKER_ST	   	3
#define SW_SERVER_ST	    	4
#define TOTAL_ST         	5
//------------------------------------------------------------------------------------

#define CHROMOSOME_NUMBER 30

//------------------------------------------------------------------------------------

#define MAX_READ_MAPPINGS 100 

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct write_batch {
  unsigned char flag;
  unsigned int size;
  unsigned int allocated_size;
  void* buffer_p;
} write_batch_t;

//-------------------------------------------------------------------------------------

write_batch_t* write_batch_new(unsigned int size, unsigned char flag);
void write_batch_free(write_batch_t* write_batch_p);

//====================================================================================

typedef struct region_batch {
  array_list_t **allocate_mapping_p;
  fastq_batch_t *unmapped_batch_p;
} region_batch_t;

void region_batch_init(array_list_t **allocate_mapping_p, fastq_batch_t *unmapped_batch_p, 
		       region_batch_t *region_batch_p);
void region_batch_free(region_batch_t *allocate_regions_p);

//====================================================================================

typedef struct sw_batch {
  unsigned int num_reads;
  array_list_t **allocate_cals_p;
  fastq_read_t **allocate_reads_p;
} sw_batch_t;

void sw_batch_init(unsigned int num_reads, array_list_t **allocate_cals_p, 
		   fastq_read_t **allocate_reads_p, sw_batch_t *sw_batch_p);
void sw_batch_free(sw_batch_t *sw_data_p);

//====================================================================================

typedef struct pair_mng {
  int pair_mode;
  size_t min_distance;
  size_t max_distance;
} pair_mng_t;

pair_mng_t *pair_mng_new(int pair_mode, size_t min_distance, size_t max_distance);
void pair_mng_free(pair_mng_t *p);


//=====================================================================================

unsigned int pack_junction(unsigned int chromosome, unsigned int strand, unsigned int start, 
			   unsigned int end, unsigned int junction_id, unsigned int num_reads, 
			   char* buffer_p);


//=====================================================================================
//=====================================================================================

#define UNKNWON_ACTION 0
#define BWT_ACTION     1
#define SEEDING_ACTION 2
#define CAL_ACTION     3
#define PAIR_ACTION    4
#define SW_ACTION      5

//------------------------------------------------------------------------------------

typedef struct aligner_batch {
  int action;
  int all_targets;
  size_t num_targets;
  size_t num_allocated_targets;  
  size_t num_mapping_lists; // = num_reads in fastq_batch

  size_t num_done;
  size_t num_to_do;

  fastq_batch_t *fq_batch;
  size_t *targets;
  array_list_t **mapping_lists;
  char *status;
} aligner_batch_t;

aligner_batch_t *aligner_batch_new(fastq_batch_t *fq_batch);
void aligner_batch_free(aligner_batch_t *p);

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif // BUFFERS_H
