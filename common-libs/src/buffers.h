#ifndef BUFFERS_H
#define BUFFERS_H

#include "containers/array_list.h"
#include "containers/list.h"

#include "bioformats/bam-sam/alignment.h"
#include "aligners/bwt/bwt.h"
#include "aligners/bwt/genome.h"
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

#define SINGLE_END_MODE 0
#define PAIRED_END_MODE 1
#define MATE_PAIR_MODE  2

//------------------------------------------------------------------------------------

#define PAIR_1 1
#define PAIR_2 2
#define PAIR_1_2 3

//------------------------------------------------------------------------------------

#define UNKNOWN_ITEM    0
#define READ_ITEM       1
#define KL_ITEM         2
#define SEED_ITEM       3
#define SPLIT_KL_ITEM   4
#define SW_ITEM         5
#define WRITE_ITEM      6

//------------------------------------------------------------------------------------

#define MISMATCH_FLAG        0
#define MATCH_FLAG           1
#define SPLICE_EXACT_FLAG    2
#define SPLICE_EXTEND_FLAG   3

//------------------------------------------------------------------------------------

#define NORMAL_MODE 0
#define SEED_MODE   1

//------------------------------------------------------------------------------------

#define INIT_BWT_INDEX         	0
#define BWT_SERVER	   	1
#define REGION_SEEKER		2
#define CAL_SEEKER	   	3
#define SW_SERVER	    	4
#define FREE_MAIN         	5
#define MAIN_INDEX         	6

//------------------------------------------------------------------------------------

#define BWT_SERVER_ST	   	0
#define REGION_SEEKER_ST	1
#define CAL_SEEKER_ST	   	2
#define SW_SERVER_ST	    	3
#define TOTAL_ST         	4
//------------------------------------------------------------------------------------

#define CHROMOSOME_NUMBER 30

//------------------------------------------------------------------------------------

#define MAX_READ_MAPPINGS 100 

//====================================================================================
//  structures and prototypes
//====================================================================================

/**
 * @brief Structure for store in files all the process data.
 * 
 * Structure for store in files all data process by each pipeline phase.
 */
typedef struct write_batch {
  unsigned char flag;           /**< type of data stored*/
  unsigned int size;            /**< actual size of the batch (in bytes)*/
  unsigned int allocated_size;  /**< maximum size of the batch (in bytes)*/
  void* buffer_p;               /**< buffer to store data*/
} write_batch_t;

//-------------------------------------------------------------------------------------

/**
 * @brief  Constructor for the @a write_batch_t structure.
 * @param  allocate_size maximum size of the batch (in bytes)
 * @param  flag type of data stored
 * @return Pointer to the new structure.
 * 
 * Constructor function that allocates memory for
 * the batch_writer structure and initialize it.
 */
write_batch_t* write_batch_new(unsigned int allocate_size, unsigned char flag);

/**
 * @brief  Destrcutor for the @a write_batch_t structure.
 * @param  write_batch_p[out] pointer to the structure to free
 * 
 * @a write_batch_t destructor that frees the memory previosly
 * allocated by constructor @a write_batch_new
 */
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
/*
typedef struct sw_batch {
  unsigned int num_reads;
  array_list_t **allocate_cals_p;
  fastq_read_t **allocate_reads_p;
} sw_batch_t;

void sw_batch_init(unsigned int num_reads, array_list_t **allocate_cals_p, 
		   fastq_read_t **allocate_reads_p, sw_batch_t *sw_batch_p);
void sw_batch_free(sw_batch_t *sw_data_p);
*/
//====================================================================================

typedef struct pair_mng {
  int pair_mode;
  size_t min_distance;
  size_t max_distance;
} pair_mng_t;

pair_mng_t *pair_mng_new(int pair_mode, size_t min_distance, size_t max_distance);
void pair_mng_free(pair_mng_t *p);

//=====================================================================================
//=====================================================================================

#define UNKNWON_ACTION 0
#define BWT_ACTION     1
#define SEEDING_ACTION 2
#define CAL_ACTION     3
#define PAIR_ACTION    4
#define SW_ACTION      5

//------------------------------------------------------------------------------------

typedef struct mapping_batch {
  int action;
  size_t num_targets;
  size_t num_allocated_targets;

  size_t num_to_do;

  array_list_t *fq_batch;
  size_t *targets;
  array_list_t **mapping_lists;
  char *status; //TODO: ¿?
  pair_mng_t *pair_mng;
} mapping_batch_t;

mapping_batch_t *mapping_batch_new(array_list_t *fq_batch, pair_mng_t *pair_mng);
void mapping_batch_free(mapping_batch_t *p);

//=====================================================================================
//=====================================================================================

/**
 * @brief Structure for store all process data in @a region_seeker_server.
 * 
 * Structure for store process data in @a region_seeker_server, contains an 
 * array list with all mappings found for all seeds in each read and a
 * batch with all reads unmapped.
 */
typedef struct cal_batch {
  array_list_t **allocate_mapping;  /**< array list with all mappings found for all seeds in each read*/
  fastq_batch_t *unmapped_batch;   /**< batch with all reads unmapped*/
} cal_batch_t;

//------------------------------------------------------------------------------------

/**
 * @brief  Constructor for the @a cal_batch_t structure.
 * @param  allocate_mapping_p array list with all mappings found for all seeds in each read
 * @param  unmapped_batch_p batch with all reads unmapped
 * @return Pointer to the new structure
 * 
 * Constructor function that allocates memory for
 * the cal_batch_t structure and initialize it.
 */
cal_batch_t* cal_batch_new(array_list_t **allocate_mapping, fastq_batch_t *unmapped_batch);

/**
 * @brief  Destrcutor for the @a cal_batch_t structure.
 * @param  allocate_cal_p[out] pointer to the structure to free
 * 
 * @a cal_batch_t destructor that frees the memory previosly
 * allocated by constructor @a cal_batch_new
 */
void cal_batch_free(cal_batch_t *allocate_cal);

//====================================================================================

/**
 * @brief Structure for store all process data in @a cal_seeker_server.
 * 
 * Structure for store process data in @a cal_seeker_server, store all cals found
 * for each read unmapped.
 */
typedef struct sw_batch {
  unsigned int num_reads;          /**< number of reads allocate in the batch*/
  array_list_t **allocate_cals_p;  /**< array list that store all cals found for each read */
  fastq_read_t **allocate_reads_p; /**< array that store for each read the header, the sequence and the quality*/
}sw_batch_t;

/**
 * @brief  Constructor for the @a sw_batch_t structure.
 * @param  num_reads number of reads that will be store
 * @param  allocate_cals_p array list that store all cals found for each read
 * @param  allocate_reads_p array that store for each read the header, the sequence and the quality
 * @return Pointer to the new structure
 * 
 * Constructor function that allocates memory for
 * the sw_batch_t structure and initialize it.
 */
sw_batch_t* sw_batch_new(unsigned int num_reads, array_list_t **allocate_cals_p, fastq_read_t **allocate_reads_p);

/**
 * @brief  Destrcutor for the @a sw_batch_t structure.
 * @param  sw_batch_p[out] pointer to the structure to free
 * 
 * @a sw_batch_t destructor that frees the memory previosly
 * allocated by constructor @a sw_batch_new
 */
void sw_batch_free(sw_batch_t *sw_batch_p);

//=====================================================================================

/**
 * @brief  Store in @a buffer_p all data information for one splice junction.
 * @param  chromosome chromosome where the splice junction was found
 * @param  strand strand where the splice junction was found
 * @param  start splice junction start coordinate 
 * @param  end splice junction end coordinate
 * @param  junction_id splice junction identification
 * @param  num_reads number of reads that support the splice junction
 * @param  buffer_p buffer to store splice junction data information  
 * @return Number of bytes stored in @a buffer_p
 * 
 * Store all data information for each splice junction found. It allocates all data 
 * in one buffer.
 */
unsigned int pack_junction(unsigned int chromosome, unsigned int strand, 
			   size_t start, size_t end, 
			   size_t junction_id, size_t num_reads, 
			   char* buffer_p);

//======================================================================================

#endif // BUFFERS_H
