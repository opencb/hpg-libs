
#ifndef BAM_DATA_BATCH_H
#define BAM_DATA_BATCH_H

#include <pthread.h>
#include <stdint.h>

#include "bam_file.h"
#include "commons.h"

#define MAX_CIGAR_POSITIONS 10


/* **************************************
 *      	Structures    		*
 * *************************************/

/**
* @brief Structure for bam core data
*
* Structure for storing bam core data
*/
typedef struct bam_data_core {
    int map_quality;		/**< Map quality. */
    short int strand;		/**< Strand (1: forward, 0: reverse). */
    short int chromosome;	/**< Chromosome. */
    int start_coordinate;	/**< Start coordinate of the alignment. */
    int alignment_length;	/**< Length of the alignment. */    
    int id_seq_index;		/**< Index of the associated sequence in the batch. */
    int cigar_index;		/**< Index of the associated cigar in the batch. */
    short int paired_end;	/**< Paired end flag (1: paired end, 0: single end). */
} bam_data_core_t;

/**
* @brief Batch of bam data 
*
* Batch of bam data 
* Bam data is a subset of interest of bam1_t structure
*/
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


/* **************************************
 *      	Functions    		*
 * *************************************/

/**
*  @brief Creates a new bam data batch
*  @param num_alignments number of alignments
*  @return pointer to the created bam_data_batch
*  
*  Creates and returns a new bam data batch
*/
bam_data_batch_t* bam_data_batch_new(int num_alignments);

/**
*  @brief Inits a given bam data batch
*  @param bam_data_batch_p pointer to the bam_data_batch
*  @param bam_batch_p pointer to the bam_batch
*  @return pointer to the initialized bam_data_batch
*  
*  Inits a given bam data batch with the information of a bam batch
*/
bam_data_batch_t* bam_data_batch_init(bam_data_batch_t* bam_data_batch_p, bam_batch_t* bam_batch_p);

/**
*  @brief Inits a given bam data batch with data of a single chromosome
*  @param bam_data_batch_p pointer to the bam_data_batch
*  @param bam_batch_p pointer to the bam_batch
*  @param start_alignment index of the first alignment to load
*  @return pointer to the initialized bam_data_batch
*  
*  Inits a given bam data batch with data of a single chromosome from a given start alignment
*/
bam_data_batch_t* bam_data_batch_by_chromosome_init(bam_data_batch_t* bam_data_batch_p, bam_batch_t* bam_batch_p, int start_alignment);

/**
*  @brief Frees a given bam data batch
*  @param bam_data_batch_p pointer to the bam_data_batch
*  @return void
*  
*  Frees a given bam data batch
*/
void bam_data_batch_free(bam_data_batch_t* bam_data_batch_p);

/**
*  @brief Gets bam_data_core from a given bam data batch
*  @param bam_data_batch_p pointer to the bam_data_batch
*  @return pointer to the bam_data_core
*  
*  Gets bam_data_core from a given bam data batch
*/
bam_data_core_t* bam_data_batch_get_core(bam_data_batch_t* bam_data_batch_p);

/**
*  @brief Gets cigar data from a given bam data batch
*  @param bam_data_batch_p pointer to the bam_data_batch
*  @return pointer to the cigar data in uint32_t type
*  
*  Gets cigar data from a given bam data batch
*/
uint32_t* bam_data_batch_get_cigar_data(bam_data_batch_t* bam_data_batch_p);

/**
*  @brief Gets sequence data from a given bam data batch
*  @param bam_data_batch_p pointer to the bam_data_batch
*  @return pointer to the sequence data
*  
*  Gets sequence data from a given bam data batch
*/
char* bam_data_batch_get_id_seq_data(bam_data_batch_t* bam_data_batch_p);

#endif /* BAM_DATA_BATCH_H */


