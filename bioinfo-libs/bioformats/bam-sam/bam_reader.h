
#ifndef BAM_READER_H
#define BAM_READER_H

#include <stdio.h>
#include <pthread.h>

#include "bam_file.h"
#include "chrom_alignments.h"
#include "commons.h"
#include "list.h"
#include "system_utils.h"

#define SEQUENTIAL_MODE  0
#define CHROMOSOME_MODE  1
#define LIST_INSERT_MODE 2

#define NO_SORT 0
#define SORT_BY_POSITION  1
#define SORT_BY_ID   2

#define MIN_FREE_MEMORY_ALLOWED_FOR_BAM_READ_KB  200

/* **************************************
 *      	Structures    		*
 * *************************************/

/**
* @brief BAM reader server
*
* BAM reader server structure
*/
typedef struct bam_reader {
    int alive;					/**< Flag of thread alive. */
    int mode;					/**< Mode of reading (sequential, by chromosome, inserting in batch list). */
    int sort;					/**< Flag for sorting. */
    int chromosome;				/**< Chromosome to read (can be ALL chromosomes). */
    int max_estimated_alignments;		/**< Max. estimated alignments. */
    size_t batch_size;				/**< Size of the batch. */
    int base_quality;				/**< Base quality for quality normalization. */
    pthread_mutex_t alive_lock;			/**< Lock for alive variable. */
    pthread_t thread;				/**< Thread structure. */

    bam_file_t* bam_file_p;			/**< BAM file handler. */
    alignments_list_t* alignments_list_p;	/**< List of alignments. */
    list_t* bam_batch_list_p;			/**< BAM batch list. */
    list_t* bam_data_batch_list_p;		/**< BAM data batch list. */
} bam_reader_t;

/* **************************************
 *      	Functions    		*
 * *************************************/

/**
*  @brief Creates a new bam reader
*  @param filename BAM file name to read
*  @param batch_size size of the batch
*  @param base_quality base quality for quality values normalization
*  @param alignments_list_p pointer to the alignments_list to store alignments
*  @param mode sequential, by chromosome or insert-in-list mode
*  @param sort flag for sorting
*  @param chromosome chromosome to read
*  @return pointer to the created bam reader
*  
*  Creates and returns a new bam reader that stores alignments in an alignment list
*/
bam_reader_t* bam_reader_new(char* filename, size_t batch_size, int base_quality, alignments_list_t* alignments_list_p, int mode, int sort, int chromosome);

/**
*  @brief Creates a new bam reader
*  @param filename BAM file name to read
*  @param batch_size size of the batch
*  @param base_quality base quality for quality values normalization
*  @param bam_data_batch_list_p pointer to the bam_data_batch_list to store alignments
*  @param mode sequential, by chromosome or insert-in-list mode
*  @return pointer to the created bam reader
*  
*  Creates and returns a new bam reader that stores alignments by batch in a bam_data_batch_list
*/
bam_reader_t* bam_reader_by_batch_new(char * filename, size_t batch_size, int base_quality, list_t* bam_data_batch_list_p, int mode);

/**
*  @brief Frees a bam reader
*  @param reader_p pointer to the reader
*  @return void
*  
*  Frees a bam reader and its inner resources
*/
void bam_reader_free(bam_reader_t* reader_p);

/**
*  @brief Starts a bam reader
*  @param reader_p pointer to the reader
*  @return void
*  
*  Starts a bam reader that has been previously created and initialized
*/
void bam_reader_start(bam_reader_t* reader_p);

/**
*  @brief Waits for the bam reader thread until it stops
*  @param reader_p pointer to the reader
*  @return total alignments read
*  
*  Waits for the bam reader thread until it stops
*/
unsigned int bam_reader_join(bam_reader_t* reader_p);

/**
*  @brief Gets the alive state of the reader
*  @param reader_p pointer to the reader
*  @return 0 if dead, 1 if alive
*  
*  Gets the alive state of the reader
*/
int bam_reader_get_alive(bam_reader_t* reader_p);

/**
*  @brief Sets the alive state of the reader
*  @param reader_p pointer to the reader
*  @param alive alive state (0 if dead, 1 if alive)
*  @return void
*  
*  Sets the alive state of the reader
*/
void bam_reader_set_alive(bam_reader_t* reader_p, int alive);

/* **********************************************
 *      	Thread functions    		*
 * *********************************************/

/**
*  @brief pthread implementation of the bam reader (sequential)
*  @param param_p parameters of the reader
*  @return void*
*  
*  Function with pthread implementation of the bam reader for sequential read
*/
void* bam_reader_sequential_thread_function(void* param_p);

/**
*  @brief pthread implementation of the bam reader (by chromosome)
*  @param param_p parameters of the reader
*  @return void*
*  
*  Function with pthread implementation of the bam reader for single chromosome read
*/
void* bam_reader_chromosome_thread_function(void* param_p);

/**
*  @brief pthread implementation of the bam reader (list insert)
*  @param param_p parameters of the reader
*  @return void*
*  
*  Function with pthread implementation of the bam reader that stores bam data batches in a list
*/
void* bam_reader_list_insert_thread_function(void* param_p);

/**
*  @brief pthread implementation of the bam reader (list insert, by chromosome)
*  @param param_p parameters of the reader
*  @return void*
*  
*  Function with pthread implementation of the bam reader that stores bam data batches in a list 
*  selecting only content from one chromosome
*/
void* bam_reader_list_insert_by_chromosome_thread_function(void* param_p);

#endif
