#ifndef BAM_WRITER_H
#define BAM_WRITER_H

#include <stdio.h>
#include <pthread.h>

#include "bam/bam.h"

#include "commons/commons.h"
#include "commons/system_utils.h"

#include "bam_commons.h"
#include "bam_file.h"
#include "chrom_alignments.h"

#define SEQUENTIAL_MODE 0
#define CHROMOSOME_MODE 1

/* **************************************
 *      	Structures    		*
 * *************************************/

/**
* @brief BAM writer server
*
* BAM writer server structure
*/
typedef struct bam_writer {
    int alive;					/**< Flag of thread alive. */
    int mode;					/**< Mode of writing (sequential or by chromosome). */
    int chromosome;				/**< Chromosome to write (can be ALL chromosomes). */
    pthread_mutex_t alive_lock;			/**< Lock for alive variable. */
    pthread_t thread;				/**< Thread structure. */

    bam_file_t* bam_file_p;			/**< BAM file handler. */
    alignments_list_t* alignments_list_p;	/**< List of alignments. */
    int* bam_write_alignment_count;		/**< Number of written alignments. */
} bam_writer_t;

/* **************************************
 *      	Functions    		*
 * *************************************/

/**
*  @brief Creates a new bam writer
*  @param filename BAM file name to write
*  @param alignments_list_p pointer to the alignments_list to store alignments
*  @param bam_header_p pointer to the bam_header of the BAM to write
*  @param mode write mode (sequential or by chromosome)
*  @param chromosome chromosome to write
*  @return pointer to the created bam writer
*  
*  Creates and returns a new bam writer
*/
bam_writer_t* bam_writer_new(char* filename, alignments_list_t* alignments_list_p, bam_header_t* bam_header_p, int mode, int chromosome);

/**
*  @brief Frees a bam writer
*  @param writer_p pointer to the writer
*  @param all flag to free inner structures (0: not free, 1: free)
*  @return void
*  
*  Frees a bam writer and its inner resources
*/
void bam_writer_free(bam_writer_t* writer_p, int all);

/**
*  @brief Starts a bam writer
*  @param writer_p pointer to the writer
*  @return void
*  
*  Starts a bam writer that has been previously created and initialized
*/
void bam_writer_start(bam_writer_t* writer_p);

/**
*  @brief Waits for the bam writer thread until it stops
*  @param writer_p pointer to the writer
*  @return total alignments written
*  
*  Waits for the bam writer thread until it stops
*/
unsigned int bam_writer_join(bam_writer_t* writer_p);

/**
*  @brief Gets the alive state of the writer
*  @param writer_p pointer to the writer
*  @return 0 if dead, 1 if alive
*  
*  Gets the alive state of the writer
*/
int bam_writer_get_alive(bam_writer_t* writer_p);

/**
*  @brief Sets the alive state of the writer
*  @param writer_p pointer to the writer
*  @param alive alive state (0 if dead, 1 if alive)
*  @return void
*  
*  Sets the alive state of the writer
*/
void bam_writer_set_alive(bam_writer_t* writer_p, int alive);

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
void* bam_writer_sequential_thread_function(void* param_p);

/**
*  @brief pthread implementation of the bam reader (by chromosome)
*  @param param_p parameters of the reader
*  @return void*
*  
*  Function with pthread implementation of the bam reader for single chromosome read
*/
void* bam_writer_chromosome_thread_function(void* param_p);

#endif  /* BAM_WRITER_H */
