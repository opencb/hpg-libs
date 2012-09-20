#ifndef FASTQ_READ_H
#define FASTQ_READ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ID_MAX_BYTES		(1 * 1024 * 1024)
#define INDEX_MAX		100000
#define SEQUENCE_MAX_BYTES	(2 * 1024 * 1024)

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Fastq read structure 
* 
* Structure for storing a fastq read
*/
typedef struct fastq_read {
    char *id;			/**< Id of the read. */
    char *sequence;		/**< Sequence of nts. */
    char *quality;		/**< Qualities. */
} fastq_read_t;

//--------------------------------------------

fastq_read_t *fastq_read_new(char *id, char *sequence, char *quality);
void fastq_read_free(fastq_read_t *fq_read);

//--------------------------------------------
//--------------------------------------------

/**
* @brief Fastq read statistics
* 
* Structure for storing statistics from a fastq read
*/
typedef struct fastq_read_stats {
    char *id;			/**< Id of the read. */
    short int seq_length;	/**< Sequence length. */
    float quality_average;	/**< Mean quality average. */
} fastq_read_stats_t;

/* **************************************
 *  		Functions		*
 * *************************************/

float fastq_quality_average(fastq_read_t fq_read_t);
void fastq_nt_quality_average(fastq_read_t fq_read_t);
void fastq_nt_count(fastq_read_t fq_read_t);
void fastq_length(fastq_read_t fq_read_t);
void fastq_remove_barcode(fastq_read_t fq_read_t, char *barcode);
void fastq_stats(fastq_read_t fq_read_t, fastq_read_stats_t *read_stats_t);

#endif	/*  FASTQ_READ_H  */
