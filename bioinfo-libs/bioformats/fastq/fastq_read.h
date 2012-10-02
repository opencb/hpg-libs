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

    int length;
} fastq_read_t;

/**
* @brief Fastq read statistics
* 
* Structure for storing statistics from a fastq read
*/
typedef struct fastq_read_stats {
    int length;
    float quality_average;	/**< Mean quality average. */
	int Ns;
	int N_out_quality;

	int num_A;
	int num_T;
	int num_C;
	int num_G;
	int num_N;
} fastq_read_stats_t;


/* **************************************
 *  		Functions		*
 * *************************************/


fastq_read_t *fastq_read_new(char *id, char *sequence, char *quality);

void fastq_read_free(fastq_read_t *fq_read);



fastq_read_stats_t *fastq_read_stats_new();

void fastq_read_stats_init(int length, float qual_average, int Ns, int N_out_qual, fastq_read_stats_t *read_stats);

void fastq_read_stats_free(fastq_read_stats_t *read_stats);

int fastq_read_stats(fastq_read_t *fq_read, fastq_read_stats_t *read_stats);

int fastq_reads_stats(array_list_t *fq_reads, array_list_t *reads_stats);




float fastq_quality_average(fastq_read_t fq_read_t);
void fastq_nt_quality_average(fastq_read_t fq_read_t);
void fastq_nt_count(fastq_read_t fq_read_t);
void fastq_length(fastq_read_t fq_read_t);
void fastq_remove_barcode(fastq_read_t fq_read_t, char *barcode);

#endif	/*  FASTQ_READ_H  */
