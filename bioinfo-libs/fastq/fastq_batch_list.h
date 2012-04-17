
#ifndef _FASTQ_BATCH_LIST_H
#define _FASTQ_BATCH_LIST_H

#include <stdio.h>
#include <pthread.h>

#include "fastq_batch.h"
#include "fastq_file.h"

//=====================================================
// structures
//=====================================================

typedef struct fastq_batch_list_item {
	int id;
	fastq_batch_t *batch_p;
	struct fastq_batch_list_item *prev_p;
	struct fastq_batch_list_item *next_p;
} fastq_batch_list_item_t;


typedef struct fastq_batch_list {
	int length;
	int length_by_source_id[MAX_NUM_PRODUCERS];
	int producers;
	fastq_batch_list_item_t *first_p;
	fastq_batch_list_item_t *last_p;
	pthread_mutex_t lock;
} fastq_batch_list_t;


//=====================================================
// functions
//=====================================================

void fastq_batch_list_item_free(fastq_batch_list_item_t* fastq_batch_list_item_p, int all);

void fastq_batch_list_init(fastq_batch_list_t* fastq_batch_list_p, int producers);

void fastq_batch_list_insert(fastq_batch_list_item_t* fastq_batch_list_item_p, fastq_batch_list_t* fastq_batch_list_p);

fastq_batch_list_item_t* fastq_batch_list_remove(fastq_batch_list_t* fastq_batch_list_p);

int fastq_batch_list_length(fastq_batch_list_t* fastq_batch_list_p);

int fastq_batch_list_length_by_source_id(fastq_batch_list_t* list_p, int source_id);

int fastq_batch_list_get_producers(fastq_batch_list_t* fastq_batch_list_p);

int fastq_batch_list_incr_producers(fastq_batch_list_t* fastq_batch_list_p);

int fastq_batch_list_decr_producers(fastq_batch_list_t* fastq_batch_list_p);

void fastq_batch_list_print(fastq_batch_list_t* fastq_batch_list_p);

void fastq_batch_list_items_free(fastq_batch_list_t* list_p);

#endif /* FASTQ_BATCH_LIST_H */
