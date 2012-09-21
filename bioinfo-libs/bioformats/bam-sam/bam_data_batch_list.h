#ifndef BAM_DATA_BATCH_LIST_H
#define BAM_DATA_BATCH_LIST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "bam_data_batch.h"

#define NUM_PROCESSED_TIMES_TO_REMOVE_ITEM  2

//=====================================================
// structures
//=====================================================

// typedef struct bam_data_batch_list_item {
//   int id;
//   int num_alignments;
//   //bam1_t** alignments_p;
//   bam_data_batch_t* batch_p;
//   struct bam_data_batch_list_item* prev_p;
//   struct bam_data_batch_list_item* next_p;
// } bam_data_batch_list_item_t;
//
// typedef struct bam_data_batch_list {
//   int length;
//   int producers;
//   pthread_mutex_t lock;
//   bam_data_batch_list_item_t* first_p;
//   bam_data_batch_list_item_t* last_p;
// } bam_data_batch_list_t;

//=====================================================
// functions
//=====================================================


// bam_data_batch_list_item_t* bam_data_batch_list_item_new(int id, int num_alignments, bam_batch_t* bam_batch_p, bam_data_batch_t* batch_p);
// void bam_data_batch_list_item_free(bam_data_batch_list_item* bam_data_batch_list_item_p, bool all);
//
// void bam_data_batch_list_init(bam_data_batch_list_t* bam_data_batch_list_p, int producers);
// void bam_data_batch_list_insert(bam_data_batch_list_item_t* bam_data_batch_list_item_p, bam_data_batch_list_t* bam_data_batch_list_p);
// bam_data_batch_list_item_t* bam_data_batch_list_remove(bam_data_batch_list_t* bam_data_batch_list_p);
// int bam_data_batch_list_length(bam_data_batch_list_t* bam_data_batch_list_p);
// int bam_data_batch_list_get_producers(bam_data_batch_list_t* bam_data_batch_list_p);
// int bam_data_batch_list_incr_producers(bam_data_batch_list_t* bam_data_batch_list_p);
// int bam_data_batch_list_decr_producers(bam_data_batch_list_t* bam_data_batch_list_p);
// void bam_data_batch_list_print(bam_data_batch_list_t* bam_data_batch_list_p);
//
// void bam_data_batch_list_items_free(bam_data_batch_list_t* list_p);

#endif /* BAM_DATA_BATCH_LIST_H */
