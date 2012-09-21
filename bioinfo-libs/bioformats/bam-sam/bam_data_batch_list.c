#include "bam_data_batch_list.h"

//=====================================================
// functions to manage bam data batch list item
//=====================================================


// bam_data_batch_list_item_t* bam_data_batch_list_item_new(int id, int num_alignments, bam_batch_t* bam_batch_p, bam_data_batch_t* bam_data_batch_p) {
//   bam_data_batch_list_item_t* item_p;
//
//   item_p = (bam_data_batch_list_item_t*) malloc(sizeof(bam_data_batch_list_item_t));
//   item_p->id = id;
//   item_p->num_alignments = num_alignments;
//   //item_p->alignments_p = bam_batch_p->alignments_p;
//   item_p->batch_p = bam_data_batch_p;
//
//   return item_p;
// }

//-----------------------------------------------------
// bam_data_batch_list_item_free
//
// free memory for this structure,
// the input parameter 'all' indicates free memory
// associated to the pointers in his structure
//-----------------------------------------------------

// void bam_data_batch_list_item_free(bam_data_batch_list_item_t* item_p, bool all) {
//
//   if (all) {
//       if (item_p->batch_p != NULL) {
//       bam_data_batch_free(item_p->batch_p);
//     }
//   }
//
//   free(item_p);
// }

//=====================================================
// functions to manage bam data batch list
//=====================================================

//-----------------------------------------------------
// bam_data_batch_list_init
//
// init to zero the object and initialize the lock
//-----------------------------------------------------

// void bam_data_batch_list_init(bam_data_batch_list_t* list_p, int producers) {
//    memset(list_p, 0, sizeof(bam_data_batch_list_t));
//    list_p->lock = PTHREAD_MUTEX_INITIALIZER;
//    list_p->producers = producers;
// }

//-----------------------------------------------------
// bam_data_batch_list_insert
//
// insert a batch object in the end of the list,
// according to fifo order,
//-----------------------------------------------------

// void bam_data_batch_list_insert(bam_data_batch_list_item_t* item_p, bam_data_batch_list_t* list_p) {
//
//   if (list_p==NULL) return;
//
//   pthread_mutex_lock(&list_p->lock);
//
//   if (list_p->first_p==NULL) {
//     //item_p->id = 0;
//     item_p->prev_p = NULL;
//     item_p->next_p = NULL;
//     list_p->first_p = item_p;
//     list_p->last_p = item_p;
//   } else {
//     list_p->last_p->next_p = item_p;
//     //item_p->id = list_p->last_p->id + 1;
//     item_p->prev_p = list_p->last_p;
//     item_p->next_p = NULL;
//     list_p->last_p = item_p;
//   }
//
//   list_p->length++;
//
//   pthread_mutex_unlock(&list_p->lock);
// }

//-----------------------------------------------------
// bam_data_batch_list_remove
//
// remove the first batch object in the begining
// of the list,
// according to fifo order,
//-----------------------------------------------------

// bam_data_batch_list_item_t* bam_data_batch_list_remove(bam_data_batch_list_t* list_p) {
//
//   if (list_p==NULL) return NULL;
//
//   pthread_mutex_lock(&list_p->lock);
//
//   // just get the first element, and if is not null
//   // update the first-element pointer
//   //
//   bam_data_batch_list_item_t* item_p = list_p->first_p;
//
//   if (item_p!=NULL) {
//     list_p->first_p = item_p->next_p;
//     list_p->length--;
//   }
//
//   pthread_mutex_unlock(&list_p->lock);
//
//   return item_p;
// }

//-----------------------------------------------------
// bam_data_batch_list_length
//
// returns list length
//-----------------------------------------------------

// int bam_data_batch_list_length(bam_data_batch_list_t* list_p) {
//
//   int length = 0;
//
//   if (list_p==NULL) return length;
//
//   pthread_mutex_lock(&list_p->lock);
//   length = list_p->length;
//   pthread_mutex_unlock(&list_p->lock);
//
//   return length;
// }

//-----------------------------------------------------
// bam_data_batch_list_get_producers
//
// insert a batch object in the end of the list,
// according to fifo order,
//-----------------------------------------------------

// int bam_data_batch_list_get_producers(bam_data_batch_list_t* list_p) {
//
//   int producers = 0;
//
//   if (list_p==NULL) return producers;
//
//   pthread_mutex_lock(&list_p->lock);
//   producers = list_p->producers;
//   pthread_mutex_unlock(&list_p->lock);
//
//   return producers;
// }

//-----------------------------------------------------
// bam_data_batch_list_incr_producers
//
// insert a batch object in the end of the list,
// according to fifo order,
//-----------------------------------------------------

// int bam_data_batch_list_incr_producers(bam_data_batch_list_t* list_p) {
//
//   if (list_p==NULL) return 0;
//
//   pthread_mutex_lock(&list_p->lock);
//   list_p->producers++;
//   pthread_mutex_unlock(&list_p->lock);
//
//   return list_p->producers;
// }

//-----------------------------------------------------
// bam_data_batch_list_decr_producers
//
// decrement in 1 the number of producers in the list
//
//-----------------------------------------------------

// int bam_data_batch_list_decr_producers(bam_data_batch_list_t* list_p) {
//
//   if (list_p==NULL) return 0;
//
//   pthread_mutex_lock(&list_p->lock);
//   list_p->producers--;
//   pthread_mutex_unlock(&list_p->lock);
//
//   return list_p->producers;
// }

//-----------------------------------------------------
// print batch list
//
//
//-----------------------------------------------------

// void bam_data_batch_list_print(bam_data_batch_list_t* list_p) {
//   int i;
//
//   if (list_p==NULL) return;
//
//   pthread_mutex_lock(&list_p->lock);
//
//   printf("Number of items: %i\n", list_p->length);
//
//   bam_data_batch_list_item_t* item_p = list_p->first_p;
//
//   while(item_p!=NULL) {
//     printf("batch id: %i", item_p->id);
//     printf("\talignments: %i\n", item_p->num_alignments);
//
//     for(int i=0; i<item_p->num_alignments; i++) {
//  //printf("i: %d, strand: %i\n", i, item_p->batch_p[i].strand);
//  //printf("i: %d, sequence: %i\n", i, item_p->sequence[i]);
//  //printf("i: %d, cigar: %i\n", i, item_p->cigar[i]);
//     }
//
//     item_p = item_p->next_p;
//   }
//
//   pthread_mutex_unlock(&list_p->lock);
// }

//-----------------------------------------------------
// bam_data_batch_list_items_free
//
// free list itmes
//
//-----------------------------------------------------

// void bam_data_batch_list_items_free(bam_data_batch_list_t* list_p) {
//   bam_data_batch_list_item_t* item_p;
//
//   while ((item_p = bam_data_batch_list_remove(list_p)) != NULL) {
//     bam_data_batch_list_item_free(item_p, true);
//   }
// }

