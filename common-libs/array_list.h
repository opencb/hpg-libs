#ifndef ARRAY_LIST_H
#define ARRAY_LIST_H

#include <stdio.h>
#include <pthread.h>

#include <log.h>

//=====================================================
// structures
//=====================================================

typedef struct array_list_item {
  int id;
  int type;
  void* data_p;
} array_list_item_t;

//-----------------------------------------------------

typedef struct array_list {
  size_t capacity;
  size_t size;

  int writers;
  int inserting;
  int removing;

  pthread_mutex_t lock;
  pthread_cond_t condition;

  list_item_t* first_p;
  list_item_t* last_p;
} array_list_t;


/**
 * array_list items functions
 */
list_item_t* array_list_item_new(int id, int type, void* data_p);

void list_item_free(list_item_t* item_p);


/**
 * array_list functions
 */

void array_list_new();

void array_list_init(size_t initial_capacity, array_list_t* array_list_p);

void array_list_free(list_t* list_p, void* (*data_callback) (void* data));

int array_list_insert_item(list_item_t* item_p, list_t* list_p);

list_item_t* array_list_remove_item(list_t* list_p);

int list_insert_item_async(list_item_t* item_p, list_t* list_p);
list_item_t* list_remove_item_async(list_t* list_p);

int array_list_size(array_list_t* list_p);
int array_list_size(array_list_t* list_p);

int list_set_writers(int writers, list_t* list_p);
int list_get_writers(list_t* list_p);
int list_incr_writers(list_t* list_p);
int list_decr_writers(list_t* list_p);

void array_list_print(list_t * list_p);

void **list_to_array(list_t *list_p);


#endif /* ARRAY_LIST_H */
