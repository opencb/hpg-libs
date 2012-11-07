/*
 * linked_list.c
 *
 *  Created on: Nov 7, 2012
 *      Author: imedina
 */

#include "linked_list.h"


linked_list_t* linked_list_new(int SYNC_MODE) {
	linked_list_t *linked_list_p = (linked_list_t*) malloc(sizeof(linked_list_t));
	linked_list_p->size = 0;
	linked_list_p->mode = SYNC_MODE;
	linked_list_p->compare_fn = compare_items;
	linked_list_p->first = NULL;
	linked_list_p->last = NULL;

	pthread_mutex_init(&(linked_list_p->lock), NULL);

	return linked_list_p;
}

void linked_list_free(linked_list_t* linked_list_p, void (*data_callback) (void* data)) {
	if(linked_list_p != NULL) {
		linked_list_clear(linked_list_p, data_callback);
		free(linked_list_p);
	}
}

linked_list_item_t* list_item_new(void *item) {
	linked_list_item_t *linked_list_item_p = (linked_list_item_t*) malloc(sizeof(linked_list_item_t));

	linked_list_item_p->prev = NULL;
	linked_list_item_p->next = NULL;

	return linked_list_item_p;
}



size_t linked_list_size(linked_list_t *linked_list_p) {
	size_t length;

	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		length = linked_list_p->size;

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}
	}else{
		length = ULONG_MAX;
	}

	return length;
}

size_t linked_list_index_of(void *item, linked_list_t *linked_list_p) {
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		size_t i = 0;
		linked_list_item_t *curr_item = linked_list_p->first;
		while(curr_item != NULL) {
			if(linked_list_p->compare_fn(curr_item, item) == 0) {
				return i;
			}
			i++;
		}

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}
	}
	return ULONG_MAX;
}

int linked_list_contains(void *item, linked_list_t *linked_list_p) {
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		linked_list_item_t *curr_item = linked_list_p->first;
		while(curr_item != NULL) {
			if(linked_list_p->compare_fn(curr_item, item) == 0) {
				return 1;
			}
		}

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}
	}
	return 0;
}

int linked_list_clear(linked_list_t *linked_list_p, void (*data_callback) (void* data)) {
	if(linked_list_p != NULL) {
		linked_list_item_t *curr_item = linked_list_p->first;
		linked_list_item_t *next_item;
		while(curr_item != NULL) {
			if(data_callback != NULL) {
				next_item = curr_item->next;
				data_callback(curr_item);
				curr_item = next_item;
			}
		}
		// Set default parameters
		linked_list_p->size = 0;
		linked_list_p->first = NULL;
		linked_list_p->last = NULL;

		return 1;
	}
	return 0;
}



int linked_list_insert(void* item_p, linked_list_t *linked_list_p) {
	return 0;
}

int linked_list_insert_first(void* item_p, linked_list_t *linked_list_p) {
	return 0;
}

int linked_list_insert_last(void* item_p, linked_list_t *linked_list_p) {
	return 0;
}

int linked_list_insert_at(size_t index, void* item_p, linked_list_t *linked_list_p) {
	return 0;
}

int linked_list_insert_all(void** item_p, size_t num_items, linked_list_t *linked_list_p) {
	return 0;
}

int linked_list_insert_all_at(size_t index, void** item_p, size_t num_items, linked_list_t* linked_list_p) {
	return 0;
}



void* linked_list_remove(void *item, linked_list_t *linked_list_p) {
	return NULL;
}

void* linked_list_remove_first(void *item, linked_list_t *linked_list_p) {
	return NULL;
}

void* linked_list_remove_last(void *item, linked_list_t *linked_list_p) {
	return NULL;
}

void* linked_list_remove_at(size_t index, linked_list_t *linked_list_p) {
	return NULL;
}

void** linked_list_remove_range(size_t start, size_t end, linked_list_t* linked_list_p) {
	return NULL;
}



void* linked_list_get(size_t index, linked_list_t *linked_list_p) {
	if(linked_list_p != NULL && index >= 0 && index <= linked_list_p->size) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		size_t i = 0;
		linked_list_item_t *curr_item = linked_list_p->first;
		while(curr_item != NULL && i < index) {
			curr_item = curr_item->next;
			i++;
		}
		return curr_item;

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}

	}

	return NULL;
}

void* linked_list_get_first(linked_list_t *linked_list_p) {
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		return linked_list_p->first;

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}

	}

	return NULL;
}

void* linked_list_get_last(linked_list_t *linked_list_p) {
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		return linked_list_p->last;

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}

	}

	return NULL;
}

linked_list_t* linked_list_sublist(size_t start, size_t end, linked_list_t *linked_list_p, linked_list_t *sublist) {
	return NULL;
}

void* linked_list_set(size_t index, void* new_item, linked_list_t *linked_list_p) {
	return NULL;
}



void linked_list_print(linked_list_t *linked_list_p) {
	printf("[");
	if(linked_list_p != NULL) {
		size_t i = 0;
		linked_list_item_t *curr_item = linked_list_p->first;
		while(curr_item != NULL && i < linked_list_p->size) {
			printf("%s,", (char*)curr_item);
			curr_item = curr_item->next;
			i++;
		}
	}
	printf("]");
}

// void **linked_list_to_array(linked_list_t *linked_list_p) {

static int compare_items(const void *item1, const void *item2) {
	return item1 != item2;
}


int linked_list_swap(const int pos1, const int pos2, linked_list_t *linked_list_p) {
	return 0;
}


void linked_list_set_flag(int flag, linked_list_t *linked_list_p) {

}

int linked_list_get_flag(linked_list_t *linked_list_p) {
	return 0;
}

