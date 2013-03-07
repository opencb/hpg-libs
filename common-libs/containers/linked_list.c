/*
 * linked_list.c
 *
 *  Created on: Nov 7, 2012
 *      Author: imedina, hmartinez
 */

#include "linked_list.h"

void print_item(void *item) {
  printf("%d->", (int)item);
}

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

linked_list_item_t* linked_list_item_new(void *item) {
	linked_list_item_t *linked_list_item_p = (linked_list_item_t*) malloc(sizeof(linked_list_item_t));

	linked_list_item_p->prev = NULL;
	linked_list_item_p->next = NULL;
	linked_list_item_p->item = item;

	return linked_list_item_p;
}

void linked_list_item_free(linked_list_item_t *linked_list_item, void (*data_callback) (void* data)) {
	if(data_callback) {
		data_callback(linked_list_item->item);
	}

	free(linked_list_item);
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
    size_t ret_value = ULONG_MAX;
    
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		size_t i = 0;
		linked_list_item_t *curr_item = linked_list_p->first;
		while(curr_item != NULL) {
			if(linked_list_p->compare_fn(curr_item, item) == 0) {
				ret_value = i;
                break;
			}
			i++;
		}

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}
	}
	return ret_value;
}

int linked_list_contains(void *item, linked_list_t *linked_list_p) {
    size_t ret_value = 0;
    
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		linked_list_item_t *curr_item = linked_list_p->first;
		while(curr_item != NULL) {
			if(linked_list_p->compare_fn(curr_item, item) == 0) {
                ret_value = 1;
                break;
			}
		}

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}
	}
	return ret_value;
}

int linked_list_clear(linked_list_t *linked_list_p, void (*data_callback) (void* data)) {
	if(linked_list_p != NULL) {
		linked_list_item_t *curr_item = linked_list_p->first;
		linked_list_item_t *next_item;
		while(curr_item != NULL) {
		  //if(data_callback != NULL) {
				next_item = curr_item->next;
				//data_callback(curr_item);
				linked_list_item_free(curr_item, data_callback);
				curr_item = next_item;
				//	}
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
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		linked_list_item_t* list_item = linked_list_item_new(item_p);
		if(!linked_list_p->first) {
			linked_list_p->first = list_item;
			linked_list_p->last = list_item;
		} else {
			list_item->next = linked_list_p->first;
			linked_list_p->first->prev = list_item;
			linked_list_p->first = list_item;
		}

		linked_list_p->size++;

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}
		return 1;
	}
	return 0;
}

int linked_list_insert_first(void* item_p, linked_list_t *linked_list_p) {
	return linked_list_insert(item_p, linked_list_p);
}

int linked_list_insert_last(void* item_p, linked_list_t *linked_list_p) {
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		linked_list_item_t* list_item = linked_list_item_new(item_p);

		if(!linked_list_p->last) {
			list_item->next = NULL;
			list_item->prev = NULL;
			linked_list_p->first = list_item;
			linked_list_p->last = list_item;
		} else {
			list_item->prev = linked_list_p->last;
			list_item->next = NULL;
			linked_list_p->last->next = list_item;
			linked_list_p->last = list_item;
		}

		linked_list_p->size++;

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}
		return 1;
	}
	return 0;
}

int linked_list_insert_at(size_t index, void* item_p, linked_list_t *linked_list_p) {
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		linked_list_item_t* list_item = linked_list_item_new(item_p);
		linked_list_item_t* list_item_aux = linked_list_p->first;

		for(size_t i = 0; i < index && list_item_aux; i++) {
			list_item_aux = list_item_aux->next;
		}
		
		list_item->prev = list_item_aux->prev;
		list_item->next = list_item_aux;

		if (list_item_aux->prev) {
			list_item_aux->prev->next = list_item;
		}

		list_item_aux->prev = list_item;

		if(linked_list_p->first == list_item_aux) {
			linked_list_p->first = list_item;
		}

		linked_list_p->size++;

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}
		return 1;
	}
	return 0;

}

int linked_list_insert_all(void** item_p, size_t num_items, linked_list_t *linked_list_p) {
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}
		
		int ret_value = 1;

		for (size_t i = 0; i < num_items; i++) {
			if(!linked_list_insert_first(item_p[i], linked_list_p)) {
                ret_value = 0;
                break;
            }
		}

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}
		
		return ret_value;
	}
	
	return 0;
}

int linked_list_insert_all_at(size_t index, void** item_p, size_t num_items, linked_list_t* linked_list_p) {
	return 0;
}


void* linked_list_remove(void *item, linked_list_t *linked_list_p) {
	return NULL;
}

void* linked_list_remove_first(linked_list_t *linked_list_p) {
    void *item = NULL;

	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		if (linked_list_p->first) {
			linked_list_item_t *list_item = linked_list_p->first;
			linked_list_p->first = linked_list_p->first->next;
			/*printf("%i->%i\n", (int)linked_list_p->first->item,
	     (int)linked_list_p->first->next->item);
			 */
			 linked_list_p->first->prev = NULL;
			 linked_list_p->size--;

			 item = list_item->item;
			 linked_list_item_free(list_item, NULL);
		}

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}

	}

	return item;

}

void* linked_list_remove_last(linked_list_t *linked_list_p) {
    void *item = NULL;

	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		if (linked_list_p->last) {
			linked_list_item_t *list_item = linked_list_p->last;
			linked_list_p->last = linked_list_p->last->prev;
			linked_list_p->last->next = NULL;
			linked_list_p->size--;

			item = list_item->item;
			linked_list_item_free(list_item, NULL);
		}

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}
	}

	return item;

}

void* linked_list_remove_at(size_t index, linked_list_t *linked_list_p) {
    void *item = NULL;
    
	if((linked_list_p != NULL) && index >= 0 && index < linked_list_p->size) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		linked_list_item_t* list_item = linked_list_p->first;

		for(size_t i = 0; i < index && list_item; i++) {
			list_item = list_item->next;
		}

		if (list_item) {
			list_item->prev->next = list_item->next;
			list_item->next->prev = list_item->prev;
			linked_list_p->size--;

			item = list_item->item;
			linked_list_item_free(list_item, NULL);
		}

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}

	}
	return item;
}

void** linked_list_remove_range(size_t start, size_t end, linked_list_t* linked_list_p) {
	return NULL;
}



void* linked_list_get(size_t index, linked_list_t *linked_list_p) {
    void *ret_value = NULL;
    
	if(linked_list_p != NULL && index >= 0 && index < linked_list_p->size) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		size_t i = 0;
		linked_list_item_t *curr_item = linked_list_p->first;
		while(curr_item != NULL && i < index) {
			curr_item = curr_item->next;
			i++;
		}

		if(curr_item) { ret_value = curr_item->item; }

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}

	}

	return ret_value;
}

void* linked_list_get_first(linked_list_t *linked_list_p) {
    void *ret_value = NULL;
    
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

        if (linked_list_p->size > 0) {
            ret_value = linked_list_p->first->item;
        }

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}

	}

	return ret_value;
}

void* linked_list_get_last(linked_list_t *linked_list_p) {
    void *ret_value = NULL;
    
	if(linked_list_p != NULL) {
		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&linked_list_p->lock);
		}

		if (linked_list_p->size > 0) {
            ret_value = linked_list_p->last->item;
        }

		if(linked_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&linked_list_p->lock);
		}

	}

	return ret_value;
}

linked_list_t* linked_list_sublist(size_t start, size_t end, linked_list_t *linked_list_p, linked_list_t *sublist) {
	return NULL;
}

void* linked_list_set(size_t index, void* new_item, linked_list_t *linked_list_p) {
	return NULL;
}



void linked_list_print(linked_list_t *linked_list_p, void (*data_callback) (void* data)) {  
	if(linked_list_p != NULL) {
		printf("(%lu)[", linked_list_p->size);
		size_t i = 0;
		linked_list_item_t *curr_item = linked_list_p->first;
		while(curr_item != NULL && i < linked_list_p->size) {
			if (data_callback) {
				data_callback(curr_item->item);
			}else {
				printf("%s,", (char*)curr_item->item);
			}
			curr_item = curr_item->next;
			i++;
		}
		printf("]\n");
	}else {
		printf("The actual linked list is NULL\n");
	}

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


/*********************************************************************/
/***             LinkedList Iterator Implementation               ***/
/*******************************************************************/

linked_list_iterator_t* linked_list_iterator_new(linked_list_t *linked_list_p) {
	assert(linked_list_p);
	linked_list_iterator_t *iterator_p = (linked_list_iterator_t *)malloc(sizeof(linked_list_iterator_t));
	iterator_p->linked_list_p = linked_list_p;
	iterator_p->curr_pos = linked_list_p->first;
	return iterator_p;
}

linked_list_iterator_t* linked_list_iterator_init(linked_list_t *linked_list_p, linked_list_iterator_t *iterator_p) {
	assert(linked_list_p && iterator_p);
	iterator_p->linked_list_p = linked_list_p;
	iterator_p->curr_pos = linked_list_p->first;
	return iterator_p;
}

void linked_list_iterator_free(linked_list_iterator_t *iterator_p) {
	free(iterator_p);
}

void* linked_list_iterator_curr(linked_list_iterator_t *iterator_p) {
	if (iterator_p && iterator_p->curr_pos) { return iterator_p->curr_pos->item; }
	else { return NULL; }
}

linked_list_item_t* linked_list_iterator_list_item_curr(linked_list_iterator_t *iterator_p) {
  if (iterator_p && iterator_p->curr_pos) { return iterator_p->curr_pos; }
  else { return NULL; }
}

void* linked_list_iterator_next(linked_list_iterator_t *iterator_p) {
	void* item = NULL;

	if (iterator_p && iterator_p->curr_pos) {
		item = iterator_p->curr_pos->item;
		//if (iterator_p->curr_pos->next) {
		iterator_p->curr_pos = iterator_p->curr_pos->next; 
		  //}
	}

	return item;
}

void* linked_list_iterator_prev(linked_list_iterator_t *iterator_p) {
  void* item = NULL;

  if (iterator_p && iterator_p->curr_pos) {
      item = iterator_p->curr_pos->item;
  }
  if (iterator_p->curr_pos->prev) { 
      iterator_p->curr_pos = iterator_p->curr_pos->prev; 
  }
  
  return item;

}

void* linked_list_iterator_last(linked_list_iterator_t *iterator_p) {
  void* item = NULL; 

  if (iterator_p && iterator_p->curr_pos) {
    item = iterator_p->curr_pos->item;
  }
  if (iterator_p->linked_list_p) { 
      iterator_p->curr_pos = iterator_p->linked_list_p->last; 
  }
  
  return item;
}


void* linked_list_iterator_first(linked_list_iterator_t *iterator_p) {
  void* item = NULL; 

  if (iterator_p && iterator_p->curr_pos) {
      item = iterator_p->curr_pos->item;
  }
  if (iterator_p->linked_list_p) { 
      iterator_p->curr_pos = iterator_p->linked_list_p->first; 
  }
  
  return item;
}


int linked_list_iterator_insert(void *item, linked_list_iterator_t *iterator_p) {
    linked_list_item_t* list_item = linked_list_item_new(item);
    list_item->next = iterator_p->curr_pos;
    
    if (iterator_p->curr_pos) {
        list_item->prev = iterator_p->curr_pos->prev;
        
        if (iterator_p->curr_pos->prev) {
            iterator_p->curr_pos->prev->next = list_item;
        }
        
        iterator_p->curr_pos->prev = list_item;

        if (iterator_p->curr_pos == iterator_p->linked_list_p->first) {
            iterator_p->linked_list_p->first = list_item;
        }
    } else if (iterator_p->linked_list_p->first == NULL) {
        iterator_p->linked_list_p->first = list_item;
        iterator_p->linked_list_p->last = list_item;
    } else {
        // Insert at the end
        list_item->prev = iterator_p->linked_list_p->last;
        iterator_p->linked_list_p->last->next = list_item;
        iterator_p->linked_list_p->last = list_item;
        iterator_p->curr_pos = list_item;
    }
    
    iterator_p->linked_list_p->size++;
    return 1;
}

void* linked_list_iterator_remove(linked_list_iterator_t *iterator_p) {
    if (iterator_p->curr_pos) {
        linked_list_item_t *list_item = iterator_p->curr_pos;
        void *item = list_item->item;

        if (iterator_p->curr_pos == iterator_p->linked_list_p->first) {
            /*************** ITERATOR IS IN FIRST POSITION *************/
            iterator_p->linked_list_p->first = iterator_p->curr_pos->next;
            if (iterator_p->linked_list_p->first) {
                //THE LIST HAS MORE THAN ONE ELEMENT
                iterator_p->linked_list_p->first->prev = NULL;
            } else {
                //THE LIST HAS ONE ELEMENT
                iterator_p->linked_list_p->last = NULL;
            }
            iterator_p->curr_pos = iterator_p->linked_list_p->first;
        } else if (iterator_p->curr_pos == iterator_p->linked_list_p->last) {
            /*************** ITERATOR IS IN LAST POSITION *************/
            iterator_p->linked_list_p->last = iterator_p->curr_pos->prev;
            iterator_p->linked_list_p->last->next = NULL;
            iterator_p->curr_pos = NULL;
        } else { 
            /*************** ITERATOR IS IN THE MIDDLE *************/
            iterator_p->curr_pos->prev->next = iterator_p->curr_pos->next;
            iterator_p->curr_pos->next->prev = iterator_p->curr_pos->prev;
            iterator_p->curr_pos = iterator_p->curr_pos->next;
        }
        iterator_p->linked_list_p->size--;

        linked_list_item_free(list_item, NULL);

        return item;
    }

    return NULL;
}

/******************************************************************/
