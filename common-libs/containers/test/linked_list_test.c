#include "linked_list_test.h"

void print_item(void *item) {
  printf("%d->", (int)item);
}

int main(int argc, char **argv) {
  printf("********************** LINKED LIST TEST START *****************************\n");

  printf("----------------------INSERT AND REMOVE TEST --------------------------\n");

  linked_list_t* list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  printf("Insert From 0 to 5:\n");
  for (int i = 5; i >= 0; i--) {
    linked_list_insert((void *)i, list);
  }

  printf("List status is:\n");
  linked_list_print(list, (void *)print_item);
  
  printf("\n");
  printf("What is the element that contains %dº position? %d\n", 3, (int)linked_list_get(3, list));
  printf("What is the First element? %d\n", (int)linked_list_get_first(list));
  printf("What is the Last  element? %d\n", (int)linked_list_get_last(list));
    
  void* item = linked_list_remove_first(list);
  printf("Remove First Item. It is %d\n", item);
  
  item = linked_list_remove_last(list);
  printf("Remove Last Item. It is %d\n", item);

  item = linked_list_remove_at(2, list);
  printf("Remove item at 2 position. It is %d\n", item);
  
  printf("List status is:\n");
  linked_list_print(list, (void *)print_item);  
  
  printf("And Now insert %i at 2 position.\n", item);
  linked_list_insert_at(2, item, list);

  printf("List status is:\n");
  linked_list_print(list, (void *)print_item);  

  printf("And Now insert 0 to First position and 5 to End position.\n");
  linked_list_insert_first((void *)0, list);
  linked_list_insert_last((void *)5, list);
  printf("List status is:\n");
  linked_list_print(list, (void *)print_item);
  
  printf("----------------------INSERT AND REMOVE TEST END --------------------------\n");

  printf("-------------------------- ITERATORS TEST ---------------------------\n");

  printf("Linked list initial status:\n");
  linked_list_print(list, (void *)print_item);

  linked_list_iterator_t* iterator = linked_list_iterator_new(list);
  printf("The iterator are in %d. And move it to next position\n", 
	 linked_list_iterator_curr(iterator));
  linked_list_iterator_next(iterator);
  
  printf("The iterator are now in %d. And move it to prev position\n", 
	 linked_list_iterator_curr(iterator));
  linked_list_iterator_prev(iterator);

  //linked_list_print(list, (void *)print_item);
  printf("The iterator are in %d. If the iterator is moved out of first element or last element, it is not moved.\n", 
	 linked_list_iterator_curr(iterator));

  linked_list_iterator_prev(iterator);
  printf("The iterator are in %d\n", 
	 linked_list_iterator_curr(iterator));
  
  linked_list_iterator_insert(8, iterator);
  printf("Linked list actual status:\n");
  linked_list_print(list, (void *)print_item);

  printf("-------------------------- ITERATORS TEST END ---------------------------\n");

  printf("**********************  LINKED LIST TEST END  *****************************\n");
}
