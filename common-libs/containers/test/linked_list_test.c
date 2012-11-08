#include "linked_list_test.h"

void print_item(void *item) {
  printf("%d->", (int)item);
}

int main(int argc, char **argv) {
  printf("********************** LINKED LIST TEST START *****************************\n");
  linked_list_t* list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  printf("Insert From 0 to 5\n");
  for (int i = 5; i >= 0; i--) {
    printf("Insert %d\n", i);
    linked_list_insert((void *)i, list);
  }

  printf("List status is:\n");
  linked_list_print(list, (void *)print_item);
  printf("What is the element that contains %dº position? %d\n", 3, (int)linked_list_get(3, list));
  printf("What is the First element? %d\n", (int)linked_list_get_first(list));
  printf("What is the Last  element? %d\n\n", (int)linked_list_get_last(list));
  
  printf("Now, we remove the first and last items. The first was %d and last was %d\n", (int)linked_list_remove_first(list), (int)linked_list_remove_last(list));
  printf("What is the First and the Last element now? The first is %d and the last is %d\n", (int)linked_list_get_first(list), (int)linked_list_get_last(list));

  printf("**********************  LINKED LIST TEST END  *****************************\n");
}
