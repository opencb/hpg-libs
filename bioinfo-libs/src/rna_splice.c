#include "rna_splice.h"

unsigned int junction_id = 0;
size_t total_splice = 0;


int node_compare(node_element_splice_t* a, unsigned int b) {
 if(a->splice_start == b){ 
   return 0;
 }else if(a->splice_start < b){
   return -1;
 }else{
   return 1;
 }
}

node_element_splice_t* node_copy(unsigned int b) {
 node_element_splice_t *node = (node_element_splice_t *)malloc(sizeof(node_element_splice_t));

 //assert(node != NULL);
 
 node->maximum_allocate_ends = 10;
 node->number_allocate_ends = 0;
 node->allocate_ends = (splice_end_t **)malloc(node->maximum_allocate_ends * sizeof(splice_end_t *));
 if (node->allocate_ends == NULL){ exit(-1); }
 
 node->splice_start = b;

 return node;
}

void node_free(node_element_splice_t* a) {  
  for (unsigned int i = 0; i < a->number_allocate_ends; i++) {
    free_splice_end(a->allocate_ends[i]);
  }

  free(a->allocate_ends);
  free(a);
}


node_element_splice_t* insert_end_splice(splice_end_t *splice_end_p, node_element_splice_t *element_p){
  element_p->allocate_ends[element_p->number_allocate_ends] = splice_end_p;
  element_p->number_allocate_ends++;

  if(element_p->number_allocate_ends >= element_p->maximum_allocate_ends){
    element_p->maximum_allocate_ends = element_p->maximum_allocate_ends * 2; 
    element_p->allocate_ends = (splice_end_t **) realloc (element_p->allocate_ends, element_p->maximum_allocate_ends * sizeof (splice_end_t *));
    if(element_p->allocate_ends == NULL){exit(-1);}
  }

  return element_p;
}

node_element_splice_t* search_and_insert_end_splice(unsigned int chromosome, unsigned char strand, 
						    unsigned int end, unsigned int splice_start, 
						    unsigned int splice_end, 
						    node_element_splice_t *element_p){
  unsigned int i;

  if(element_p->splice_start_extend > splice_start){
    element_p->splice_start_extend = splice_start;
  }	  
  
  for(i = 0; i < element_p->number_allocate_ends; i++){
    if( (element_p->allocate_ends[i]->end == end) && (element_p->allocate_ends[i]->strand == strand) ) { 

      element_p->allocate_ends[i]->reads_number++;
      
      if(element_p->allocate_ends[i]->splice_end_extend < splice_end){
	element_p->allocate_ends[i]->splice_end_extend = splice_end;
      }  
	  
      return element_p;
    }
  }
  
  splice_end_t *splice_end_p = new_splice_end(strand, end, splice_end, splice_end_p);
  
  return insert_end_splice(splice_end_p, element_p);
}


allocate_splice_elements_t* init_allocate_splice_elements(allocate_splice_elements_t* chromosomes_avls_p){
   int i;
   for(i = 0; i < CHROMOSOME_NUMBER; i++){
     chromosomes_avls_p[i].avl_splice = cp_avltree_create_by_option(COLLECTION_MODE_NOSYNC | 
								    COLLECTION_MODE_COPY   |
								    COLLECTION_MODE_DEEP, 
								    (cp_compare_fn) node_compare, 
								    (cp_copy_fn) node_copy, 
								    (cp_destructor_fn)node_free, 
								    (cp_copy_fn) node_copy, 
								    (cp_destructor_fn)node_free);
     
     if(chromosomes_avls_p[i].avl_splice == NULL){exit(-1);}
     pthread_mutex_init(&(chromosomes_avls_p[i].mutex), NULL);
   }
   
   return chromosomes_avls_p;
}

allocate_splice_elements_t* allocate_new_splice(unsigned int chromosome, unsigned char strand, 
						unsigned int end, unsigned int start, 
						unsigned int splice_start, unsigned int splice_end, 
						allocate_splice_elements_t* chromosome_avls_p){
  node_element_splice_t *node;

  node = (node_element_splice_t *)cp_avltree_get(chromosome_avls_p[chromosome].avl_splice, start);

  if(node == NULL) {
    node = cp_avltree_insert(chromosome_avls_p[chromosome].avl_splice, start, start);
    node->splice_start_extend = splice_start;
  }
  
  node = search_and_insert_end_splice(chromosome, strand, end, splice_start, splice_end, node);

  return chromosome_avls_p;
}

splice_end_t* new_splice_end(unsigned char strand, unsigned int end, 
			     unsigned int splice_end, splice_end_t *splice_end_p){

  splice_end_p = (splice_end_t *)malloc(sizeof(splice_end_t));

  if(splice_end_p == NULL){exit(-1);}

  splice_end_p->strand = strand;
  splice_end_p->end = end;
  splice_end_p->splice_end_extend = splice_end;
  splice_end_p->reads_number = 1;
  
  return splice_end_p;
}


void free_splice_end(splice_end_t *splice_end_p){
  free(splice_end_p);
}


allocate_buffers_t * process_avlnode_in_order(cp_avlnode *node, unsigned int chromosome, 
					      list_t* write_list_p, unsigned int write_size,   allocate_buffers_t *allocate_batches){
  
  
  if (node->left) { 
    allocate_batches = process_avlnode_in_order(node->left, chromosome, write_list_p, write_size, allocate_batches);
  }
  
  allocate_batches = process_avlnode_ends_in_order((node_element_splice_t *)node->value, chromosome, write_list_p, write_size,
						   allocate_batches);
  
  if (node->right) { 
    allocate_batches = process_avlnode_in_order(node->right, chromosome,  write_list_p, write_size, allocate_batches);
  }
  
  return allocate_batches;
  // return exact_splice_write_p;

}

allocate_buffers_t* process_avlnode_ends_in_order(node_element_splice_t *node, unsigned int chromosome,
					     list_t* write_list_p, unsigned int write_size, allocate_buffers_t *allocate_batches) {
  int i;
  char strand[2] = {'+', '-'};
  list_item_t* item_p = NULL;
  unsigned int bytes_exact, bytes_extend;
  allocate_batches->write_exact_sp;
  //  write_batch_t* extend_splice_write_p = write_batch_new(write_size, SPLICE_EXTEND_FLAG);

  for(i = 0; i < node->number_allocate_ends; i++){
    if(( allocate_batches->write_exact_sp->size + 100) > write_size) {
      item_p = list_item_new(0, WRITE_ITEM,  allocate_batches->write_exact_sp);
      list_insert_item(item_p, write_list_p);
      allocate_batches->write_exact_sp = write_batch_new(write_size, SPLICE_EXACT_FLAG);
    } 
    if(( allocate_batches->write_extend_sp->size + 100) > write_size) {
      item_p = list_item_new(0, WRITE_ITEM,  allocate_batches->write_extend_sp);
      list_insert_item(item_p, write_list_p);
      allocate_batches->write_extend_sp = write_batch_new(write_size, SPLICE_EXTEND_FLAG);
    } 

    bytes_exact = pack_junction(chromosome, node->allocate_ends[i]->strand, node->splice_start, node->allocate_ends[i]->end, junction_id, node->allocate_ends[i]->reads_number, &(allocate_batches->write_exact_sp->buffer_p[allocate_batches->write_exact_sp->size])); 
    
    bytes_extend = pack_junction(chromosome, node->allocate_ends[i]->strand, node->splice_start_extend, node->allocate_ends[i]->splice_end_extend, junction_id, node->allocate_ends[i]->reads_number, &(allocate_batches->write_extend_sp->buffer_p[allocate_batches->write_extend_sp->size])); 
    
    allocate_batches->write_exact_sp->size += bytes_exact;
    allocate_batches->write_extend_sp->size += bytes_extend;
    
    total_splice += node->allocate_ends[i]->reads_number;
    junction_id++;
  }
  return allocate_batches;
  //return exact_splice_write_p;
}

/*
void process_avlnode_in_order(cp_avlnode *node, unsigned int chromosome, 
			      list_t* write_list_p, unsigned int write_size){
  if (node->left) { 
    process_avlnode_in_order(node->left, chromosome, write_list_p, write_size);
  }

  process_avlnode_ends_in_order((node_element_splice_t *)node->value, chromosome, write_list_p, write_size);
  
  if (node->right) { 
    process_avlnode_in_order(node->right, chromosome,  write_list_p, write_size);
  }
}

void process_avlnode_ends_in_order(node_element_splice_t *node, unsigned int chromosome,
				   list_t* write_list_p, unsigned int write_size) {
  int i;
  char strand[2] = {'+', '-'};
  list_item_t* item_p = NULL;
  unsigned int bytes_exact, bytes_extend;
  write_batch_t* exact_splice_write_p  = write_batch_new(write_size, SPLICE_EXACT_FLAG);
  write_batch_t* extend_splice_write_p = write_batch_new(write_size, SPLICE_EXTEND_FLAG);

  
  //Order Ends
  /* splice_end_t *aux;
  for(unsigned int j = 0; j < node->number_allocate_ends; j++){
      for(unsigned int z = j; z < node->number_allocate_ends; z++){
	if (node->allocate_ends[j]->strand > node->allocate_ends[z]->strand) {
	  aux = node->allocate_ends[j];
	  node->allocate_ends[j] =  node->allocate_ends[z];
	  node->allocate_ends[z] = aux;
	}else if (node->allocate_ends[j]->end > node->allocate_ends[z]->end) {
	    aux = node->allocate_ends[j];
	    node->allocate_ends[j] =  node->allocate_ends[z];
	    node->allocate_ends[z] = aux;
	}
      }
      }
  

  for(i = 0; i < node->number_allocate_ends; i++){
    // More than one read
    if((exact_splice_write_p->size + 100) > write_size) {
      item_p = list_item_new(0, WRITE_ITEM, exact_splice_write_p);
      list_insert_item(item_p, write_list_p);
      exact_splice_write_p = write_batch_new(write_size, SPLICE_EXACT_FLAG);
    } 
    
    if((extend_splice_write_p->size + 100) > write_size) {
      item_p = list_item_new(0, WRITE_ITEM, extend_splice_write_p);
      list_insert_item(item_p, write_list_p);
      exact_splice_write_p = write_batch_new(write_size, SPLICE_EXTEND_FLAG);
    } 
    

    bytes_exact = pack_junction(chromosome, node->allocate_ends[i]->strand, node->splice_start, node->allocate_ends[i]->end, junction_id, node->allocate_ends[i]->reads_number, &(exact_splice_write_p->buffer_p[exact_splice_write_p->size])); 
    
    bytes_extend = pack_junction(chromosome, node->allocate_ends[i]->strand, node->splice_start_extend, node->allocate_ends[i]->splice_end_extend, junction_id, node->allocate_ends[i]->reads_number, &(extend_splice_write_p->buffer_p[extend_splice_write_p->size])); 
    
    exact_splice_write_p->size += bytes_exact;
    extend_splice_write_p->size += bytes_extend;
    
    total_splice += node->allocate_ends[i]->reads_number;
    junction_id++;
  }
  
  if(exact_splice_write_p != NULL) {
    list_item_t* item_p = NULL;
    if(exact_splice_write_p->size > 0) {
      item_p = list_item_new(0, WRITE_ITEM, exact_splice_write_p);
      list_insert_item(item_p, write_list_p);
      printf("Insert 1\n");
    } else {
      write_batch_free(exact_splice_write_p);
    }
  }
  
  if(extend_splice_write_p != NULL) {
    list_item_t* item_p = NULL;
    if(extend_splice_write_p->size > 0) {
      item_p = list_item_new(0, WRITE_ITEM, extend_splice_write_p);
      list_insert_item(item_p, write_list_p);
      printf("Insert 2\n");
    } else {
      write_batch_free(extend_splice_write_p);
    }
  }
  
}
*/


void process_and_free_chromosome_avls(allocate_splice_elements_t *chromosome_avls, 
				      list_t* write_list_p, unsigned int write_size) {
  int c;
  allocate_buffers_t *allocate_batches = (allocate_buffers_t *)malloc(sizeof(allocate_buffers_t));
  write_batch_t *exact_splice_write_p;
  write_batch_t *extend_splice_write_p;

  for(c = 0; c < CHROMOSOME_NUMBER; c++){
    if(chromosome_avls[c].avl_splice->root != NULL) {
      allocate_batches->write_exact_sp  = write_batch_new(write_size, SPLICE_EXACT_FLAG);
      allocate_batches->write_extend_sp  = write_batch_new(write_size, SPLICE_EXTEND_FLAG);
      //allocate_batches->write_extend_sp  = write_batch_new(1000, SPLICE_EXTEND_FLAG);

      allocate_batches = process_avlnode_in_order(chromosome_avls[c].avl_splice->root, c, write_list_p, write_size, allocate_batches);
      
      exact_splice_write_p = allocate_batches->write_exact_sp;
      extend_splice_write_p = allocate_batches->write_extend_sp;
      
      if(exact_splice_write_p != NULL) {
	list_item_t* item_p = NULL;
	if(exact_splice_write_p->size > 0) {
	  item_p = list_item_new(0, WRITE_ITEM, exact_splice_write_p);
	  list_insert_item(item_p, write_list_p);
	} else {
	  write_batch_free(exact_splice_write_p);
	}
      }

      if(extend_splice_write_p != NULL) {
	list_item_t* item_p = NULL;
	if(extend_splice_write_p->size > 0) {
	  item_p = list_item_new(0, WRITE_ITEM, extend_splice_write_p);
	  list_insert_item(item_p, write_list_p);
	  } else {
	  write_batch_free(extend_splice_write_p);
	}
      }
      
    }//end IF chromosome splice not NULL
    cp_avltree_destroy(chromosome_avls[c].avl_splice);
  }
  
  free(allocate_batches);

  if (statistics_on) { 
    statistics_set(TOTAL_ST, 3, total_splice, statistics_p);
  }
  
  list_decr_writers(write_list_p);
}

//====================================================================================================

void cp_avlnode_print_new(cp_avlnode *node, int level){
  int i;
  if (node->right) cp_avlnode_print_new(node->right, level + 1);
  for (i = 0; i < level; i++) printf("  . ");
  printf("(%d) [%i => %i]", node->balance, ((node_element_splice_t *)node->key)->splice_start, ((node_element_splice_t *)node->value)->splice_start);
  node_list_print((node_element_splice_t *)node->value);
  if (node->left) cp_avlnode_print_new(node->left, level + 1);
}

void cp_avlnode_print_in_order(cp_avlnode *node){
    
  if (node->left) cp_avlnode_print_in_order(node->left);
  printf("(%d) [%i => %i]", node->balance, ((node_element_splice_t *)node->key)->splice_start, ((node_element_splice_t *)node->value)->splice_start);
  node_list_print((node_element_splice_t *)node->value);
  if (node->right) cp_avlnode_print_in_order(node->right);
}

void node_list_print(node_element_splice_t *node){
  int i;
  printf("::(%i)Ends{", node->number_allocate_ends);
  for(i = 0; i < node->number_allocate_ends; i++)
    printf("|%i-%i|#", node->allocate_ends[i]->end, node->allocate_ends[i]->reads_number);
    printf("}\n");
}
