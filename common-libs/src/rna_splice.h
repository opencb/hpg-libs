#ifndef RNA_SPLICE_H
#define RNA_SPLICE_H

#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include "cprops/avl.h"
#include "containers/list.h"
#include "buffers.h"
#include "timing.h"

typedef struct allocate_buffers {
  write_batch_t *write_exact_sp;
  write_batch_t *write_extend_sp;
  FILE *fd_extend;
  FILE *fd_exact;
}allocate_buffers_t;


/**
 * @brief Structure for store splice junctions end information.
 * 
 * Structure for store the strand, end and extend end genome coordinate.  
 */
typedef struct splice_end{
  size_t end;  /**< exact end genome coordinate */
  unsigned char strand;  /**< genome strand */
  size_t reads_number;  /**< number of reads with the same start and end coordiante */
  size_t splice_end_extend;  /**< extend end genome coordinate */
}splice_end_t;

/**
 * @brief  Reserve memory for a new @a splice_end_t and initialize it.
 * @param  strand genome strand
 * @param  end exact end genome coordinate
 * @param  splice_end extend end genome coordinate 
 * @param  splice_end_p @a splice_end_t structure 
 * @return value of new node
 * 
 * Reserve memory for a new @a splice_end_t structure and initialize their fields. 
 */
splice_end_t* new_splice_end(unsigned char strand, size_t end, 
			     size_t splice_end);


void free_splice_end(splice_end_t *splice_end_p);

/**
 * @brief Structure for one node of avl tree.
 *
 * Structure for store all splice junctions data information. It represents 
 * one node of the avl tree.   
 */
typedef struct node_element_splice{
  size_t splice_start;  /**< exact start genome coordinate */
  size_t splice_start_extend;  /**< splice_start_extend extend end genome coordinate */
  size_t maximum_allocate_ends; /**< maximum_allocate_ends maximum number for allocate splice junction */
  size_t number_allocate_ends;  /**< number of ends are allocate in the array */
  splice_end_t **allocate_ends;  /**< array for allocate ends */
}node_element_splice_t;

/**
 * @brief Structure for store splice junction. 
 * 
 * Structure that contains one avl tree for store the splice junctions found in
 * a particular chromosome and one mutex for insertions synchronized.
 */
typedef struct allocate_splice_elements{
  cp_avltree *avl_splice; /**< avl tree for store splice junctions */
  pthread_mutex_t mutex;  /**< mutex for insertions synchronized */
}allocate_splice_elements_t;

/**
 * @brief  Compare one node of avl with one unsigned int.
 * @param  a node of avl tree
 * @param  b unsigned int number to compare with the node 
 * @return 0 if the two values are the same. 1 if a start splice are major and -1 in otherwise
 *
 * Compare the splice start field of @a node_element_splice_t with the second parameter. 
 */
int node_compare(node_element_splice_t* a, size_t b);

/**
 * @brief  Reserve memory for new node and initialize it.
 * @param  b splice start junction value
 * @return The copy of new element
 * 
 * Compare the splice start field of @a node_element_splice_t with the second parameter. 
 */
node_element_splice_t* node_copy(size_t b);

/**
 * @brief  Insert a new end splice junction.
 * @param  splice_end_p new splice junction end for store in the ends array
 * @return Structure with splice junctions ends result 
 * 
 * Insert the new splice junction end into the array and then checks if it is full. 
 * If it is full reserve more memory for it. 
 */
node_element_splice_t* insert_end_splice(splice_end_t *splice_end_p, 
					 node_element_splice_t *element_p);


/**
 * @brief  First insert a new start splice junction and finally insert a new end splice junction.
 * @param  chromosome splice junction chromosome location
 * @param  strand splice junction strand location
 * @param  end exact end genome coordinate 
 * @param  splice_start extend start genome coordinate
 * @param  splice_end extend end genome coordinate
 * @param  element_p structure for store splice junctions ends
 * @return Structure with splice junctions ends result
 * 
 * First insert the extend splice start in the correct position and then increments the number
 * of splice junctions that cover this coordinates and insert the splice junction end coordinate. 
 */
node_element_splice_t* search_and_insert_end_splice(unsigned int chromosome, unsigned char strand, 
						    size_t end , size_t splice_start, 
						    size_t splice_end, 
						    node_element_splice_t *element_p);

/**
 * @brief  Initializer of structures for store splice junctions.
 * @param  chromosomes_splice_avls_p structure for store avls trees 
 * @return Structure for store avls trees result 
 * 
 * Firt initialize and reserve memory for each avl need and then initialize each mutex. 
 */
allocate_splice_elements_t* init_allocate_splice_elements(allocate_splice_elements_t* chromosomes_splice_avls_p);


/**
 * @brief  Insert new splice junction in the corresponding avl tree.
 * @param  chromosome chromosome genome location
 * @param  strand strand genome location
 * @param  end exact genome end location
 * @param  start exact genome start location
 * @param  splice_start extend start genome location
 * @param  splice_end extend end genome location
 * @param  chromosomes_splice_avls_p structure for store one avl tree for each chromosome
 * @return Structure for store avls trees result
 * 
 * First select the corresponding avl tree for the splice junction chromosome. Then search in the 
 * splice junction start is already hosted. Finally insert the splice junction end calling to @a 
 * search_and_insert_end_splice. 
 */
allocate_splice_elements_t* allocate_new_splice(unsigned int chromosome, unsigned char strand, 
						size_t end, size_t start, 
						size_t splice_start, size_t splice_end, 
						allocate_splice_elements_t* chromosomes_splice_avls_p);

/**
 * @brief  Recursive function for process all nodes of one avl tree.
 * @param  node one node of the avl tree 
 * @param  chromosome chromosome genome location
 * @param  write_list_p list for store write batches
 * @param  write_size size of one write batch 
 * 
 * Recursive function that process first the left node of the avl tree, then de root node
 * and finally the right node. 
 */
allocate_buffers_t* process_avlnode_in_order(cp_avlnode *node, unsigned int chromosome, list_t* write_list_p, unsigned int write_size,  allocate_buffers_t* buffers_write_p);

/**
 * @brief  Sequential function for process all splice junctions ends allocate in one node.
 * @param  node one node of the avl tree
 * @param  chromosome chromosome genome location
 * @param  write_list_p list for store write batches
 * @param  write_size size of one write batch 
 * 
 * Sequential function that process all splice junctions end stored in one node. For each end
 * calls the function @a pack_junction for generate the buffer with all information corresponding
 * with the splice junction found, then with all it is created a write batch. 
 */
allocate_buffers_t* process_avlnode_ends_in_order(node_element_splice_t *node, unsigned int chromosome, list_t* write_list_p, unsigned int write_size, allocate_buffers_t *buffers_write_p);


/**
 * @brief  Sequential function for process each chromosome avl.
 * @param  chromosome_avls structure for store avls trees
 * @param  write_list_p list for store write batches
 * @param  write_size size of one write batch 
 * 
 * Sequential function for process each chromosome avl tree and free.
 */
void write_chromosome_avls(allocate_splice_elements_t *chromosome_avls, 
			   list_t* write_list_p, char *extend_sp, char *exact_sp, 
			   unsigned int write_size);

//====================================================================================================

void cp_avlnode_print_new(cp_avlnode *node, int level);

void cp_avlnode_print_in_order(cp_avlnode *node);

void node_list_print(node_element_splice_t *node);

#endif
