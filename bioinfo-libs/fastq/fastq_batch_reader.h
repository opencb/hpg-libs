
#ifndef FASTQ_BATCH_READER_H
#define FASTQ_BATCH_READER_H

#include <stdio.h>
#include <stdio.h>
#include <pthread.h>

#include "fastq_file.h"
#include "fastq_read.h"
#include "fastq_batch_list.h"
#include "list.h"

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Fastq batch reader 
* 
* Structure for implementing a fastq batch reader server
*/
typedef struct fastq_batch_reader {
  int alive;				/**< Alive indicator of the server. */
  pthread_mutex_t alive_lock;		/**< Alive variable lock. */

  int eof;				/**< End Of File indicator. */
  pthread_mutex_t eof_lock;		/**< End Of File indicator lock. */

  size_t batch_size;			/**< Fastq batch size (in bytes). */
  int batch_list_max_length;		/**< Fastq batch list maximum length. */

  pthread_t thread;			/**< pthread element. */

  int source_id;			/**< Source id of the fastq file (pair1 or pair2). */

  fastq_file_t* fastq_file_p;		/**< Pointer to the fastq file object. */
  fastq_batch_list_t* batch_list_p;	/**< Pointer to the fastq batch list. */
  list_t* qc_batch_list_p;		/**< Pointer to the qc batch list. */
} fastq_batch_reader_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Creates a new fastq batch reader
*  @param filename name of the fastq file
*  @param source_id source id of the fastq file (pair1 or pair2)
*  @param batch_list_p pointer to the fastq batch list
*  @param batch_size size of the fastq batch
*  @param qc_batch_list_p pointer to the qc batch list
*  @param batch_list_max_length Fastq batch list maximum length
*  @return fastq_batch_reader_t pointer to the reader
*  
*  Free a fastq batch list item (all = 0 -> fastq batch not freed, all = 1 -> fastq batch freed)
*/
fastq_batch_reader_t* fastq_batch_reader_new(char* filename, int source_id, fastq_batch_list_t* batch_list_p, size_t batch_size, list_t* qc_batch_list_p, int batch_list_max_length);

/**
*  @brief Frees a fastq batch reader
*  @param[in,out] reader_p pointer to the fastq batch reader
*  @return void
*  
*  Frees a fastq batch reader 
*/
void fastq_batch_reader_free(fastq_batch_reader_t* reader_p);

/**
*  @brief Starts a fastq batch reader
*  @param[in,out] reader_p pointer to the fastq batch reader
*  @return void
*  
*  Starts a fastq batch reader 
*/
void fastq_batch_reader_start(fastq_batch_reader_t* reader_p);

/**
*  @brief Joins/waits a fastq batch reader
*  @param reader_p pointer to the fastq batch reader
*  @return void
*  
*  Joins/waits until the pthread of the fastq batch reader is finished
*/
unsigned int fastq_batch_reader_join(fastq_batch_reader_t* reader_p);

/**
*  @brief Joins/waits a fastq batch reader
*  @param reader_p pointer to the fastq batch reader
*  @return void
*  
*  Joins/waits until the pthread of the fastq batch reader is finished
*/
fastq_batch_list_item_t* fastq_batch_reader_next_batch(fastq_batch_reader_t* reader_p);

/**
*  @brief Getter of the eof variable of the reader
*  @param reader_p pointer to the fastq batch reader
*  @return void
*  
*  Getter method of the eof variable of the reader
*/
int fastq_batch_reader_get_eof(fastq_batch_reader_t* reader_p);

/**
*  @brief Setter of the eof variable of the reader
*  @param reader_p pointer to the fastq batch reader
*  @param eof value of the eof variable to set
*  @return void
*  
*  Setter method of the eof variable of the reader
*/
void fastq_batch_reader_set_eof(fastq_batch_reader_t* reader_p, int eof);

/**
*  @brief Getter of the alive variable of the reader
*  @param reader_p pointer to the fastq batch reader
*  @return void
*  
*  Getter method of the alive variable of the reader
*/
int fastq_batch_reader_get_alive(fastq_batch_reader_t* reader_p);

/**
*  @brief Setter of the alive variable of the reader
*  @param reader_p pointer to the fastq batch reader
*  @param alive value of the alive variable to set
*  @return void
*  
*  Setter method of the alive variable of the reader
*/
void fastq_batch_reader_set_alive(fastq_batch_reader_t* reader_p, int alive);

/**
*  @brief pthread implementation of the fastq batch reader
*  @param param_p 
*  @return void
*  
*  Function with pthread implementation of the fastq batch reader
*/
void* fastq_batch_reader_thread_function(void* param_p);


#endif	/*  FASTQ_BATCH_READER_H  */
