#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

//#include "list.h"
//#include "timing.h"

//#include "fastq_file.h"
//#include "fastq_ex_batch.h"
#include "fastq_batch_reader.h"

/*-----------------------------------------------------
 * init structure for the batch reader
 *----------------------------------------------------*/

void fastq_batch_reader_input_init(char* filename, int batch_size, list_t* list_p, fastq_batch_reader_input_t* input_p) {
  input_p->filename = filename;
  input_p->batch_size = batch_size;
  input_p->list_p = list_p;
}

/*-----------------------------------------------------
 * this functions reads read from disk, and save them
 * into a list to be processed by the next thread
 *----------------------------------------------------*/

void fastq_batch_reader(fastq_batch_reader_input_t* input_p) {
  
  unsigned int total_reads = 0;

  /*struct timespec ts;
  ts.tv_sec = 1;
  ts.tv_nsec = 0;*/

  char* filename = input_p->filename;
  int batch_size = input_p->batch_size;
  list_t* list_p = input_p->list_p;

  printf("fastq_ex_batch_reader (%i): START, for file %s\n", omp_get_thread_num(), filename);
		
  int num_reads = 0, num_batches = 0;
  fastq_batch_t* batch_p;
  list_item_t *item_p = NULL;

  fastq_file_t* file_p = fastq_fopen(filename);
  
  while (1) {
    //if (time_on) { timing_start(FASTQ_READER_INDEX, 0, timing_p); }
    // allocationg memory for the current fastq batch
    //
    batch_p = fastq_batch_new(batch_size);
    //batch_p->source_id = reader_p->source_id;
      
    // read reads from file
    //
    num_reads = fastq_fread_batch_max_size(batch_p, batch_size, file_p);
    total_reads	+= num_reads;

    //if (time_on) { stop_timer(t1_read, t2_read, read_time); }
    
    // if no reads, free memory and go out....
    //
    if (num_reads==0) {
      fastq_batch_free(batch_p);
      break;
    }

    // otherwise, create a new batch object..
    // and insert this batch to the corresponding list
    //
    item_p = list_item_new(num_batches, READ_ITEM, batch_p);
    //if (time_on) { timing_stop(FASTQ_READER_INDEX, 0, timing_p); }

    list_insert_item(item_p, list_p);
  
    //printf("fastq_ex_batch_reader: reading batch %i done !!\n", num_batches);
    num_batches++;
    
    //if (num_batches==2) break;
  } // end of batch loop
  
  list_decr_writers(list_p);
  fastq_fclose(file_p);

  printf("fastq_ex_batch_reader: END, %i total reads (%i batches), for file %s\n", total_reads, num_batches, filename);
}

/*-----------------------------------------------------
 *----------------------------------------------------*/


















/* ******************************************************
 *    		Function implementations  		*
 * ******************************************************

fastq_batch_reader_t* fastq_batch_reader_new(char* filename, int source_id, fastq_batch_list_t* fastq_batch_list_p, size_t batch_size, list_t* qc_batch_list_p, int batch_list_max_length) {
    fastq_batch_reader_t* fastq_batch_reader_p = (fastq_batch_reader_t*) calloc(1, sizeof(fastq_batch_reader_t));

    fastq_batch_reader_p->fastq_file_p = fastq_fopen_mode(filename, "r");    // open the input file
    fastq_batch_reader_p->source_id = source_id;

    fastq_batch_reader_p->batch_size = batch_size;
    fastq_batch_reader_p->batch_list_max_length = batch_list_max_length;

    fastq_batch_reader_p->eof = 0;
    pthread_mutex_init(&(fastq_batch_reader_p->eof_lock), NULL);

    fastq_batch_reader_p->alive = 1;
    pthread_mutex_init(&(fastq_batch_reader_p->alive_lock), NULL);

    fastq_batch_reader_p->batch_list_p = fastq_batch_list_p;
    fastq_batch_list_incr_producers(fastq_batch_list_p);

    fastq_batch_reader_p->qc_batch_list_p = qc_batch_list_p;

    return fastq_batch_reader_p;
}

void fastq_batch_reader_free(fastq_batch_reader_t* fastq_batch_reader_p) {
    // close the input file and exit
    fastq_fclose(fastq_batch_reader_p->fastq_file_p);
    free(fastq_batch_reader_p);
}

void fastq_batch_reader_start(fastq_batch_reader_t* fastq_batch_reader_p) {
    // create and launch pthread
    pthread_create(&(fastq_batch_reader_p->thread), NULL, fastq_batch_reader_thread_function, (void*) fastq_batch_reader_p);
}

unsigned int fastq_batch_reader_join(fastq_batch_reader_t* fastq_batch_reader_p) {
    char log_message[50];
    sprintf(log_message, "THREAD-READ: waiting for batch server thread....\n");
    LOG_DEBUG(log_message);

    void* r;
    pthread_join(fastq_batch_reader_p->thread, &r);
    return (uintptr_t) r;
}

fastq_batch_list_item_t* fastq_batch_reader_next_batch(fastq_batch_reader_t* fastq_batch_reader_p) {
    return fastq_batch_list_remove(fastq_batch_reader_p->batch_list_p);
}

void fastq_batch_reader_set_eof(fastq_batch_reader_t* fastq_batch_reader_p, int eof) {
    pthread_mutex_lock(&(fastq_batch_reader_p->eof_lock));
    fastq_batch_reader_p->eof = eof;
    pthread_mutex_unlock(&(fastq_batch_reader_p->eof_lock));
}

int fastq_batch_reader_get_eof(fastq_batch_reader_t* fastq_batch_reader_p) {
    int eof;
    pthread_mutex_lock(&(fastq_batch_reader_p->eof_lock));
    eof = fastq_batch_reader_p->eof;
    pthread_mutex_unlock(&(fastq_batch_reader_p->eof_lock));
    return eof;
}

void fastq_batch_reader_set_alive(fastq_batch_reader_t* fastq_batch_reader_p, int alive) {
    pthread_mutex_lock(&(fastq_batch_reader_p->alive_lock));
    fastq_batch_reader_p->alive = alive;
    pthread_mutex_unlock(&(fastq_batch_reader_p->alive_lock));
}

int fastq_batch_reader_get_alive(fastq_batch_reader_t* fastq_batch_reader_p) {
    int alive;
    pthread_mutex_lock(&(fastq_batch_reader_p->alive_lock));
    alive = fastq_batch_reader_p->alive;
    pthread_mutex_unlock(&(fastq_batch_reader_p->alive_lock));
    return alive;
}

void* fastq_batch_reader_thread_function(void* param_p) {
    unsigned int total_reads = 0;
    char log_message[100];

    unsigned int usecs = 1000; // micro seconds to wait (1 ms.)

    fastq_batch_reader_t* fastq_batch_reader_p = (fastq_batch_reader_t*) param_p;    // cast void* param_p to fastq_batch_reader_t

    sprintf(log_message, "Thread-READ: START, for file %s\n", fastq_batch_reader_p->fastq_file_p->filename);
    LOG_DEBUG(log_message);

    fastq_batch_reader_set_eof(fastq_batch_reader_p, 0);
    fastq_batch_reader_set_alive(fastq_batch_reader_p, 1);

    int num_reads = 0, num_batchs = 0;
    fastq_batch_t* fastq_batch_p;
    fastq_batch_list_item_t *item_p = NULL;

    while (1) {        
        LOG_DEBUG("Thread-READ: Before allocating memory for fastq_batch_p...");
        fastq_batch_p = fastq_batch_new();    // allocationg memory for the current fastq batch

        fastq_batch_init(fastq_batch_p, fastq_batch_reader_p->batch_size);
        fastq_batch_p->source_id = fastq_batch_reader_p->source_id;

        while ((fastq_batch_list_length(fastq_batch_reader_p->batch_list_p) >= fastq_batch_reader_p->batch_list_max_length) ||
               (fastq_batch_reader_p->batch_list_p->length_by_source_id[fastq_batch_reader_p->source_id] > (fastq_batch_reader_p->batch_list_max_length / 2)) ||
               (fastq_batch_reader_p->qc_batch_list_p->length > (fastq_batch_reader_p->batch_list_max_length))) {

            // Delay for a bit
            sched_yield();
            usleep(usecs);
        }

        if (time_flag) {
            start_timer(t1_read);
        }

        // read reads from file
        num_reads = fastq_fread_batch_max_size(fastq_batch_p, fastq_batch_reader_p->batch_size, fastq_batch_reader_p->fastq_file_p);
        total_reads += num_reads;

        // if there is no reads exit the loop
        if (num_reads == 0) {
            fastq_batch_free(fastq_batch_p);
            break;
        }

        // otherwise, create a new batch object
        item_p = (fastq_batch_list_item_t*) malloc(sizeof(fastq_batch_list_item_t));
        item_p->batch_p = fastq_batch_p;
        item_p->id = num_batchs;

        // fastq_batch_print(fastq_batch_p, reader_fd);

        // insert this batch to the list (synchronization is implemented in the structure, see list.c)
        fastq_batch_list_insert(item_p, fastq_batch_reader_p->batch_list_p);

        sprintf(log_message, "Thread-READ: ....reading batch %i (%i reads), done !!!!\n", num_batchs, num_reads);
        LOG_DEBUG(log_message);

        num_batchs++;

        if (verbose) {
            stop_timer(t1_read, t2_read, read_time);
        }

    } // end of batch loop

    fastq_batch_reader_set_eof(fastq_batch_reader_p, 1);
    fastq_batch_reader_set_alive(fastq_batch_reader_p, 0);
    fastq_batch_list_decr_producers(fastq_batch_reader_p->batch_list_p);

    sprintf(log_message, "Thread-READ: END (%i reads)\n", total_reads);
    LOG_DEBUG(log_message);

    pthread_exit((void*)total_reads);
}
*/