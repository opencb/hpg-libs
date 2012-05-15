
#include "bam_data_batch_list.h"
#include "bam_data_batch.h"
#include "bam_reader.h"

#include "bam.h"
#include "bam_commons.h"
#include "bam_qc_batch.h"
#include "chrom_alignments.h"
#include "sort_thrust.h"



//=====================================================
// functions to manage bam read server
//
//=====================================================

//----------------------------------------------
// bam reader new (sequential & chromosome)
//--------------------------------------------

bam_reader_t* bam_reader_new(char* filename, size_t batch_size, int base_quality, alignments_list_t* alignments_list_p, int mode, int sort, int chromosome) {
  
  bam_reader_t* reader_p = (bam_reader_t*) malloc(sizeof(bam_reader_t));
  
  // open the input file....
  //
  reader_p->bam_file_p = bam_fopen(filename);	// read mode by default
  reader_p->mode = mode;
  reader_p->sort = sort;
  reader_p->chromosome = chromosome;
  reader_p->max_estimated_alignments = get_max_estimated_alignments_by_chromosome(filename);
  reader_p->batch_size = batch_size;
  reader_p->base_quality = base_quality;
    
  reader_p->alive = 1;
  //reader_p->alive_lock = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_init(&(reader_p->alive_lock), NULL);

  reader_p->alignments_list_p = alignments_list_p;
  reader_p->bam_data_batch_list_p = NULL;
 
  return reader_p;
}

//----------------------------------------------
// bam reader new (insert into list of batches)
//--------------------------------------------

bam_reader_t* bam_reader_by_batch_new(char * filename, size_t batch_size, int base_quality, list_t* bam_data_batch_list_p, int mode) {
  
  bam_reader_t* reader_p = (bam_reader_t*) malloc(sizeof(bam_reader_t));
  
  // open the input file....
  //
  reader_p->bam_file_p = bam_fopen(filename);	// read mode by default
  reader_p->mode = mode;
  reader_p->sort = -1;
  reader_p->chromosome = -1;
  reader_p->max_estimated_alignments = get_max_estimated_alignments_by_chromosome(filename);
  reader_p->batch_size = batch_size;
  reader_p->base_quality = base_quality;
    
  reader_p->alive = 1;
  //reader_p->alive_lock = PTHREAD_MUTEX_INITIALIZER;
  pthread_mutex_init(&(reader_p->alive_lock), NULL);
  
  reader_p->alignments_list_p = NULL;
  reader_p->bam_data_batch_list_p = bam_data_batch_list_p;
  
  //bam_data_batch_list_incr_producers(bam_data_batch_list_p);
  //list_incr_writers(bam_data_batch_list_p);		//initialized when creating the list
  
  return reader_p;
}

//----------------------------------------------
// bam reader free
//--------------------------------------------

void bam_reader_free(bam_reader_t* reader_p) {

  // close the input file and exiting....
  //
  //bam_file_close(reader_p->bam_file_p);
  list_item_t* item_p;

  if (reader_p->mode == LIST_INSERT_MODE) {
    //bam_data_batch_list_items_free(reader_p->bam_data_batch_list_p);
    for (int i = 0; i < reader_p->bam_data_batch_list_p->length; i++) {
      item_p = list_remove_item(reader_p->bam_data_batch_list_p);
      list_item_free(item_p);
    }
  } else {
    free(reader_p->alignments_list_p);
  }

  free(reader_p);
}

//----------------------------------------------
// bam reader start
//----------------------------------------------

void bam_reader_start(bam_reader_t* reader_p) {

  // create and launch pthread....
  //
  //printf("launching bam read server thread....\n");

  if (reader_p->mode == CHROMOSOME_MODE) {
      pthread_create(&(reader_p->thread), NULL, bam_reader_chromosome_thread_function, (void*) reader_p);
  } else if (reader_p->mode == SEQUENTIAL_MODE) {
      pthread_create(&(reader_p->thread), NULL, bam_reader_sequential_thread_function, (void*) reader_p);
  } else if (reader_p->mode == LIST_INSERT_MODE) {
      pthread_create(&(reader_p->thread), NULL, bam_reader_list_insert_thread_function, (void*) reader_p);
      //pthread_create(&(reader_p->thread), NULL, bam_reader_list_insert_by_chromosome_thread_function, (void*) reader_p);
  }
}

//----------------------------------------------
// bam reader join
//----------------------------------------------

unsigned int bam_reader_join(bam_reader_t* reader_p) {
  void* r;
  pthread_join(reader_p->thread, &r);
  return (uintptr_t) r;
}

//----------------------------------------------
// 'secure' set and get for alive
//--------------------------------------------

void bam_reader_set_alive(bam_reader_t* reader_p, int alive) {
  pthread_mutex_lock(&(reader_p->alive_lock));
  reader_p->alive = alive;
  pthread_mutex_unlock(&(reader_p->alive_lock));
}

int bam_reader_get_alive(bam_reader_t* reader_p) {
  int alive;
  pthread_mutex_lock(&(reader_p->alive_lock));
  alive = reader_p->alive;
  pthread_mutex_unlock(&(reader_p->alive_lock));
  return alive;
}

//-----------------------------------------------------
// bam reader thread function
//
// this thread reads aligments from a BAM file in disk,
// and save them classifiend by chromosome to be sorted 
// by the gpu
// 
//-----------------------------------------------------

//----------------------------------------------
// reader with one file pass per chromosome
//----------------------------------------------

void* bam_reader_sequential_thread_function(void* param_p) {

  bam_reader_t* reader_p = (bam_reader_t*) param_p;

  bam_reader_set_alive(reader_p, 1);
  
  char log_message[200];
  sprintf(log_message, "Thread-READ: START, for file %.150s and chromosome %i\n", reader_p->bam_file_p->filename, reader_p->chromosome);
  LOG_DEBUG(log_message);
 
  int num_alignments = 0, total_alignments = 0;
  int current_chromosome = 0;
  bam_batch_t* batch_p;
  
  if (reader_p->chromosome == ALL_CHROMOSOMES) {
    current_chromosome = 0;
  } else {
    current_chromosome = reader_p->chromosome;
  }
  
  alignments_list_new_chrom_alignment(current_chromosome, reader_p->max_estimated_alignments, reader_p->alignments_list_p);

  while (1) {
    //printf("Thread-READ: Before reading bam batch for chromosome %i...\n", current_chromosome);
    
    batch_p = bam_batch_new(reader_p->batch_size, SINGLE_CHROM_BATCH);
         
    if (time_flag) { start_timer(t1_read); }
    
    num_alignments = bam_fread_max_size_by_chromosome(batch_p, reader_p->batch_size, current_chromosome, reader_p->bam_file_p);
    total_alignments += num_alignments;

    if (time_flag) { stop_timer(t1_read, t2_read, read_time); }

    //printf("chromosome %i : num_alignments: %i, allocated alignments = %i\n", current_chromosome, num_alignments, batch_p->allocated_alignments);     

    if (batch_p->num_alignments == 0) {

        bam_fclose(reader_p->bam_file_p);
	reader_p->bam_file_p = bam_fopen(reader_p->bam_file_p->filename);	
	
	// sorting step
	if (time_flag) { start_timer(t1_sort); }	
	if (reader_p->sort == SORT_BY_POSITION) {
	  sort_alignments_by_position(reader_p->alignments_list_p, current_chromosome);
	} else if (reader_p->sort == SORT_BY_ID) {
	  
	}
	
	if (time_flag) { stop_timer(t1_sort, t2_sort, sort_time); }
	
	current_chromosome++;	
	bam_batch_free(batch_p, 1);

	if ((reader_p->chromosome == ALL_CHROMOSOMES) && (current_chromosome == NUM_OF_CHROMOSOMES)) { 
	  break;
	} else if ((reader_p->chromosome != ALL_CHROMOSOMES) && (current_chromosome == (reader_p->chromosome + 1))) { 
	  break;
	}
	
	alignments_list_new_chrom_alignment(current_chromosome, reader_p->max_estimated_alignments, reader_p->alignments_list_p);
	
    } else {
	alignments_list_insert_batch(batch_p, reader_p->alignments_list_p);

	bam_batch_free(batch_p, 0);  
    }   

    while (get_free_memory() < MIN_FREE_MEMORY_ALLOWED_FOR_BAM_READ_KB)  {
      //printf("Thread-READ: go to sleep for a while...\n");
      sched_yield();      
      usleep(10000);      
    }

  }

  bam_fclose(reader_p->bam_file_p);

  if (reader_p->chromosome == ALL_CHROMOSOMES) {
    for(int i=0; i<NUM_OF_CHROMOSOMES; i++) {
      chrom_alignments_set_complete(alignments_list_get_chrom_alignment(i, reader_p->alignments_list_p), 1);
    }
  } else {
    chrom_alignments_set_complete(alignments_list_get_chrom_alignment(reader_p->chromosome, reader_p->alignments_list_p), 1);
  }
 
  bam_reader_set_alive(reader_p, 0);
  bam_reader_alive=0;
  
  LOG_DEBUG("Thread-READ: END \n");

  pthread_exit((void*)total_alignments);
}

//----------------------------------------------
// reader with one file pass per chromosome
//----------------------------------------------

void* bam_reader_chromosome_thread_function(void* param_p) {
 
  bam_reader_t* reader_p = (bam_reader_t*) param_p;

  bam_reader_set_alive(reader_p, 1);
  
  char log_message[200];
  sprintf(log_message, "Thread-READ: START, for file %.150s and chromosome %i\n", reader_p->bam_file_p->filename, reader_p->chromosome);
  LOG_DEBUG(log_message);
 
  int num_alignments = 0, total_alignments = 0;
  int current_chromosome = 0;
  bam_batch_t* batch_p;
  
  for (int i=0; i<NUM_OF_CHROMOSOMES; i++) {
    alignments_list_new_chrom_alignment(i, reader_p->max_estimated_alignments, reader_p->alignments_list_p);
  }

  while (1) {

    //printf("Thread-READ: Before reading bam batch...\n");

    batch_p = bam_batch_new(reader_p->batch_size, MULTIPLE_CHROM_BATCH);
         
    if (time_flag) { start_timer(t1_read); }
    
    num_alignments = bam_fread_max_size(batch_p, reader_p->batch_size, reader_p->base_quality, reader_p->bam_file_p);
    total_alignments += num_alignments;
    //printf("chromosome %i : num_alignments: %i, allocated alignments = %i\n", current_chromosome, num_alignments, batch_p->allocated_alignments);     
    
    if (time_flag) { stop_timer(t1_read, t2_read, read_time); }
 
    if (batch_p->num_alignments == 0) {
    
	// sorting step
	//once for each chromosome
	//
	if (time_flag) { start_timer(t1_sort); }
	
	if (reader_p->sort == SORT_BY_POSITION) {
	  sort_alignments_by_position(reader_p->alignments_list_p, reader_p->chromosome);
	}
	
	if (time_flag) { stop_timer(t1_sort, t2_sort, sort_time); }
	
	bam_batch_free(batch_p, 1);
	break;

    } else {
	
	alignments_list_insert_batch(batch_p, reader_p->alignments_list_p);
	
	bam_batch_free(batch_p, 0);
    }   

    while (get_free_memory() < MIN_FREE_MEMORY_ALLOWED_FOR_BAM_READ_KB)  {
      //printf("Thread-READ: go to sleep for a while...\n");
      sched_yield();      
      usleep(10000);
    }
    
  }

  bam_fclose(reader_p->bam_file_p);

  for(int i=0; i<NUM_OF_CHROMOSOMES; i++) {
    chrom_alignments_set_complete(alignments_list_get_chrom_alignment(i, reader_p->alignments_list_p), 1);
  }
  
  bam_reader_set_alive(reader_p, 0);
  bam_reader_alive=0;
  
  LOG_DEBUG("Thread-READ: END \n");

  pthread_exit((void*)total_alignments);
}

//----------------------------------------------
// reader with batch insertion in list
//----------------------------------------------

void* bam_reader_list_insert_thread_function(void* param_p) {
  
  if (time_flag) { start_timer(t1_reader_server); }
  if (time_flag) { start_timer(t1_active_reader); }
 
  bam_reader_t* reader_p = (bam_reader_t*) param_p;

  bam_reader_set_alive(reader_p, 1);
  
  char log_message[200];
  sprintf(log_message, "Thread-READ: START, for file %.150s\n", reader_p->bam_file_p->filename);
  LOG_DEBUG(log_message);
  
  int num_alignments = 0, total_alignments = 0;
  int current_chromosome = 0;
  int counter = 0;
  bam_batch_t* bam_batch_p;
  bam_data_batch_t* bam_data_batch_p;
  //bam_data_batch_list_item* bam_data_batch_list_item_p;
  list_item_t* bam_data_batch_list_item_p;
  uint8_t* prev_seq = NULL;
  int* prev_seq_length = NULL;
  int* prev_seq_start_coordinate = NULL;
 
  while (1) {
   
    //printf("Thread-READ: Before reading bam batch...\n");
    
    bam_batch_p = bam_batch_new(reader_p->batch_size, MULTIPLE_CHROM_BATCH); 
        
    if (time_flag) { start_timer(t1_read); }
    
    num_alignments = bam_fread_max_size_no_duplicates(bam_batch_p, reader_p->batch_size, reader_p->base_quality, reader_p->bam_file_p, prev_seq, prev_seq_length, prev_seq_start_coordinate);
    total_alignments += num_alignments;
    //printf("bam_batch_p->allocated_alignments: %i, bam_batch_p->num_alignments: %i\n", bam_batch_p->allocated_alignments, bam_batch_p->num_alignments);
    
    if (time_flag) { stop_timer(t1_read, t2_read, read_time); }

    // exit condition, no more alignments in bam file
    //
    if (num_alignments == 0) {
      bam_batch_free(bam_batch_p, 1);
      break;
    }

    //if (time_flag) { start_timer(t1_read); }    
    bam_data_batch_p = bam_data_batch_new(bam_batch_p->num_alignments); 	//bam_batch_p->num_alignments is still 0!!!!
    bam_data_batch_p = bam_data_batch_init(bam_data_batch_p, bam_batch_p);
    //if (time_flag) { stop_timer(t1_read, t2_read, read_time); }
    
    //if (time_flag) { start_timer(t1_read); }
    //bam_data_batch_list_item_p = bam_data_batch_list_item_new(counter, num_alignments, bam_batch_p, bam_data_batch_p);     
    //bam_data_batch_list_insert(bam_data_batch_list_item_p, reader_p->bam_data_batch_list_p);
    bam_data_batch_list_item_p = list_item_new(counter, 0, bam_data_batch_p);
    list_insert_item(bam_data_batch_list_item_p, reader_p->bam_data_batch_list_p);    
    
    //if (time_flag) { stop_timer(t1_read, t2_read, read_time); }
    //bam_data_batch_list_print(reader_p->bam_data_batch_list_p);

    bam_batch_free(bam_batch_p, 1);
    
    //if (time_flag) { stop_timer(t1_read, t2_read, read_time); }
    
    while (get_free_memory() < MIN_FREE_MEMORY_ALLOWED_FOR_BAM_READ_KB)  {
      printf("Thread-READ: go to sleep for a while...\n");
      sched_yield();      
      usleep(10000);
    }
    
    counter++;    
  }
  
  bam_fclose(reader_p->bam_file_p);
  
  //bam_data_batch_list_decr_producers(reader_p->bam_data_batch_list_p);
  list_decr_writers(reader_p->bam_data_batch_list_p);
  bam_reader_set_alive(reader_p, 0);
  bam_reader_alive=0;
  
  if (time_flag) { stop_timer(t1_reader_server, t2_reader_server, reader_server_time); }
  LOG_DEBUG("Thread-READ: END \n");

  pthread_exit((void*)total_alignments);
}

//----------------------------------------------
// reader with batch insertion in list
// 
// sequential read, one pass, each batch only
// contains data from a single chromosome
//----------------------------------------------

void* bam_reader_list_insert_by_chromosome_thread_function(void* param_p) {
  
  if (time_flag) { start_timer(t1_reader_server); }
  
  bam_reader_t* reader_p = (bam_reader_t*) param_p;

  bam_reader_set_alive(reader_p, 1);
  
  char log_message[200];
  sprintf(log_message, "Thread-READ: START, for file %.150s\n", reader_p->bam_file_p->filename);		
  LOG_DEBUG(log_message);
 
  int num_alignments = 0, total_alignments = 0;
  int current_chromosome = 0;
  int counter = 0;
  bam_batch_t* bam_batch_p;
  bam_data_batch_t* bam_data_batch_p;
  //bam_data_batch_list_item* bam_data_batch_list_item_p;
  list_item_t* bam_data_batch_list_item_p;

  while (1) {

    //printf("Thread-READ: Before reading bam batch...\n");
    
    bam_batch_p = bam_batch_new(reader_p->batch_size, MULTIPLE_CHROM_BATCH); 
        
    if (time_flag) { start_timer(t1_read); }
    
    num_alignments = bam_fread_max_size(bam_batch_p, reader_p->batch_size, reader_p->base_quality, reader_p->bam_file_p);
    total_alignments += num_alignments;
    //printf("bam_batch_p->allocated_alignments: %i, bam_batch_p->num_alignments: %i\n", bam_batch_p->allocated_alignments, bam_batch_p->num_alignments);
    
    if (time_flag) { stop_timer(t1_read, t2_read, read_time); }

    // exit condition, no more alignments in bam file
    //
    if (num_alignments == 0) {
      break;
    }
    
    int start_alignment = 0;
   
    while (1) {
 
      bam_data_batch_p = bam_data_batch_new(bam_batch_p->num_alignments);
      bam_data_batch_p = bam_data_batch_by_chromosome_init(bam_data_batch_p, bam_batch_p, start_alignment);

      start_alignment += bam_data_batch_p->num_alignments;

      if (start_alignment == num_alignments) break;
      
      //bam_data_batch_list_item_p = bam_data_batch_list_item_new(counter, bam_data_batch_p->num_alignments, bam_batch_p, bam_data_batch_p);     
      //bam_data_batch_list_insert(bam_data_batch_list_item_p, reader_p->bam_data_batch_list_p);
      bam_data_batch_list_item_p = list_item_new(counter, 0, bam_data_batch_p);
      list_insert_item(bam_data_batch_list_item_p, reader_p->bam_data_batch_list_p);
      
      counter++;    
    }
    //bam_data_batch_list_print(reader_p->bam_data_batch_list_p);

    bam_batch_free(bam_batch_p, 1);
    //bam_data_batch_free(bam_data_batch_p);

    while (get_free_memory() < MIN_FREE_MEMORY_ALLOWED_FOR_BAM_READ_KB)  {
      //printf("Thread-READ: go to sleep for a while...\n");
      sched_yield();      
      usleep(10000);
    }

  }
      
  bam_fclose(reader_p->bam_file_p);
  
  //bam_data_batch_list_decr_producers(reader_p->bam_data_batch_list_p);
  list_decr_writers(reader_p->bam_data_batch_list_p);
  bam_reader_set_alive(reader_p, 0);
  bam_reader_alive=0;

  if (time_flag) { stop_timer(t1_reader_server, t2_reader_server, reader_server_time); }
  LOG_DEBUG("Thread-READ: END \n");

  pthread_exit((void*)total_alignments);
}

/**
 *  @brief Generic BAM reader
 *  @param param_p structure with reader parameters *  
 *  @return void
 *  
 *  Generic BAM reader fills a bam batch list with batches of the indicated size containing bam1_t alignments
 */

// void* bam_reader_thread_function(void* param_p) {  
//   if (time_flag) { start_timer(t1_reader_server); }
//   if (time_flag) { start_timer(t1_active_reader); }
//  
//   bam_reader_t* reader_p = (bam_reader_t*) param_p;
// 
//   bam_reader_set_alive(reader_p, 1);
//   
//   char log_message[200];
//   sprintf(log_message, "Thread-READ: START, for file %.150s\n", reader_p->bam_file_p->filename);
//   LOG_DEBUG(log_message);
//   
//   int num_alignments = 0, total_alignments = 0;
//   int current_chromosome = 0;
//   int counter = 0;
//   bam_batch_t* bam_batch_p;
//   list_item_t* bam_batch_list_item_p;
//   uint8_t* prev_seq = NULL;
//   int* prev_seq_length = NULL;
//   int* prev_seq_start_coordinate = NULL;
//  
//   while (true) {
//    
//     //printf("Thread-READ: Before reading bam batch...\n");
//     
//     bam_batch_p = bam_batch_new(reader_p->batch_size, MULTIPLE_CHROM_BATCH); 
//         
//     if (time_flag) { start_timer(t1_read); }
//     
//     num_alignments = bam_fread_max_size_no_duplicates(bam_batch_p, reader_p->batch_size, reader_p->base_quality, reader_p->bam_file_p, prev_seq, prev_seq_length, prev_seq_start_coordinate);
//     total_alignments += num_alignments;
//     
//     if (time_flag) { stop_timer(t1_read, t2_read, read_time); }
// 
//     //exit condition, no more alignments in bam file
//     if (num_alignments == 0) {
//       break;
//     }
//     
//     //if (time_flag) { start_timer(t1_read); }
//     bam_batch_list_item_p = list_item_new(counter, 0, bam_batch_p);
//     list_insert_item(bam_batch_list_item_p, reader_p->bam_batch_list_p);
//     //if (time_flag) { stop_timer(t1_read, t2_read, read_time); }
//     
//     while (get_free_memory() < MIN_FREE_MEMORY_ALLOWED_FOR_BAM_READ_KB)  {
//       //printf("Thread-READ: go to sleep for a while...\n");
//       sched_yield();      
//       usleep(10000);
//     }
//     
//     counter++;    
//   }
//   
//   bam_fclose(reader_p->bam_file_p);
//   
//   //bam_data_batch_list_decr_producers(reader_p->bam_data_batch_list_p);
//   list_decr_writers(reader_p->bam_batch_list_p);
//   bam_reader_set_alive(reader_p, 0);
//   bam_reader_alive=0;
//   
//   if (time_flag) { stop_timer(t1_reader_server, t2_reader_server, reader_server_time); }
//   LOG_DEBUG("Thread-READ: END \n");
// 
//   pthread_exit((void*)total_alignments);
// }