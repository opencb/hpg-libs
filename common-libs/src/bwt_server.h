#ifndef BWT_SERVER_H
#define BWT_SERVER_H

#include "commons/commons.h"
#include "containers/array_list.h"

#include "buffers.h"
//#include "fastq_ex_batch.h"
//#include "bam_file.h"
//#include "timing.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct bwt_server_input {
  //int mode;
  size_t batch_size;
  size_t write_size;
  bwt_optarg_t *bwt_optarg;
  list_t* read_list;
  list_t* pair_list;
  list_t* write_list;
  list_t* unmapped_read_list;
  bwt_index_t *bwt_index;
} bwt_server_input_t;


//------------------------------------------------------------------------------------

void bwt_server_input_init(list_t *read_list, unsigned int batch_size, 
			   bwt_optarg_t *bwt_optarg, bwt_index_t *bwt_index, 
			   list_t *pair_list, list_t *write_list, size_t write_size,  
			   list_t* unmapped_read_list, bwt_server_input_t* input);

//------------------------------------------------------------------------------------
// main function
//------------------------------------------------------------------------------------


void bwt_server_cpu(bwt_server_input_t* input);

/*
//void bwt_server_cuda(bwt_server_input_t* input_p);
void server_read_kl(list_t* kl_list_p);
void server_read_split(list_t* split_list_p);
void server_cal_seeker(list_t* kl_list_p);*/
//void bw_server_input_init(list_t* read_list_p, list_t* kl_list_p, pthread_mutex_t* lock_p, gpu_t* gpu_p, int num_gpus, int num_threads, int mode, bw_server_input_t* input_p);
//void bw_server_input_init_omp(list_t* read_list_p, list_t* kl_list_p, int mode, bw_server_input_omp_t* input_p);

//====================================================================================
// apply_bwt
//====================================================================================

void apply_bwt(bwt_server_input_t* input, aligner_batch_t *batch);

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif  // BWT_SERVER_H
