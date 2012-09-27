
#ifndef BW_SERVER_H
#define BW_SERVER_H

#include "buffers.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

//------------------------------------------------------------------------------------

void bw_server_input_init(list_t* read_list_p, list_t* kl_list_p, pthread_mutex_t* lock_p, gpu_t* gpu_p, int num_gpus, int num_threads, int mode, bw_server_input_t* input_p);

//------------------------------------------------------------------------------------
// main function
//------------------------------------------------------------------------------------

void bw_server_cuda(bw_server_input_t* input_p);
void bwt_server_omp(bwt_server_input_omp_t* input_p);
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif  // BW_SERVER_H
