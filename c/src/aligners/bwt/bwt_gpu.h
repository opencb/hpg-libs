#ifndef BWT_GPU_H
#define BWT_GPU_H

#include "bwt.h"
#include "gpu.h"

#define NORMAL_MODE 0
#define SEED_MODE 1

//-----------------------------------------------------------------------------
// general functions
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------

/*size_t bwt_map_exact_seed_batch_gpu(fastq_batch_t *batch,
                                    bwt_optarg_t *bwt_optarg,
                                    cal_optarg_t *cal_optarg,
                                    bwt_index_t *index,
                                    gpu_context_t *gpu_context,
                                    array_list_t **mapping_list);


size_t bwt_map_exact_batch_gpu(fastq_batch_t *batch,
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index,
			       gpu_context_t *gpu_context,
			       fastq_batch_t *unmapped_batch,
			       array_list_t *mapping_list);

size_t bwt_map_exact_seeds_seq_gpu(char *seq, size_t seed_size, 
				   size_t min_seed_size,
				   bwt_optarg_t *bwt_optarg, 
				   bwt_index_t *index, 
				   array_list_t *mapping_list);
*/

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#endif // BWT_GPU_H
