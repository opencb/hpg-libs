
extern "C" {
#include "gpu.h"
}

#include "BW_gpu.cuh"

//-----------------------------------------------------------------------------
// GPU context functions
//-----------------------------------------------------------------------------

gpu_context_t *gpu_context_new(unsigned int id, unsigned int num_threads,
			       bwt_index_t *index) {
  gpu_context_t *context = (gpu_context_t*) calloc(1, sizeof(gpu_context_t));

  // initializes
  context->id = id;
  context->num_threads = num_threads;

  context->d_nWe = NULL; 
  context->d_nWe_size = 0;

  context->d_We = NULL;
  context->d_We_size = 0;

  context->d_k = NULL;
  context->d_k_size = 0;

  context->d_l = NULL;
  context->d_l_size = 0;

  // cudaSetDevice
  gpu_set_device(id);

  // copyVectorGPU
  gpu_copy_vector(&context->d_C,   &index->h_C,   sizeof(size_t));
  gpu_copy_vector(&context->d_C1,  &index->h_C1,  sizeof(size_t));
  gpu_copy_vector(&context->d_rC,  &index->h_rC,  sizeof(size_t));
  gpu_copy_vector(&context->d_rC1, &index->h_rC1, sizeof(size_t));

  // copyCompMatrixGPU
  gpu_copy_matrix(&context->d_O, &index->h_O);
  reverseStrandO(&context->d_rO, &context->d_O);

  return context;
}

//-----------------------------------------------------------------------------

void gpu_context_free(gpu_context_t *context) {
  if (context != NULL) return;

  // free memory                                                  
  if (context->d_We != NULL) {
    cudaFree(context->d_We);
    //manageCudaError();
  }


  if (context->d_nWe != NULL) {
    cudaFree(context->d_nWe);
    //manageCudaError();
  }

  if (context->d_k != NULL) {
    cudaFree(context->d_k);
    //manageCudaError();
  }

  if (context->d_l != NULL) {
    cudaFree(context->d_l);
    //manageCudaError();
  }

  free(context);
}

//-----------------------------------------------------------------------------
// set device
//-----------------------------------------------------------------------------

extern "C" void gpu_set_device(unsigned int id) {
     cudaSetDevice(id);
}

//-----------------------------------------------------------------------------
// copỳ vector and matrix to device
//-----------------------------------------------------------------------------

extern "C" void gpu_copy_vector(vector *dest, vector *src, size_t size) {
  copyVectorGPU(dest, src, size);
}

extern "C" void gpu_copy_matrix(comp_matrix *dest, comp_matrix *src) {
  copyCompMatrixGPU(dest, src);
}

//-----------------------------------------------------------------------------
// reallocate memory in device
//-----------------------------------------------------------------------------

extern "C" void gpu_reallocate_memory(size_t *new_size, size_t *old_size, void **data) {
  if (*old_size < *new_size) {
    if (*data != NULL) {
      cudaFree(*data);
      //manageCudaError();
    }
    cudaMalloc((void**) data, *new_size);
    //manageCudaError();
    *old_size = *new_size;
  } 
}

//-----------------------------------------------------------------------------
// gpu get k and l values
//-----------------------------------------------------------------------------

extern "C" void gpu_get_kl_values(size_t seed_size, size_t min_seed_size, 
       	   			  size_t num_max_seeds,
				  size_t num_reads, size_t seqs_size,
				  char *seqs, size_t *indices,
				  gpu_context_t *context, 
				  size_t *k_values, size_t *l_values,
				  unsigned char mode) {
  
  unsigned int num_threads = context->num_threads;
  unsigned int num_blocks = (num_reads / num_threads) + ((num_reads % num_threads == 0) ? 0 : 1);
  
  // allocate memory in device (GPU) for inputs
  gpu_reallocate_memory(&seqs_size, &context->d_We_size, (void **) &context->d_We);

  size_t indices_size = (num_reads + 1) * sizeof(size_t);
  gpu_reallocate_memory(&indices_size, &context->d_nWe_size, (void **) &context->d_nWe);
  size_t kl_size;

  if (mode == SEED_MODE) {
    kl_size = num_reads * (2 * num_max_seeds) * sizeof(size_t);
  }else {  
    kl_size = num_reads * sizeof(size_t);
  }

  //printf("KL %i\n", num_reads * (2 * num_max_seeds));

  gpu_reallocate_memory(&kl_size, &context->d_k_size, (void **) &context->d_k);
  gpu_reallocate_memory(&kl_size, &context->d_l_size, (void **) &context->d_l);
  
  // copy from host to device (CPU -> GPU)
  cudaMemcpy(context->d_We, seqs, seqs_size, cudaMemcpyHostToDevice);
  //manageCudaError();
  cudaMemcpy(context->d_nWe, indices, indices_size, cudaMemcpyHostToDevice);
  //manageCudaError();

  // call CUDA kernels and copy gpu results to cpu and insert them to list to be processed    

  // searching with d_O (normal strand)
  if (mode == SEED_MODE) {
     //printf("SEED MODE :\n Num Reads: %i\n Seed Size: %i\n Min Seed Size: %i\n Num Max Seeds: %i\n", num_reads, seed_size, min_seed_size, num_max_seeds);

     BWExactSearchBackwardGPUSeedsWrapper_ex(num_blocks, num_threads, context->d_We, context->d_nWe, 
     					context->d_k, context->d_l, 0, context->d_O.siz-2, 
				        &context->d_C, &context->d_C1, &context->d_O, num_reads, 
					seed_size, min_seed_size, num_max_seeds);
  }else {
     BWExactSearchBackwardGPUWrapper_ex(num_blocks, num_threads, context->d_We, context->d_nWe, 
     					context->d_k, context->d_l, 0, context->d_O.siz-2, 
				        &context->d_C, &context->d_C1, &context->d_O, num_reads);
  }

  // copy back the results (GPU -> CPU);
  cudaMemcpy(k_values, context->d_k, kl_size, cudaMemcpyDeviceToHost);
  //manageCudaError();
  cudaMemcpy(l_values, context->d_l, kl_size, cudaMemcpyDeviceToHost);
  //manageCudaError();


  // searching with d_rO (reverse strand)
  if (mode == SEED_MODE) {
     //printf("SEED MODE\n");
     BWExactSearchForwardGPUSeedsWrapper_ex(num_blocks, num_threads, context->d_We, 
					    context->d_nWe,	
					    context->d_k, context->d_l, 0, 
					    context->d_O.siz-2, &context->d_rC, 
					    &context->d_rC1, &context->d_rO, num_reads, 
					    seed_size, min_seed_size, num_max_seeds);
     // copy back the results (GPU -> CPU);
     cudaMemcpy(k_values + (num_reads * (2 * num_max_seeds)), context->d_k, kl_size, cudaMemcpyDeviceToHost);
     //manageCudaError();
     cudaMemcpy(l_values + (num_reads * (2 * num_max_seeds)), context->d_l, kl_size, cudaMemcpyDeviceToHost);
     //manageCudaError();
     
  }else {
     BWExactSearchForwardGPUWrapper_ex(num_blocks, num_threads, context->d_We, context->d_nWe,	
     				       context->d_k, context->d_l, 0, context->d_O.siz-2, 
                                       &context->d_rC, &context->d_rC1, &context->d_rO, num_reads);
     // copy back the results (GPU -> CPU);
     cudaMemcpy(k_values + num_reads, context->d_k, kl_size, cudaMemcpyDeviceToHost);
     //manageCudaError();
     cudaMemcpy(l_values + num_reads, context->d_l, kl_size, cudaMemcpyDeviceToHost);
     //manageCudaError();
     
  }

}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
