#ifndef GPU
#define GPU

#ifdef __cplusplus
extern "C" {
#endif

#include "BW_io.h"
#include "bwt_gpu.h"


//-----------------------------------------------------------------------------
// main structures
//-----------------------------------------------------------------------------

typedef struct gpu_context {
    // device id
  unsigned int id;

  // number of threads
  unsigned int num_threads;
  
  // BWT indices
  vector d_C, d_C1, d_rC, d_rC1;
  comp_matrix d_O, d_rO;

  // allocated GPU  memory and size
  size_t d_We_size;
  size_t d_nWe_size;
  size_t d_k_size;
  size_t d_l_size;

  char *d_We;
  unsigned int *d_nWe;
  size_t *d_k;
  size_t *d_l;
} gpu_context_t;


gpu_context_t *gpu_context_new(unsigned int id, unsigned int num_threads,
			       bwt_index_t *index);
void gpu_context_free(gpu_context_t *context);

//-----------------------------------------------------------------------------
// set device
//-----------------------------------------------------------------------------

void gpu_set_device(unsigned int id);

//-----------------------------------------------------------------------------
// copỳ vector and matrix to device
//-----------------------------------------------------------------------------

void gpu_copy_vector(vector *dest, vector *src, size_t size);
void gpu_copy_matrix(comp_matrix *dest, comp_matrix *src);

//-----------------------------------------------------------------------------
// reallocate memory in device
//-----------------------------------------------------------------------------

void gpu_reallocate_memory(size_t *new_size, size_t *old_size, void **data);

//-----------------------------------------------------------------------------
// gpu get k and l values
//-----------------------------------------------------------------------------

void gpu_get_kl_values(size_t num_reads, size_t seqs_size,
		       char *seqs, size_t *indices,
		       gpu_context_t *context, 
		       size_t *k_values, size_t *l_values);
  

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

#endif // GPU



