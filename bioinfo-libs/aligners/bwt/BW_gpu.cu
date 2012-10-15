/*
    bwa_gpu a set of tools which allow short sequence alignment using the Burrows-Wheeler
    transform usign both CPU and GPU approaches.
    Copyright (C) 2011  Jose Salavert Torres, Ignacio Blanquer Espert,
                        Andres Tomas Dominguez, Vicente Hernandez Garcia,
	 		Ignacio Medina Castello, Joaquin Tarraga Gimenez,
			Joaquin Dopazo Blazquez

    Contact e-mail: josator@fiv.upv.es, iblanque@dsic.upv.es, atomas@dsic.upv.es,
                    vhernand@dsic.upv.es, imedina@cipf.es, jtarraga@cipf.es,
                    jdopazo@cipf.es

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include"BW_io.h"
#include"BW_gpu.cuh"

#include <cuda.h>
//#include "cuPrintf.cu"

#if defined VECTOR_O_32BIT_COMPRESSION

__device__ unsigned int getOcompValueGPU(size_t n, size_t m, comp_matrix O) {

  size_t pos, desp;

  pos  = m / 32;
  desp = m % 32;

  return O.desp[n][pos] + __popc( O.count[n][pos] << (32 - (desp + 1)) );

}

#elif defined VECTOR_O_64BIT_COMPRESSION

__device__ unsigned int getOcompValueGPU(size_t n, size_t m, comp_matrix O) {

  size_t pos, desp;

  pos  = m / 64;
  desp = m % 64;

  return O.desp[n][pos] + __popcll( O.count[n][pos] << (64 - (desp + 1)) );

}

#endif

void readCompMatrixGPU(comp_matrix *matrix, const char *directory, const char *name) {

  size_t err=0;
  cudaError_t error;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".desp");

  fp  = fopen(path,  "rb+");
  checkFileOpen(fp, path);

  err = fread(&matrix->siz,    sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  err = fread(&matrix->n_desp, sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  err = fread(&matrix->m_desp, sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  for (size_t i=0; i<matrix->n_desp; i++) {
    cudaMallocHost((void**) &matrix->desp[i], matrix->m_desp * sizeof(unsigned int));
    manageCudaError();
    err = fread(matrix->desp[i], sizeof(unsigned int), matrix->m_desp, fp);
    checkFileRead(err, matrix->m_desp, path);
  }

  fclose(fp);

#if defined VECTOR_O_32BIT_COMPRESSION

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".count");

  fp  = fopen(path,  "rb+");
  checkFileOpen(fp, path);

  err = fread(&matrix->n_count,   sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  err = fread(&matrix->m_count,   sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  for (size_t i=0; i<matrix->n_count; i++){
    cudaMallocHost((void**) &matrix->count[i], matrix->m_count * sizeof(unsigned int));
    manageCudaError();
    err = fread(matrix->count[i], sizeof(unsigned int), matrix->m_count, fp);
    checkFileRead(err, matrix->m_count, path);
  }

  fclose(fp);

#elif defined VECTOR_O_64BIT_COMPRESSION

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".count");

  fp  = fopen(path,  "rb+");
  checkFileOpen(fp, path);

  err = fread(&matrix->n_count, sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  err = fread(&matrix->m_count, sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  for (size_t i=0; i<matrix->n_count; i++){

    cudaMallocHost((void**) &matrix->count[i], matrix->m_count * sizeof(unsigned long long));
    manageCudaError();
    err = fread(matrix->count[i], sizeof(unsigned long long), matrix->m_count, fp);
    checkFileRead(err, matrix->m_count, path);

  }

  fclose(fp);

#endif

}

void copyVectorGPU(vector *device, vector *host, size_t data_size) {

  cudaError_t error;

  device->n = host->n;
  //  printf("device->n = %d\n", device->n);
  cudaMalloc((void**) &device->vector,  device->n * data_size);
  //  printf("cudaMalloc done !!\n");
  manageCudaError();
  cudaMemcpy(device->vector, host->vector, device->n * data_size, cudaMemcpyHostToDevice);
  //  printf("cudaMemcpy done !!\n");
  manageCudaError();

}

void copyCompMatrixGPU(comp_matrix *device, comp_matrix *host) {

  cudaError_t error;

  device->siz    = host->siz;
  device->n_desp = host->n_desp;
  device->m_desp = host->m_desp;

  for (size_t i=0; i<device->n_desp; i++) {
    cudaMalloc((void**) &device->desp[i], device->m_desp * sizeof(unsigned int));
    manageCudaError();
    cudaMemcpy(device->desp[i], host->desp[i], host->m_desp * sizeof(unsigned int), cudaMemcpyHostToDevice);
    manageCudaError();
  }

#if defined   VECTOR_O_32BIT_COMPRESSION
  
  device->n_count = host->n_count;
  device->m_count = host->m_count;

  for (size_t i=0; i<device->n_count; i++) {
    cudaMalloc((void**) &device->count[i], device->m_count * sizeof(unsigned int));
    manageCudaError();
    cudaMemcpy(device->count[i], host->count[i], host->m_count * sizeof(unsigned int), cudaMemcpyHostToDevice);
    manageCudaError();
  }

#elif defined VECTOR_O_64BIT_COMPRESSION

  device->n_count = host->n_count;
  device->m_count = host->m_count;

  for (size_t i=0; i<device->n_count; i++) {
    cudaMalloc((void**) &device->count[i], device->m_count * sizeof(unsigned long long));
    manageCudaError();
    cudaMemcpy(device->count[i], host->count[i], host->m_count * sizeof(unsigned long long), cudaMemcpyHostToDevice);
    manageCudaError();
  }

#endif
}

void freeCompMatrixGPUHost(comp_matrix *matrix) {

for (size_t i=0; i<matrix->n_desp; i++) {
    cudaFreeHost(matrix->desp[i]);
#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION
    cudaFreeHost(matrix->count[i]);
#endif
  }

}

void freeCompMatrixGPUDevice(comp_matrix *matrix) {

  for (size_t i=0; i<matrix->n_desp; i++) {
    cudaFree(matrix->desp[i]);
#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION
    cudaFree(matrix->count[i]);
#endif
  }

}

__global__ void BWExactSearchBackwardGPU(char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, unsigned int *C, unsigned int *C1, comp_matrix O) {

  int i;
  BWiterationVariablesGPU();
  size_t k2, l2;
  unsigned int offset  = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ unsigned int Cshared[4];
  __shared__ unsigned int C1shared[4];
  
  if (threadIdx.x<4) {
    Cshared[threadIdx.x]  = C[threadIdx.x];
    C1shared[threadIdx.x] = C1[threadIdx.x];
  }

  __syncthreads();

  k2 = k_ini; l2 = l_ini;

  for (i=nW[offset]-1; (k2<=l2) && (i>=0); i--)
    BWiterationGPU(k2, l2, k2, l2, W[offset*MAXLINE+i], Cshared, C1shared, O);

  k[offset] = k2;
  l[offset] = l2;

}

__global__ void BWExactSearchForwardGPU(char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, unsigned int *C, unsigned int *C1, comp_matrix O) {

  int i;
  BWiterationVariablesGPU();
  size_t k2, l2;
  unsigned int offset  = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ unsigned int Cshared[4];
  __shared__ unsigned int C1shared[4];

  if (threadIdx.x<4) {
    Cshared[threadIdx.x] = C[threadIdx.x];
    C1shared[threadIdx.x] = C1[threadIdx.x];
  }

  __syncthreads();

  k2 = k_ini;  l2 = l_ini;

  for (i=0; (k2<=l2) && (i<nW[offset]); i++)
    BWiterationGPU(k2, l2, k2, l2, W[offset*MAXLINE+i], Cshared, C1shared, O);

  k[offset] = k2;
  l[offset] = l2;

}

//-----------------------------------------------------------------------------
// _ex functions, for extended fucntions, 
//
// These functions adds a new parameter num_reads, in addition, the nW vector
// contains indices to the start of each read in W vector, instead of the length
// of the read as the original functions ("no _ex") do
//
//-----------------------------------------------------------------------------

__global__ void BWExactSearchBackwardGPU_ex(char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, unsigned int *C, unsigned int *C1, comp_matrix O, size_t num_reads) {

  int i, len;
  BWiterationVariablesGPU();
  size_t k2, l2, index;
  unsigned int offset  = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ unsigned int Cshared[4];
  __shared__ unsigned int C1shared[4];
  
  if (threadIdx.x<4) {
    Cshared[threadIdx.x]  = C[threadIdx.x];
    C1shared[threadIdx.x] = C1[threadIdx.x];
  }

  __syncthreads();

  if (offset < num_reads) {

    k2 = k_ini; l2 = l_ini;
    index = nW[offset];
    len = nW[offset + 1] - index - 1;

    for (i = len - 1; (k2 <= l2) && (i >= 0); i--) {
      BWiterationGPU(k2, l2, k2, l2, W[index + i], Cshared, C1shared, O);
    }
    k[offset] = k2; //k_ini; //k2;
    l[offset] = l2; //l_ini; //l2;
  }
}

//-----------------------------------------------------------------------------

__global__ void BWExactSearchForwardGPU_ex(char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, unsigned int *C, unsigned int *C1, comp_matrix O, size_t num_reads) {

  int i, len;
  BWiterationVariablesGPU();
  size_t k2, l2, index;
  unsigned int offset  = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ unsigned int Cshared[4];
  __shared__ unsigned int C1shared[4];

  if (threadIdx.x<4) {
    Cshared[threadIdx.x] = C[threadIdx.x];
    C1shared[threadIdx.x] = C1[threadIdx.x];
  }

  __syncthreads();

  if (offset < num_reads) {
    k2 = k_ini;  l2 = l_ini;
    index = nW[offset];
    len = nW[offset + 1] - index - 1;

    for (i = 0; (k2 <= l2) && (i < len); i++) {
      BWiterationGPU(k2, l2, k2, l2, W[index + i], Cshared, C1shared, O);
    }

    k[offset] = k2; //index; // k_ini; //k2;
    l[offset] = l2; //len; // l_ini; //l2;
  }
}

//-----------------------------------------------------------------------------
// end of _ex functions
//-----------------------------------------------------------------------------


//IN PROGRESS:
__global__ void BWExactFinalResultsBackwardGPU(char *W, unsigned int *C, unsigned int *C1, comp_matrix O, results_list rl_prev, results_list rl_next, unsigned int chunk_size, int *stack_size) {

  BWiterationVariablesGPU();
  size_t k, l;
  int start, pos, pos_start, end;
  unsigned read_index, read_offset;
  unsigned int offset  = blockIdx.x * blockDim.x + threadIdx.x;
  result *r_prev, *r_next;

  __shared__ unsigned int Cshared[4];
  __shared__ unsigned int C1shared[4];

  if (offset==0)
    *stack_size=0;

  if (threadIdx.x<4) {
    Cshared[threadIdx.x] = C[threadIdx.x];
    C1shared[threadIdx.x] = C1[threadIdx.x];
  }

  __syncthreads();

  r_prev = &rl_prev.list[offset];

  start      = r_prev->start;
  pos        = r_prev->pos;
  end        = r_prev->end;
  k          = r_prev->k;
  l          = r_prev->l;
  read_index = r_prev->read_index;
  read_offset = /*read_index*/offset*MAXLINE;

  pos_start = pos - chunk_size;
  if (pos_start < start) pos_start = start;

  for(; pos>=pos_start; pos--) {
    BWiterationGPU(k, l, k, l, (size_t)W[read_offset + pos], Cshared, C1shared, O);
    if (k > l) {
      pos=start-1; break;
    }
  }

  r_next = &rl_next.list[/*atomicAdd(stack_size,1)*/offset];

  r_next->start = start;
  r_next->pos = pos;
  r_next->end = end;
  r_next->k = k;
  r_next->l = l;
  r_next->read_index = read_index;

}

__global__ void BWExactSearchBackwardVectorGPU(char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, unsigned int *C, unsigned int *C1, comp_matrix O) {

  int i;
  BWiterationVariablesGPU();
  size_t k2, l2;
  unsigned int offset  = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ unsigned int Cshared[4];
  __shared__ unsigned int C1shared[4];

  if (threadIdx.x<4) {
    Cshared[threadIdx.x] = C[threadIdx.x];
    C1shared[threadIdx.x] = C1[threadIdx.x];
  }

  __syncthreads();

  k2 = k_ini;  l2 = l_ini;

  for (i=nW[offset]-1; i>=0; i--) {

    BWiterationGPU(k2, l2, k2, l2, W[offset*MAXLINE + i], Cshared, C1shared, O);

    k[offset*MAXLINE+i] = k2;
    l[offset*MAXLINE+i] = l2;

    if (k2 > l2) {
      i--;
      break;
    }

  }

  for(;i>=0; i--) {
    k[offset*MAXLINE+i] = k2;
    l[offset*MAXLINE+i] = l2;
  }

}

__global__ void BWExactSearchForwardVectorGPU(char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, unsigned int *C, unsigned int *C1, comp_matrix O) {

  int i;
  BWiterationVariablesGPU();
  size_t k2, l2;
  unsigned int offset  = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ unsigned int Cshared[4];
  __shared__ unsigned int C1shared[4];

  if (threadIdx.x<4) {
    Cshared[threadIdx.x] = C[threadIdx.x];
    C1shared[threadIdx.x] = C1[threadIdx.x];
  }

  __syncthreads();

  k2 = k_ini;  l2 = l_ini;

  for (i=0; i<nW[offset]; i++) {

    BWiterationGPU(k2, l2, k2, l2, W[offset*MAXLINE + i], Cshared, C1shared, O);

    k[offset*MAXLINE+i] = k2;
    l[offset*MAXLINE+i] = l2;

  }

}

/*
__global__ void BWExactIterativeSearchGPU(char *W, int *nW, int *nWe, int *k, int *l, int k_ini, int l_ini, int *C,  int *O,  int sizO) {

  int i, b; //, pos;
  unsigned long int k2, l2;
  char val1, val2, val3, val4;
  int siz1, siz2, siz3, siz4;

  int offset  = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ int Cshared[4];

  if (threadIdx.x<4) { // Minimum data is a 32 block
    Cshared[threadIdx.x] = C[threadIdx.x];
  }

  __syncthreads();

  k2 = k_ini;
  l2 = l_ini;

  // First block of 4 bases should not be fully filled

  int real_size = nWe[offset];
  int resto = real_size % 4;
 
  b = W[offset*MAXLINECOMP + nW[offset] - 1];

  switch(resto) {

  case 0:
    val4 = ( b >> 6 ) & 3;
    siz4 = val4*sizO;
    
    k2 = Cshared[val4] + O[siz4 + k2    ] + 1;
    l2 = Cshared[val4] + O[siz4 + l2 + 1];
    
  case 3:
    val3 = ( b >> 4 ) & 3;
    siz3 = val3*sizO;

    k2 = Cshared[val3] + O[siz3 + k2    ] + 1;
    l2 = Cshared[val3] + O[siz3 + l2 + 1];
    
  case 2:
    val2 = ( b >> 2 ) & 3;
    siz2 = val2*sizO;

    k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
    l2 = Cshared[val2] + O[siz2 + l2 + 1];

  case 1:
    val1 = ( b      ) & 3;
    siz1 = val1*sizO;

    k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
    l2 = Cshared[val1] + O[siz1 + l2 + 1];

  }

  __syncthreads();

  for (i=nW[offset]-2; (k2<=l2) && (i>=0); i--) {

    b = W[offset*MAXLINECOMP + i];

    val4 = ( b >> 6 ) & 3;
    siz4 = val4*sizO;

    k2 = Cshared[val4] + O[siz4 + k2    ] + 1;
    l2 = Cshared[val4] + O[siz4 + l2 + 1];

    val3 = ( b >> 4 ) & 3;
    siz3 = val3*sizO;

    k2 = Cshared[val3] + O[siz3 + k2    ] + 1;
    l2 = Cshared[val3] + O[siz3 + l2 + 1];
    
    val2 = ( b >> 2 ) & 3;
    siz2 = val2*sizO;

    k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
    l2 = Cshared[val2] + O[siz2 + l2 + 1];

    val1 = ( b      ) & 3;
    siz1 = val1*sizO;

    k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
    l2 = Cshared[val1] + O[siz1 + l2 + 1];

  }

  __syncthreads();

   k[offset] = k2;
   l[offset] = l2;

}

__global__ void BWExactIterativeSearchGPURev(char *W, int *nW, int *nWe, int *k, int *l, int k_ini, int l_ini, int *C,  int *O,  int sizO) {

  int i, b;//, pos;
  unsigned long int k2, l2;
  char val1, val2, val3, val4;
  int siz1, siz2, siz3, siz4;

  int offset  = blockIdx.x * blockDim.x + threadIdx.x;

  __shared__ int Cshared[4];

  if (threadIdx.x<4) { // Minimum data is a 32 block
    Cshared[threadIdx.x] = C[threadIdx.x];
  }
 
  __syncthreads();

  k2 = k_ini;
  l2 = l_ini;

  for (i=0; (k2<=l2) && (i<nW[offset]-1); i++) {

    b = W[offset*MAXLINECOMP + i];

    val1 = ( b      ) & 3;
    siz1 = val1*sizO;

    k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
    l2 = Cshared[val1] + O[siz1 + l2 + 1];

    val2 = ( b >> 2 ) & 3;
    siz2 = val2*sizO;

    k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
    l2 = Cshared[val2] + O[siz2 + l2 + 1];
    
    val3 = ( b >> 4 ) & 3;
    siz3 = val3*sizO;

    k2 = Cshared[val3] + O[siz3 + k2    ] + 1;
    l2 = Cshared[val3] + O[siz3 + l2 + 1];

    val4 = ( b >> 6 ) & 3;
    siz4 = val4*sizO;
    
    k2 = Cshared[val4] + O[siz4 + k2    ] + 1;
    l2 = Cshared[val4] + O[siz4 + l2 + 1];

  }

   __syncthreads();

  if (k2<=l2) {

    // Last block of 4 bases should not be fully filled

    int real_size = nWe[offset];
    int resto = real_size % 4;

    b = W[offset*MAXLINECOMP + nW[offset] - 1];

    switch(resto) {

    case 0:

      val1 = ( b      ) & 3;
      siz1 = val1*sizO;

      k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
      l2 = Cshared[val1] + O[siz1 + l2 + 1];

      val2 = ( b >> 2 ) & 3;
      siz2 = val2*sizO;

      k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
      l2 = Cshared[val2] + O[siz2 + l2 + 1];

      val3 = ( b >> 4 ) & 3;
      siz3 = val3*sizO;

      k2 = Cshared[val3] + O[siz3 + k2    ] + 1;
      l2 = Cshared[val3] + O[siz3 + l2 + 1];

      val4 = ( b >> 6 ) & 3;
      siz4 = val4*sizO;
    
      k2 = Cshared[val4] + O[siz4 + k2    ] + 1;
      l2 = Cshared[val4] + O[siz4 + l2 + 1];
  
      break;

    case 3:

      val1 = ( b      ) & 3;
      siz1 = val1*sizO;

      k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
      l2 = Cshared[val1] + O[siz1 + l2 + 1];

      val2 = ( b >> 2 ) & 3;
      siz2 = val2*sizO;

      k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
      l2 = Cshared[val2] + O[siz2 + l2 + 1];

      val3 = ( b >> 4 ) & 3;
      siz3 = val3*sizO;

      k2 = Cshared[val3] + O[siz3 + k2    ] + 1;
      l2 = Cshared[val3] + O[siz3 + l2 + 1];

      break;

    case 2:

      val1 = ( b      ) & 3;
      siz1 = val1*sizO;

      k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
      l2 = Cshared[val1] + O[siz1 + l2 + 1];

      val2 = ( b >> 2 ) & 3;
      siz2 = val2*sizO;

      k2 = Cshared[val2] + O[siz2 + k2    ] + 1;
      l2 = Cshared[val2] + O[siz2 + l2 + 1];

      break;

    case 1:
      val1 = ( b      ) & 3;
      siz1 = val1*sizO;

      k2 = Cshared[val1] + O[siz1 + k2    ] + 1;
      l2 = Cshared[val1] + O[siz1 + l2 + 1];

      break;
    
    }

  }

  __syncthreads();

  k[offset] = k2;
  l[offset] = l2;

}
*/

void BWExactSearchBackwardGPUWrapper_ex(unsigned int num_bloques, unsigned int tam_bloques, char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, vector *C, vector *C1, comp_matrix *O, size_t num_reads) {
  BWExactSearchBackwardGPU_ex<<<num_bloques,tam_bloques>>>(W, nW, k, l, k_ini, l_ini, C->vector, C1->vector, *O, num_reads);
}

void BWExactSearchForwardGPUWrapper_ex(unsigned int num_bloques, unsigned int tam_bloques, char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, vector *C, vector *C1, comp_matrix *O, size_t num_reads) {
  BWExactSearchForwardGPU_ex<<<num_bloques,tam_bloques>>>(W, nW, k, l, k_ini, l_ini, C->vector, C1->vector, *O, num_reads);
}




void BWExactSearchBackwardGPUWrapper(unsigned int num_bloques, unsigned int tam_bloques, char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, vector *C, vector *C1, comp_matrix *O) {
  BWExactSearchBackwardGPU<<<num_bloques,tam_bloques>>>(W, nW, k, l, k_ini, l_ini, C->vector, C1->vector, *O);
}

void BWExactSearchForwardGPUWrapper(unsigned int num_bloques, unsigned int tam_bloques, char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, vector *C, vector *C1, comp_matrix *O) {
  BWExactSearchForwardGPU<<<num_bloques,tam_bloques>>>(W, nW, k, l, k_ini, l_ini, C->vector, C1->vector, *O);
}

void BWExactFinalResultsBackwardGPUWrapper(unsigned int num_bloques, unsigned int tam_bloques, char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next, unsigned int chunk_size, int *stack_size) {
  BWExactFinalResultsBackwardGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, *O, *rl_prev, *rl_next, chunk_size, stack_size);
}

void BWExactSearchBackwardVectorGPUWrapper(unsigned int num_bloques, unsigned int tam_bloques, char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, vector *C, vector *C1, comp_matrix *O) {
  BWExactSearchBackwardVectorGPU<<<num_bloques,tam_bloques>>>(W, nW, k, l, k_ini, l_ini, C->vector, C1->vector, *O);
}

void BWExactSearchForwardVectorGPUWrapper(unsigned int num_bloques, unsigned int tam_bloques, char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, vector *C, vector *C1, comp_matrix *O) {
  BWExactSearchForwardVectorGPU<<<num_bloques,tam_bloques>>>(W, nW, k, l, k_ini, l_ini, C->vector, C1->vector, *O);
}
