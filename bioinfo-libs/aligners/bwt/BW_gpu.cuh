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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _BW_GPU_
#define _BW_GPU_

#define BWiterationVariablesGPU() size_t aux_b;

#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION

#define BWiterationGPU(k_in, l_in, k_out, l_out, b, C, C1, O)		\
  {									\
    aux_b = (b);							\
    (k_out) = (C1)[aux_b] + getOcompValueGPU(aux_b, (k_in)  , (O));	\
    (l_out) = (C)[aux_b]  + getOcompValueGPU(aux_b, (l_in)+1, (O));     \
  }

#else

#define BWiterationGPU(k_in, l_in, k_out, l_out, b, C, C1, O)		\
  {									\
    aux_b = (b);			       				\
    (k_out) = (C1)[aux_b] + (O).desp[aux_b][(k_in)  ];	                \
    (l_out) = (C)[aux_b]  + (O).desp[aux_b][(l_in)+1];	                \
  }

#endif

#define manageCudaError()						\
  error = cudaGetLastError();						\
  if (error != cudaSuccess) {						\
    fprintf(stderr, "Error kernel: %s\n", cudaGetErrorString(error));	\
    exit(1);								\
  }

void copyVectorGPU(vector *device, vector *host, size_t data_size);

void readCompMatrixGPU(comp_matrix *matrix, const char *directory, const char *name);
void copyCompMatrixGPU(comp_matrix *device, comp_matrix *host);
void freeCompMatrixGPUHost(comp_matrix *matrix);
void freeCompMatrixGPUDevice(comp_matrix *matrix);

void BWExactSearchBackwardGPUWrapper(unsigned int num_bloques, unsigned int tam_bloque, char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, vector *C, vector *C1, comp_matrix *O);

void BWExactSearchForwardGPUWrapper(unsigned int num_bloques, unsigned int tam_bloque, char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, vector *C, vector *C1, comp_matrix *O);

//-----------------------------------------------------------------------------
// _ex functions, for extended fucntions, 
//
// These functions adds a new parameter num_reads, in addition, the nW vector
// contains indices to the start of each read in W vector, instead of the length
// of the read as the original functions ("no _ex") do
//
//-----------------------------------------------------------------------------
void BWExactSearchBackwardGPUWrapper_ex(unsigned int num_bloques, unsigned int tam_bloques, char *W, 
     				        unsigned int *nW, size_t *k, size_t *l, size_t k_ini, 
					size_t l_ini, vector *C, vector *C1, comp_matrix *O, size_t num_reads);

void BWExactSearchForwardGPUWrapper_ex(unsigned int num_bloques, unsigned int tam_bloques, char *W, 
     				       unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, 
				       vector *C, vector *C1, comp_matrix *O, size_t num_reads);

void BWExactSearchBackwardGPUSeedsWrapper_ex(unsigned int num_bloques, unsigned int tam_bloques, 
     					     char *W, unsigned int *nW, size_t *k, size_t *l, 
					     size_t k_ini, size_t l_ini, vector *C, vector *C1, 
					     comp_matrix *O, size_t num_reads, size_t seed_size, 
					     size_t min_seed_size, size_t num_max_seeds);

void BWExactSearchForwardGPUSeedsWrapper_ex(unsigned int num_bloques, unsigned int tam_bloques, 
     					    char *W, unsigned int *nW, size_t *k, size_t *l, 
					    size_t k_ini, size_t l_ini, vector *C, vector *C1, 
					    comp_matrix *O, size_t num_reads, size_t seed_size, 
					    size_t min_seed_size, size_t num_max_seeds);

//-----------------------------------------------------------------------------
// end of _ex functions
//-----------------------------------------------------------------------------

void BWExactFinalResultsBackwardGPUWrapper(unsigned int num_bloques, unsigned int tam_bloques, char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next, unsigned int chunk_size, int *stack_size);

void BWExactSearchBackwardVectorGPUWrapper(unsigned int num_bloques, unsigned int tam_bloque, char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, vector *C, vector *C1, comp_matrix *O);
void BWExactSearchForwardVectorGPUWrapper(unsigned int num_bloques, unsigned int tam_bloque, char *W, unsigned int *nW, size_t *k, size_t *l, size_t k_ini, size_t l_ini, vector *C, vector *C1, comp_matrix *O);

/*
__global__ void BWExactIterativeSearchGPURev(char *W, int *nW, int *nWe, int *k, int *l, int k_ini, int l_ini, int *C,  int *O,  int sizO);
__global__ void BWExactIterativeSearchGPURev(char *W, int *nW, int *nWe, int *k, int *l, int k_ini, int l_ini, int *C,  int *O,  int sizO);
*/

#endif
