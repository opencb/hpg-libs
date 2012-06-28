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

#ifndef _BW_SEARCH_
#define _BW_SEARCH_

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "BW_io.h"

#define BWiterationVariables() size_t aux_b;

#ifndef VECTOR_O_COMPRESSION

#define BWiteration(k_in, l_in, k_out, l_out, b, C, C1, O)		\
  {									\
    aux_b = (b);							\
    (k_out) = (C1)->vector[aux_b] + (O)->count[aux_b][(k_in)];		\
    (l_out) = (C)->vector[aux_b]  + (O)->count[aux_b][(l_in)+1];	\
  }

#else

#define BWiteration(k_in, l_in, k_out, l_out, b, C, C1, O)		 \
  {									 \
    aux_b = (b);							 \
    (k_out) = (C1)->vector[aux_b] + getOcompValue(aux_b, (k_in)  , (O)); \
    (l_out) = (C)->vector[aux_b]  + getOcompValue(aux_b, (l_in)+1, (O)); \
  }

#endif

void calculateD(vector *D, byte_vector *W, vector *C, vector *C1, comp_matrix *Oi);
void BWRecursiveSearch(char *W, size_t i, unsigned int z, size_t k, size_t l, vector *D, vector *C, vector *C1, comp_matrix *O, results_list *r_list);

inline void BWExactSearchBackward(char *W, int start, int end, vector *C, vector *C1, comp_matrix *O, result *r){

  //printf("CALL Backward\n");
  BWiterationVariables();
  size_t k2, l2;
  int i;

  k2 = r->k; l2 = r->l;

  //printf("B1º -> %lu - %lu\n", k2, l2);

  for(i=end; i>=start; i--) {

    BWiteration(k2, l2, k2, l2, W[i], C, C1, O);
    //printf("B -> %lu -> %lu - %u\n", i, k2, l2);
    if (k2 > l2) break;

  }

  r->k = k2;
  r->l = l2;
  r->start = start;

}

inline void BWExactSearchForward(char *W, int start, int end, vector *C, vector *C1, comp_matrix *Oi, result *r) {
  //printf("CALL Forward\n");

  BWiterationVariables();
  size_t k2, l2;
  int i;

  k2 = r->k;  l2 = r->l;

  //printf("F1º -> %lu - %lu\n", k2, l2);

  for(i=start; i<=end; i++) {

    BWiteration(k2, l2, k2, l2, W[i], C, C1, Oi);
    //printf("F-> %lu -> %lu - %u\n", i, k2, l2);
    if (k2 > l2) break;

  }

  r->k = k2;
  r->l = l2;
  r->end = end;

}

void BWExactSearchBackwardVector(char *W, int start, int end, size_t k, size_t l, size_t **_vec_k, size_t **_vec_l, vector *C, vector *C1, comp_matrix *O);
void BWExactSearchForwardVector(char *W, int start, int end, size_t k, size_t l, size_t **_vec_ki, size_t **_vec_li, vector *C, vector *C1, comp_matrix *Oi);

#ifdef __cplusplus
extern "C" {
#endif
void BWIterativeSearch1(char *W, int start, int end, size_t *vec_k, size_t *vec_l, size_t *vec_ki, size_t *vec_li, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, results_list *r_list);
#ifdef __cplusplus
}
#endif

void BWBackwardSimpleSearch1(char *W, int start, int end, vector *C, vector *C1, comp_matrix *O, result *r, results_list *r_list);

void BWForwardSimpleSearch1(char *W, int start, int end, vector *C, vector *C1, comp_matrix *O, result *r, results_list *r_list);

#endif
