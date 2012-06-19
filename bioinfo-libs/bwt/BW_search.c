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

//TODO: Pasar res en vez de k, l, start y end, en todos los algoritmos.

#include "BW_search.h"

void calculateD(vector *D, byte_vector *W, vector *C, vector *C1, comp_matrix *Oi) {

  BWiterationVariables();
  size_t i;
  size_t k,l;
  unsigned int z;
  size_t nBlast = Oi->m_count-2; //Substract 2, vector = starts in -1

  D->n = W->n;

  k=0;
  l=nBlast;
  z=0;

  D->vector = (unsigned int *) malloc(D->n * sizeof(unsigned int));
  checkMalloc(D, "calculateD");

  for (i=0; i<W->n; i++) {

    BWiteration(k, l, k, l, W->vector[i], C, C1, Oi);

    if (k>l) {
      k=0;
      l=nBlast;
      z++;
    }

    D->vector[i] = z;

  }

}

void BWRecursiveSearch(char *W, size_t i, unsigned int z, size_t k, size_t l, vector *D, vector *C, vector *C1, comp_matrix *O, results_list *r_list) {

  BWiterationVariables();
  size_t b;
  size_t k2, l2;

  result r;

  if (z<1) { /*(z-1<0)*/

    if (i<1) { /*(i-1<0)*/

      for (b=0;b<nA;b++) {

	BWiteration(k, l, k2, l2, b, C, C1, O);

	if (k2<=l2) {

	  if ( ((size_t) W[i]) == b) {
	    new_result(&r, k2, l2, 0, 0, 0);
	    add_result(&r, r_list);
	  }

	}

      }

    } else { /*(i-1>=0)*/

      if (z < D->vector[i-1]) return;

      for (b=0; b<nA; b++) {

	BWiteration(k, l, k2, l2, b, C, C1, O);
	
	if (k2<=l2) {
	  
	  if ( ((size_t) W[i]) == b) {

	    BWRecursiveSearch(W, i-1, z, k2, l2, D, C, C1, O, r_list); //Match

	  }

	}
	
      }

    }

  } else { /* (z-1>=0) */

    if (i<1) { /*(i-1<0)*/

      if (z-1<D->vector[i]) {

	new_result(&r, k, l, 0, 0, 0);
	add_result(&r, r_list);

	for (b=0;b<nA;b++) {

	  BWiteration(k,l,k2,l2,b, C, C1, O);
	  
	  if (k2<=l2) {
	    if (((size_t)W[i]) == b) {
	      new_result(&r, k2, l2, 0, 0, 0);
	      add_result(&r, r_list);
	    } else {
	      new_result(&r, k2, l2, 0, 0, 0);
	      add_result(&r, r_list);
	    }

	  }

	}

      } else {

	new_result(&r, k, l, 0, 0, 0);
	add_result(&r, r_list);

	for (b=0;b<nA;b++) {

	  BWiteration(k,l,k2,l2,b,C,C1,O);

	  if (k2<=l2) {
	    
	    BWRecursiveSearch(W, i, z-1, k2, l2, D, C, C1, O, r_list); // Insertion

	    if (((size_t)W[i]) == b) {
	      new_result(&r, k2, l2, 0, 0, 0);
	      add_result(&r, r_list);
	    } else {
	      new_result(&r, k2, l2, 0, 0, 0);
	      add_result(&r, r_list);
	    }

	  }

	}

      }

    } else { /*(i-1>=0)*/

      if (z-1>=D->vector[i-1]) BWRecursiveSearch(W, i-1, z-1, k, l, D, C, C1, O, r_list); // Deletion

      for (b=0;b<nA;b++) {

	BWiteration(k,l,k2,l2,b,C,C1,O);

	if (k2<=l2) {

	  if (z-1>=D->vector[i]) BWRecursiveSearch(W, i, z-1, k2, l2, D, C, C1, O, r_list); //Insertion

	  if (((size_t)W[i]) == b) {

	    if (z>=D->vector[i-1]) BWRecursiveSearch(W, i-1, z, k2, l2, D, C, C1, O, r_list); //Match

	  } else {

	    if (z-1>=D->vector[i-1]) BWRecursiveSearch(W, i-1, z-1, k2, l2, D, C, C1, O, r_list); //Missmatch

	  }

	}

      }

    }

  }

}

void BWExactSearchBackwardVector(char *W, int start, int end, size_t k, size_t l, size_t **_vec_k, size_t **_vec_l, vector *C, vector *C1, comp_matrix *O) {

  if (k > l)       return;
  if (start > end) return;

  BWiterationVariables();
  size_t k2, l2;
  size_t *vec_k, *vec_l;
  int last;
  int i, j;

  last = end-start;
  vec_k = (size_t *) malloc((last+1) * sizeof(size_t));
  vec_l = (size_t *) malloc((last+1) * sizeof(size_t));

  k2 = k;
  l2 = l;

  for(i=end, j=last; i>=start; i--, j--) {

    BWiteration(k2, l2, k2, l2, (size_t)W[i], C, C1, O);

    vec_k[j] = k2;
    vec_l[j] = l2;

    if (k2 > l2) {
      i--; j--;
      break;
    }

  }

  for(;i>=start; i--, j--) {
    vec_k[j] = k2;
    vec_l[j] = l2;
  }

  *_vec_k = vec_k;
  *_vec_l = vec_l;

}

void BWExactSearchForwardVector(char *W, int start, int end, size_t k, size_t l, size_t **_vec_ki, size_t **_vec_li, vector *C, vector *C1, comp_matrix *Oi) {

  if (k > l) return;
  if (start > end) return;

  BWiterationVariables();
  size_t k2, l2;
  size_t *vec_ki, *vec_li;
  int last;
  int i, j;

  last = end-start;
  vec_ki = (size_t *) malloc((last+1) * sizeof(size_t));
  vec_li = (size_t *) malloc((last+1) * sizeof(size_t));

  k2 = k;
  l2 = l;

  for(i=start, j=0; i<=end; i++, j++) {

    BWiteration(k2, l2, k2, l2, W[i], C, C1, Oi);

    vec_ki[j] = k2;
    vec_li[j] = l2;

    if (k2 > l2) {
      i++; j++;
      break;
    }

  }

  for(; i<=end; i++, j++) {
    vec_ki[j] = k2;
    vec_li[j] = l2;
  }

  *_vec_ki = vec_ki;
  *_vec_li = vec_li;

}

void BWIterativeSearch1(char *W, int start, int end, size_t *vec_k, size_t *vec_l, size_t *vec_ki, size_t *vec_li, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, results_list *r_list) {

  BWiterationVariables();

  size_t _k, _l, _ki, _li, _k_aux, _l_aux, _ki_aux, _li_aux;
  size_t results, results_last;

  int i, j, half, n;

  result r;

  n = end - start;
  half = n / 2;

  if (vec_k[0] <= vec_l[0]) {
    new_result(&r, vec_k[0], vec_l[0], start, end, 0);
    add_result(&r, r_list);
  }

  results = vec_l[0] - vec_k[0];

  results_last = results;
  _k  = vec_k[1];
  _l  = vec_l[1];
  results = _l  - _k;

  //printf("B -> %d: %d -> %d, %u, %u\n", 0, results, results_last, _k, _l);

  if (results != results_last) {

    //printf("*B -> %d: %d -> %d, %u, %u\n", 0, results, results_last, _k, _l);

    //printf("%d: %d -> %d\n", 0, results, results_last);
    //Deletion
    new_result(&r, _k, _l, start, end, 0);
    add_mismatch(&r, 1, 4, start);
    add_result(&r, r_list);

    for (size_t b=0;b<nA;b++) {

      BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);

      //printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

      if (_k_aux > _l_aux) continue;
      //printf("*W -> %d, %d - %d\n", b, _k_aux, _l_aux);

      size_t b_w = (size_t) W[start];

      //Missmatch
      if (b!=b_w) {
	new_result(&r, _k_aux, _l_aux, start, end, 0);
	add_mismatch(&r, 2, b, start);
	add_result(&r, r_list);
      }

      //Insertion
      BWiteration(_k_aux, _l_aux, _k_aux, _l_aux, b_w, C, C1, O);

      if (_k_aux <= _l_aux) {
	new_result(&r, _k_aux, _l_aux, start, end, 0);
	add_mismatch(&r, 3, b, start);
	add_result(&r, r_list);
      }

    }

  }

  for (i=start+2, j=2; j<=half; i++, j++) {

    results_last = results;
    _k = vec_k[j];
    _l = vec_l[j];
    results = _l  - _k;

    //printf("B -> %d: %d -> %d, %u, %u\n", j-1, results, results_last, _k, _l);

    if (results == results_last) continue;

    //printf("*B -> %d: %d -> %d, %u, %u\n", j-1, results, results_last, _k, _l);

    //Deletion
    new_result(&r, _k, _l, i-1, end, 0);
    add_mismatch(&r, 1, 4, i-1);
    BWExactSearchBackward(W, start, i-2, C, C1, O, &r);
    if (r.k<=r.l) add_result(&r, r_list);

    for (size_t b=0;b<nA;b++) {

      BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);

      if (_k_aux > _l_aux) continue;

      //Mismatch
      if (b!=(size_t)W[i-1]) {
	new_result(&r, _k_aux, _l_aux, i-1, end, 0);
	add_mismatch(&r, 2, b, i-1);
	BWExactSearchBackward(W, start, i-2, C, C1, O, &r);
	if (r.k<=r.l) add_result(&r, r_list);
      }

      //Insertion
      new_result(&r, _k_aux, _l_aux, i-1, end, 0);
      add_mismatch(&r, 3, b, i-1);
      BWExactSearchBackward(W, start, i-1, C, C1, O, &r);
      if (r.k<=r.l) add_result(&r, r_list);

    }

  }

  //printf("\n");

  //TODO: Gestionar bien los errores de la busqueda al revés con Si y restas (precalcular Si con |X| - pos)
  half--;
  results = vec_li[n] - vec_ki[n];

  results_last = results;
  _ki  = vec_ki[n-1];
  _li  = vec_li[n-1];
  results = _li - _ki;

  //printf("F-> %d: %d -> %d, %u, %u\n", n, results, results_last, _ki, _li);

  if (results != results_last) {

    //printf("*F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

    //Deletion
    new_result(&r, _ki, _li, start, end, 1);
    add_mismatch(&r, 1, 4, end);
    add_result(&r, r_list);

    for (size_t b=0;b<nA;b++) {

      BWiteration(_ki, _li, _ki_aux, _li_aux, b, C, C1, Oi);

      if (_ki_aux > _li_aux) continue;

      size_t b_w = (size_t) W[end];

      //Mismatch
      if (b!=b_w) {
	new_result(&r, _ki_aux, _li_aux, start, end, 1);
	add_mismatch(&r, 2, b, end);
	add_result(&r, r_list);
      }

      //Insertion
      BWiteration(_ki_aux, _li_aux, _ki_aux, _li_aux, b_w, C, C1, Oi);

      //printf("\tI -> %d - %d\n", _ki_aux, _li_aux);

      if (_ki_aux <= _li_aux){
	new_result(&r, _ki_aux, _li_aux, start, end, 1);
	add_mismatch(&r, 3, b, end);
	add_result(&r, r_list);
      }

    }

  }

  for(i=end-2,j=n-2; j>=half; i--, j--) {

    results_last = results;
    _ki  = vec_ki[j];
    _li  = vec_li[j];
    results = _li - _ki;

    //printf("F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

    if (results == results_last) continue;

    //printf("*F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

    //TODO: Anadir contador para podar cuando se que ya he encontrado las cadenas que permite la variabilidad en este simbolo.

    //Deletion
    new_result(&r, _ki, _li, start, i+1, 1);
    add_mismatch(&r, 1, 4, i+1);
    BWExactSearchForward(W, i+2, end, C, C1, Oi, &r);
    if (r.k<=r.l) add_result(&r, r_list);

    for (size_t b=0;b<nA;b++) {

      BWiteration(_ki, _li, _ki_aux, _li_aux, b, C, C1, Oi);

      //printf("W -> %d, %d - %d\n", b, _ki_aux, _li_aux);

      if (_ki_aux > _li_aux) continue;

      //Mismatch
      if (b!= (size_t) W[i+1]) {
	new_result(&r, _ki_aux, _li_aux, start, i+1, 1);
	add_mismatch(&r, 2, b, i+1);
	BWExactSearchForward(W, i+2, end, C, C1, Oi, &r);
	if (r.k<=r.l) add_result(&r, r_list);
      }

      //Insertion
      new_result(&r, _ki_aux, _li_aux, start, i+1, 1);
      add_mismatch(&r, 3, b, i+1);
      BWExactSearchForward(W, i+1, end, C, C1, Oi, &r);
      if (r.k<=r.l) add_result(&r, r_list);

    }

  }

  //printf("\n");

}

/*
void BWIterativeSearch1bis(char *W, int start, int end, size_t *vec_k, size_t *vec_l, size_t *vec_ki, size_t *vec_li, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, results_list *r_list) {

  n = end - start;
  half = n / 2;


  BWBackwardSimpleSearch1(W, int start, int end, C, C1, O, result *res, r_list);

  half--;

  BWForwardSimpleSearch1(W, int start, int end, C, C1, Oi, result *res, r_list);

  //TODO: Borrar último resultado si es la busqueda exacta.

}
*/

void BWBackwardSimpleSearch1(char *W, int start, int end, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list) {

  BWiterationVariables();
  size_t _k,_l, _k_next, _l_next, _k_aux, _l_aux;
  size_t results, results_next;
  int i;

  result r;

  _k_next = res->k;
  _l_next = res->l;
  results_next = _l_next - _k_next;

  for(i=end; i>=start; i--) {

    _k = _k_next;
    _l = _l_next;

    //printf("%d:\n", i);

    if (_k > _l) return;

    BWiteration(_k, _l, _k_next, _l_next, W[i], C, C1, O);
    results      = results_next;
    results_next = _l_next - _k_next;
    //printf("(%lu, %lu, %lu)\t", results, _k, _l);

    if (results == results_next) continue;

    //Deletion
    new_result(&r, _k, _l, i, end, 0);
    add_mismatch(&r, 1, 4, i);
    BWExactSearchBackward(W, start, i-1, C, C1, O, &r);
    if (r.k<=r.l) add_result(&r, r_list);

    for (size_t b=0;b<nA;b++) {

      BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);

      //printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

      if (_k_aux > _l_aux) continue;

      //Mismatch
      if (b!=(size_t)W[i]) {
	new_result(&r, _k_aux, _l_aux, i, end, 0);
	add_mismatch(&r, 2, b, i);
	BWExactSearchBackward(W, start, i-1, C, C1, O, &r);
	if (r.k<=r.l) add_result(&r, r_list);
      }

      //Insertion
      new_result(&r, _k_aux, _l_aux, i, end, 0);
      add_mismatch(&r, 3, b, i);
      BWExactSearchBackward(W, start, i, C, C1, O, &r);
      if (r.k<=r.l) add_result(&r, r_list);

    }

  }
  //printf("\n");
  if (_k_next > _l_next) return;

  new_result(&r, _k_next, _l_next, start, end, 0);
  add_mismatch(&r, 0, 4, start);
  add_result(&r, r_list);

}

void BWForwardSimpleSearch1(char *W, int start, int end, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list) {

  BWiterationVariables();
  size_t _k, _l, _k_next, _l_next, _k_aux, _l_aux;
  size_t results, results_next;
  int i;

  result r;

  _k_next = res->k;
  _l_next = res->l;
  results_next = _l_next - _k_next;

  for(i=start; i<=end; i++) {

    _k = _k_next;
    _l = _l_next;

    //printf("%d:\n", i);

    if (_k > _l) return;

    BWiteration(_k, _l, _k_next, _l_next, W[i], C, C1, O);
    results      = results_next;
    results_next = _l_next - _k_next;
    if (results == results_next) continue;

    //Deletion
    new_result(&r, _k, _l, start, i, 1);
    add_mismatch(&r, 1, 4, i);
    BWExactSearchForward(W, i+1, end, C, C1, O, &r);
    if (r.k<=r.l) add_result(&r, r_list);

    for (size_t b=0;b<nA;b++) {

      BWiteration(_k, _l, _k_aux, _l_aux, b, C, C1, O);

      //printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

      if (_k_aux > _l_aux) continue;

      //Mismatch
      if (b!=(size_t)W[i]) {
	new_result(&r, _k_aux, _l_aux, start, i, 1);
	add_mismatch(&r, 2, b, i);
	BWExactSearchForward(W, i+1, end, C, C1, O, &r);
	if (r.k<=r.l) add_result(&r, r_list);
      }

      //Insertion
      new_result(&r, _k_aux, _l_aux, start, i, 1);
      add_mismatch(&r, 3, b, i);
      BWExactSearchForward(W, i, end, C, C1, O, &r);
      if (r.k<=r.l) add_result(&r, r_list);

    }

  }

  if (_k_next > _l_next) return;

  new_result(&r, _k_next, _l_next, start, end, 1);
  add_mismatch(&r, 0, 4, end);
  add_result(&r, r_list);

}
