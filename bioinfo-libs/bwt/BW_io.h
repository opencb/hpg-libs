/*
    bwa_gpu a set of tools which allow short sequence alignment using the Burrows-Wheeler transform    usign both CPU and GPU approaches.
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

#ifndef BW_IO
#define BW_IO

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include<stddef.h>
#include<malloc.h>

#define nA              4

#define MAX_MISMATCHES  1

#define INDEX_EXOME 100
#define IDMAX 3

#define MAXLINE     200
#define MAXLINECOMP  50
#define MAX_READ_GPU 256000

#define TRUE            1
#define FALSE           0

#define DELETION  1
#define MISMATCH  2
#define INSERTION 3

#define DD  -1
#define AA  0
#define CC  1
#define GG  2
#define TT  3

extern const char alph_rep[];
//const char alph_rep[] ={'A','C','G','T'};

#ifdef VERBOSE_DBG
#define printUIntMatrix(M,n,m)			          \
  printf("Matrix " #M ":\n");			          \
  for (size_t i=0; i<((size_t) (n)); i++) {		  \
    printf("%lu: ", i);					  \
    for (size_t j=0; j<((size_t) (m)); j++) {		  \
      printf("%u ", (M)[i][j]);     		          \
    }                                                     \
    printf("\n");				          \
  }
#else
#define printUIntMatrix(M,n,m);
#endif

#ifdef VERBOSE_DBG
#define print32bitUIntMatrix(M,n,siz,m)					\
  {									\
    unsigned int bit, b32;						\
    printf("Matrix " #M ":\n");						\
    for (size_t i=0; i<((size_t) (n)); i++) {				\
      printf("%lu: ", i);						\
      for (size_t j=0; j<((size_t) (m)); j++) {				\
	b32 = j / 32;							\
	bit  = j % 32;							\
	printf("%u ", ((M)[i][b32] >> (32-(bit+1))) & 1);		\
	if (bit+1 == 32) printf("\t");					\
      }									\
      printf("\n");							\
    }									\
  }
#else
#define print32bitUIntMatrix(M,n,siz,m);
#endif

#ifdef VERBOSE_DBG

#ifndef VECTOR_O_COMPRESSION

#define printCompMatrix(O)				\
  printUIntMatrix((O).count, (O).n_count, (O).m_count);

#else

#define printCompMatrix(O)						\
  print32bitUIntMatrix((O).count, (O).n_count, (O).siz_count, (O).m_count); \
  printUIntMatrix((O).desp, (O).n_desp, (O).m_desp);
  
#endif

#else
#define printCompMatrix(O);
#endif

#ifdef VERBOSE_DBG
#define printUIntVector(V,n)				\
  printf("Vector " #V ":\n");				\
  for (size_t i=0; i<((size_t) (n)); i++) {		\
    if(((int)(V)[i])==-1)				\
      printf("- ");				        \
    else                                                \
      printf("%u ", (V)[i]);				\
  }                                                     \
  printf("\n");
#else
#define printUIntVector(V,n);
#endif

#ifdef VERBOSE_DBG
#define printString(S) printf("String " #S ":\n%s\n", S);			
#else
#define printString(S);
#endif

#define checkMalloc(D, path)						 \
  if ((D)==NULL) {							 \
    fprintf(stderr, "Data structure " #D " in %s is too large", (path)); \
    exit(1);								 \
  }

#define checkFileOpen(fp, path)						\
  if (!(fp)) {								\
    fprintf(stderr, "Error opening file: %s\n", (path));		\
    exit(1);								\
  }

#define checkFileRead(err, nmemb, path)					\
  if ((err) != (nmemb)) {						\
    fprintf(stderr, "Error reading file '%s'\n", (path));		\
    exit(1);								\
  }

#define checkFileWrite(err, nmemb, path)			\
  if ((err) != (nmemb)) {					\
    fprintf(stderr, "Error writing file '%s'\n", (path));	\
    exit(1);							\
  }

#define timevars() struct timeval t1, t2;

#define tic(msg)				   \
  printf(">> " msg "\n");			   \
  fflush(stdout);				   \
  gettimeofday(&t1, NULL);

#define toc()								                 \
  gettimeofday(&t2, NULL);						                 \
  printf("<< Finished %.0f usecs\n", (t2.tv_sec-t1.tv_sec)*1e6+(t2.tv_usec-t1.tv_usec)); \
  fflush(stdout);

typedef struct {

  unsigned int *count[nA];
  size_t n_count;
  size_t m_count;

#ifdef VECTOR_O_COMPRESSION

  size_t siz_count;
  unsigned int *desp[nA];
  size_t n_desp;
  size_t m_desp;

#endif

} comp_matrix;

typedef struct {

  unsigned int *vector;
  size_t n;

} vector;

typedef struct {

  unsigned int *vector;
  size_t siz;
  size_t n;
  size_t ratio;

} comp_vector;

typedef struct {

  char *vector; //TODO: Change to signed char just in case.
  size_t n;

} byte_vector;

typedef struct {
  char chromosome[INDEX_EXOME*IDMAX];
  unsigned int start[INDEX_EXOME];
  unsigned int end[INDEX_EXOME];
  unsigned int offset[INDEX_EXOME];
  unsigned int size;
} exome;

typedef struct {
  size_t k, l;
  int start, end;
  char dir;
  char num_mismatches;
  char err_kind[MAX_MISMATCHES];
  char base[MAX_MISMATCHES];
  int position[MAX_MISMATCHES];
} result;

typedef struct {
  result *list;
  size_t n;
  size_t max_results;
  size_t read_index;
} results_list;

inline void new_results_list(results_list *r_list, size_t max_results) {

  r_list->list = (result *) malloc(max_results * sizeof(result));
  checkMalloc(r_list->list, "new_result_list");

  r_list->n=0;
  r_list->max_results = max_results;

}


inline void new_result(result *r, size_t k, size_t l, int start, int end, char dir) {
    r->k = k;
    r->l = l;
    r->start = start;
    r->end = end;
    r->dir = dir;
    r->num_mismatches = 0;
}

inline void add_mismatch(result *r, char err_kind, char base, int position) {

  if (r->num_mismatches < MAX_MISMATCHES) {

    size_t pos_mismatch = r->num_mismatches;

    r->err_kind[pos_mismatch] = err_kind;
    r->base[pos_mismatch] = base;
    r->position[pos_mismatch] = position;

    r->num_mismatches++;

  } else {

    fprintf(stderr, "Number of allowed mismatches exceeded: %d\n", r->num_mismatches);
    exit(1);

  }

}

inline void copy_result(result *dest, result *orig) {

  new_result(
	     dest,
	     orig->k,
	     orig->l,
	     orig->start,
	     orig->end,
	     orig->dir
	     );

  for(int i=0; i<orig->num_mismatches; i++) {
    add_mismatch(dest, orig->err_kind[i], orig->base[i], orig->position[i]);
  }

}

inline void add_result(result *orig, results_list *r_list) {

  if (r_list->n < r_list->max_results) {

    result *dest;
    dest = &r_list->list[r_list->n];
    r_list->n++;

    copy_result(dest, orig);

  } else {

    fprintf(stderr, "Number of allowed results exceeded: %lu", r_list->max_results);
    exit(1);

  }

}

#ifdef VECTOR_O_COMPRESSION

inline unsigned int bitcount(unsigned int uint32bit) {//TODO: Probar la llamada en C para ver si usa la instrucción máquina. //TODO: Programar compresión de 64 bits tb.

  //Parallel binary bit add
  uint32bit = uint32bit - ((uint32bit >> 1) & 0x55555555);
  uint32bit = (uint32bit & 0x33333333) + ((uint32bit >> 2) & 0x33333333);
  return (((uint32bit + (uint32bit >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;

}

inline unsigned int getOcompValue(size_t n, size_t m, comp_matrix *O) {

  size_t pos, desp;

  pos  = m / 32;
  desp = (32 - ((m+1) % 32)) % 32; //TODO cambiar compresión para no tener que restar 32 cada vez (almacenar los bits al revés).

  return O->desp[n][pos] + bitcount(O->count[n][pos] >> desp);

}

#endif

inline char getBfromO(size_t m, comp_matrix *O) {

  size_t n = O->n_count;

#ifndef VECTOR_O_COMPRESSION

  for(size_t i=0; i<n; i++)
    if ( (O->count[i][m] < O->count[i][m+1]) ) return i;

#else

  m++;

  size_t pos, desp;
  pos  = m / 32;
  desp = (32 - ((m+1) % 32)) % 32;

  for(size_t i=0; i<n; i++) {
    if ( ( (O->count[i][pos] >> desp) & 1) != 0) return i;
  }

#endif

  return -1;

}

inline unsigned int getScompValue(size_t m, comp_vector *Scomp, vector *C, comp_matrix *O) {

  size_t j;
  size_t i;

  i=m; j=0;

  char b_aux;

  while (i % Scomp->ratio) {

    b_aux = getBfromO(i, O);

    if ((int) b_aux == -1) {

      i=0;

    } else {

#ifndef VECTOR_O_COMPRESSION
      i = C->vector[(int)b_aux] + O->count[(int)b_aux][i+1/*0 is -1*/];
#else
      i = C->vector[(int)b_aux] + getOcompValue((int)b_aux, i+1, O);
#endif

    }

    j++;
  }

  return (Scomp->vector[i / Scomp->ratio] + j) % (O->m_count-1);

}

inline unsigned int getScompValueB(size_t m, comp_vector *Scomp, vector *C, comp_matrix *O, byte_vector *B) {

  size_t j;
  size_t i;

  i=m; j=0;

  char b_aux;

  while (i % Scomp->ratio) {

    b_aux = B->vector[i];

    if ((int) b_aux == -1) {

      i=0;

    } else {

#ifndef VECTOR_O_COMPRESSION
      i = C->vector[(int)b_aux] + O->count[(int)b_aux][i+1/*0 is -1*/];
#else
      i = C->vector[(int)b_aux] + getOcompValue((int)b_aux, i+1, O);
#endif

    }

    j++;
  }

  return (Scomp->vector[i / Scomp->ratio] + j) % (O->m_count-1);

}

#ifdef __cplusplus
extern "C" {
#endif

void reverseStrandC(vector *r_C, vector *s_C, vector *r_C1, vector *s_C1);
void reverseStrandO(comp_matrix *r_O, comp_matrix *s_O);
void freeCompMatrix(comp_matrix *matrix);

void readUIntVector(vector *vector, const char *directory, const char *name);
void readUIntCompVector(comp_vector *vector, const char *directory, const char *name);
void readCharVector(byte_vector *vector, const char *directory, const char *name);
void readCompMatrix(comp_matrix *matrix, const char *directory, const char *name);

void saveUIntVector(vector *vector, const char *directory, const char *name);
void saveUIntCompVector(comp_vector *vector, const char *directory, const char *name);
void saveCharVector(byte_vector *vector, const char *directory, const char *name);
void saveCompMatrix(comp_matrix *matrix, const char *directory, const char *name);

#ifdef __cplusplus
}
#endif

void initReplaceTable();

char *replaceBases(char *uncoded, char *coded, size_t length);

int nextFASTAToken(FILE *queries_file, char *original, char *contents, unsigned int *nquery, char *compressed, unsigned int *ncompress);

void load_duplicated_reference(byte_vector *Xorig, byte_vector *X, const char *path);
//unsigned load_duplicated_references(char **_X, unsigned int *max_nX, const char *path);

void revstring(char *X, size_t nX);

size_t comp4basesInChar(char *X, size_t nX, char *Y);
unsigned int binsearch(unsigned int *array, unsigned int size, size_t key);

void load_exome_file(exome *ex, const char *name);
void initialize_init_mask();
int write_results(results_list *r_list, exome* ex, comp_vector *S, comp_vector *Si, vector *C, comp_matrix *O, comp_matrix *Oi, char *mapping, int nW, int type, FILE *fp);
void free_results_list(results_list *r_list);




#endif // BW_IO
