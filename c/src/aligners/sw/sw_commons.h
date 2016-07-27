#ifndef SW_COMMONS_H
#define SW_COMMONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>
#include <omp.h>

#include "macros.h"


#ifdef TIMING
extern double *sse_matrix_t, *sse_tracking_t;
#endif // TIMING

#ifdef __AVX512__
#define SIMD_DEPTH 16
#define SIMD_ALIGNED 64
#else
#ifdef __AVX2__
#define SIMD_DEPTH 8
#define SIMD_ALIGNED 32
#else
#define SIMD_DEPTH 4
#define SIMD_ALIGNED 16
#endif
#endif


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

typedef float subst_matrix_t[128][128];

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

/**
 * @brief Structure for the Smith-Waterman algorithm parameters.
 *
 * Structure for the Smith-Waterman algorithm parameters:
 * gap penalties (gap open and extending) and the substitution
 * score matrix.
 */
typedef struct sw_optarg {
  float match;   /**< Score for match. */
  float mismatch; /**< Penalty for mismatch. */
  float gap_open;   /**< Penalty for gap openning. */
  float gap_extend; /**< Penalty for gap extending. */
  char *subst_matrix_name; /**< Substitution score matrix name. */
  subst_matrix_t subst_matrix; /**< Substitution score matrix. */
} sw_optarg_t;

//------------------------------------------------------------------------------------

/**
 * @brief Constructor for the @a sw_optarg_t structure.
 * @param gap_open Penalty for gap openning.
 * @param gap_extend Penalty for gap extending.
 * @param subst_matrix_name Substitution score matrix name.
 * @return Pointer to the new structure.
 *
 * Constructor function that allocates memory for
 * the Smith-Waterman algorithm parameters.
 */
sw_optarg_t* sw_optarg_new(float match, float mismatch,
                           float gap_open, float gap_extend, char *subst_matrix_name);

//------------------------------------------------------------------------------------

/**
 * @brief Destructor for the @a sw_optarg_t structure.
 * @param optarg_p[out] pointer to the structure to free
 *
 * @a sw_optarg_t destructor that frees the memory previously
 * allocated by the constructor @a sw_optarg_new.
 */
void sw_optarg_free(sw_optarg_t* optarg_p);

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

/**
 * @brief Output structure for Smith-Waterman algorithm.
 *
 * Output structure for the Smith-Waterman algorithm,
 * basically, it contains the pointers to the aligned sequences,
 * start positions and scores.
 */
typedef struct sw_multi_output {
  unsigned int num_queries;  /**< Num. queries. */
  char **query_map_p; /**< Pointers to the aligned target sequences (queries). */
  char **ref_map_p; /**< Pointers to the aligned reference sequences. */
  unsigned int *query_map_len_p; /**< Pointers to the aligned target sequences lengths. */
  unsigned int *ref_map_len_p; /**< Pointers to the aligned reference sequences lenghts. */
  unsigned int *query_start_p; /**< Pointers to the start positions in the target sequences (queries). */
  unsigned int *ref_start_p; /**< Pointers to the start positions in the reference sequences. */
  float *score_p;      /**< Pointers to the resulting scores. */
} sw_multi_output_t;

//------------------------------------------------------------------------------------

/**
 * @brief Constructor for the @a sw_multi_output_t structure.
 * @param num_queries Number of queries (target sequences)
 * @return Pointer to the new structure.
 *
 * @a sw_simd_output_t constructor that allocates memory for
 * the aligned sequences and lengths pointers, the number of pointers
 * depends on the depth.
 */
sw_multi_output_t* sw_multi_output_new(unsigned int num_queries);

//------------------------------------------------------------------------------------

/**
 * @brief Destructor for the @a sw_simd_output_t structure.
 * @param output_p[out] pointer to the structure to free
 *
 * @a sw_simd_output_t destructor that frees the memory previously
 * allocated by the constructor @a sw_simd_output_new.
 */
void sw_multi_output_free(sw_multi_output_t* output_p);

//------------------------------------------------------------------------------------

/**
 * @brief Saves a @a sw_simd_output_t structure.
 * @param output_p pointer to the structure storing the aligned sequences
 * @param num_alignments number of alignments to save
 * @param file_p pointer to the file where to store the alignments
 *
 * Saves num_alignments aligned sequences, in pairs: target/reference
 * into the given file
 */
void sw_multi_output_save(int num_alignments, sw_multi_output_t* output_p, FILE *file_p);
int sw_multi_output_string(int num_alignments, sw_multi_output_t* output_p, char *buf);

//------------------------------------------------------------------------------------

#endif // SW_COMMONS

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
