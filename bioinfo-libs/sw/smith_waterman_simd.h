#ifndef SMITH_WATERMAN_SIMD_H
#define SMITH_WATERMAN_SIMD_H

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <immintrin.h>

#ifdef __AVX__
#define SIMD_DEPTH 8
#define SIMD_ALIGN 32
#else
#define SIMD_DEPTH 4
#define SIMD_ALIGN 16
#endif // __AVX__


typedef float subst_matrix_t[128][128];

//====================================================================================
// Smith-Waterman structures and functions (SIMD version)
//====================================================================================

typedef struct sw_simd_context {
    float gap_open;   /**< Penalty for gap openning. */
    float gap_extend; /**< Penalty for gap extending. */

#ifdef __AVX__
     __m256 zero_simd;       /**< Register set to zero (SIMD registers). */
     __m256 gap_open_simd;   /**< Penalty for gap openning (SIMD registers). */
     __m256 gap_extend_simd; /**< Penalty for gap extending (SIMD registers). */
#else
     __m128 zero_simd;       /**< Register set to zero (SIMD registers). */
     __m128 gap_open_simd;   /**< Penalty for gap openning (SIMD registers). */
     __m128 gap_extend_simd; /**< Penalty for gap extending (SIMD registers). */
#endif // __AVX__

     int x_size;   /**< x-length of the H score matrix. */
     int y_size;   /**< y-length of the H score matrix. */
     int max_size; /**< H score matrix size, i.e, x_size * y_size. */
          
     float substitution[2]; /**< Array for storing the match and mismatch penalties. */
     subst_matrix_t *subst_matrix_p; /**< Array for storing the match and mismatch penalties. */
     
     float *E; /**< E vector. */
     float *F; /**< F vector. */
     float *H; /**< H score matrix. */

     int *compass_p; /**< Path direction pointer array, to traceback. */

     int seq_x_end[SIMD_DEPTH]; /**< x-positions of the maximum scores. */
     int ref_y_end[SIMD_DEPTH]; /**< y-positions of the maximum scores. */
     
     float *h_end[SIMD_DEPTH];  /**< pointers to the H cells with the maximum scores. */
     char *seq_end[SIMD_DEPTH]; /**< pointers to the chars of target sequences corresponding to the maximum scores. */
     char *ref_end[SIMD_DEPTH]; /**< pointers to the chars of reference sequences corresponding to the maximum scores. */
     int *compass_end[SIMD_DEPTH]; /**< pointers to the path direction pointer array corresponding to the maximum scores. */
     
     char *a_map; /**< temporary pointer to the aligned target sequence. */
     char *b_map; /**< temporary pointer to the aligned reference sequence. */    
} sw_simd_context_t;

//------------------------------------------------------------------------------------

/**
 * @brief Constructor for the @a sw_simd_context_t structure.
 * @param gap_open penalty for gap openning (e.g., 10.0)
 * @param gap_extend penalty for gap extending (e.g., 0.5)
 * @param subst_matrix_p poiter to the substitution score matrix
 * @return Pointer to the new structure.
 *
 * @a sw_context_see_t constructor that allocates memory for
 * the Smith-Waterman context for the SIMD version. The SIMD context
 * consists of public parameters (match, mismatch, gap_open and gap_extend) and
 * private parameters (E, F, H,..).
 */
sw_simd_context_t *sw_simd_context_new(float gap_open, float gap_extend, subst_matrix_t *subst_matrix_p);

//------------------------------------------------------------------------------------

/**
 * @brief Destructor for the @a sw_simd_context_t structure.
 * @param context_p[out] pointer to the structure to free
 *
 * @a sw_simd_context_t destructor that frees the memory previously
 * allocated by the constructor @a sw_context_see_t new.
 */
void sw_simd_context_free(sw_simd_context_t *context_p);

//------------------------------------------------------------------------------------

/**
 * @brief Updates the SIMD context structure.
 * @param x_size x-length of the score matrix
 * @param y_size y-length of the score matrix
 * @param context_p[out] pointer to the structure to update
 *
 * @a sw_simd_context_t destructor that frees the memory previously
 * allocated by the constructor @a sw_context_see_t new. This is
 * a private function, it should be called by the @a smithwaterman_see
 * function. Mainly, updates the E, F, H vectors sizes and allocate memory
 * for them.
 */
void sw_simd_context_update(int x_size, int y_size, sw_simd_context_t* context_p);

//------------------------------------------------------------------------------------

/**
 * @brief Performs the Smith-Waterman algorithm based-on SIMD instructions.
 * @param input_p pointer to the input sequences to align
 * @param[out] output_p pointer to the output aligned sequences
 * @param[out] context_p pointer to the SIMD context
 *
 * Based-on SIMD instrunctions, this function performs the Smith-Waterman algorithm,
 * it will perform 4 x Smith-Waterman in parallel using the XMM registers.
 */
void smith_waterman_simd(char** query_p, char** ref_p, unsigned int num_queries,
			 char **query_map_p, char **ref_map_p,
			 unsigned int *query_start_p, unsigned int *ref_start_p,
			 float *score_p, sw_simd_context_t* context_p);

#endif // SMITH_WATERMAN_SIMD_H

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
