#ifndef SW_SERVER_H
#define SW_SERVER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "commons/commons.h"
#include "commons/system_utils.h"
#include "containers/list.h"

#include "aligners/sw/smith_waterman.h"

#include "timing.h"
#include "buffers.h"
#include "genome.h"

//====================================================================================
//  Input structure for Smith-Waterman server
//====================================================================================

/**
 * @brief Input structure for the Smith-Waterman server function.
 *
 * This structure contains all required parameters by
 * the Smith-Waterman server function (@see sw_server).
 */
typedef struct sw_server_input {
     float match;      /**< Penalty for match. */
     float mismatch;   /**< Penalty for mismatch. */
     float gap_open;   /**< Penalty for gap opening. */
     float gap_extend; /**< Penalty for gap extending. */
     float min_score;  /**< Minimum normalized score (0..1) for valid alignments. */
     unsigned int flank_length; /**< Length to extend the CAL region. */
     unsigned int write_size;   /**< Size of the writing batch (to disk). */

     list_t* sw_list_p;    /**< Pointer to the list that contains the input sequences to align. */
     list_t* write_list_p; /**< Pointer to the list that contains the output aligned sequences. */

     genome_t* genome_p;   /**< Pointer to the genome structure to get the reference sequences. */
} sw_server_input_t;

//------------------------------------------------------------------------------------

/**
 * @brief Initialization function for the @a sw_server_input_t structure.
 * @param sw_list_p pointer to the list that contains the input sequences to align
 * @param write_list_p pointer to the list that contains the output aligned sequences
 * @param write_size Size of the writing batch (to disk)
 * @param match Penalty for match
 * @param mismatch Penalty for mismatch
 * @param gap_open Penalty for gap opening
 * @param gap_extend Penalty for gap extending
 * @param min_score Minimum normalized score (0..1) for valid alignments
 * @param flank_length Length to extend the CAL region
 * @param genome_p pointer to the genome structure to get the reference sequences
 * @param[out] input_p pointer to the structure to initialize
 *
 * This function takes the input parameters and initializes the @a sw_server_input_t
 * structure that will be used by the Smith-Waterman server function (@see sw_server).
 */
void sw_server_input_init(list_t* sw_list_p, list_t* write_list_p, unsigned int write_size, 
			  float match, float mismatch, float gap_open, float gap_extend, 
			  float min_score, unsigned int flank_length, genome_t* genome_p,
			  sw_server_input_t* input_p);


//====================================================================================
// apply_sw
//====================================================================================

typedef struct sw_output {
  int strand;
  size_t chromosome;
  size_t ref_start;
  size_t ref_len;
  size_t mref_len;
  size_t mquery_start;
  size_t mref_start;
  float score;
  float norm_score;
  char* mquery;
  char* mref;
} sw_output_t;

sw_output_t *sw_output_new(int strand, size_t chrom, size_t ref_start, size_t ref_len,
			   size_t mref_len, size_t mquery_start, size_t mref_start,
			   float score, float norm_score, char* mquery, char* mref);

void sw_output_free(sw_output_t *p);

//--------------------------------------------------------------------------------------

void apply_sw(sw_server_input_t* input, aligner_batch_t *batch);

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

#endif  // SW_SERVER_H
