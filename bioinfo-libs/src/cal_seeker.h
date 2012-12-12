#ifndef CAL_SEEKER_H
#define CAL_SEEKER_H

#include "buffers.h"
#include "aligners/bwt/bwt.h"

#define MAX_CALS 200
#define MAX_RNA_CALS 100

//====================================================================================
//  structures and prototypes
//====================================================================================

/**
 * @brief Structure for store all parameters needed for run @a cal_seeker_server.
 * 
 * Structure for store some configuration values and data structures like lists 
 * to insert and read data batches.
 */
typedef struct cal_seeker_input {
  unsigned batch_size;             /**< size of data batches*/
  list_t *write_list;            /**< list for store write batches*/
  list_t *regions_list;          /**< list for read batches with all regions found for each read */
  list_t *sw_list;               /**< list to store batches with all CALs found for each read */
  list_t *pair_list;
  cal_optarg_t *cal_optarg;      /**< cal seeker configuration values */
} cal_seeker_input_t;

/**
 * @brief  Initializer for the @a cal_seeker_input_t structure.
 * @param  region_list_p list for read batches with all regions found for each read
 * @param  cal_optarg_p cal seeker configuration values 
 * @param  write_list_p list for store write batches
 * @param  write_size size of write batches
 * @param  sw_list_p list to store batches with all CALs found for each read
 * @param  pair_list_p list to store batches with pairs
 *
 * 
 * Initialize all @a cal_seeker_input_t fields with the input parameters.
 */
void cal_seeker_input_init(list_t *regions_list, cal_optarg_t *cal_optarg, 
			   list_t* write_list, unsigned int write_size, 
			   list_t *sw_list, list_t *pair_list,
			   cal_seeker_input_t *input);

//--------------------------------------------------------------------------------------

/**
 * @brief CAL Seeker server make CALs from regions
 * @param input_p all configuration values and data structures needed by the server
 *
 * CAL Seeker server extract batches of reads with their regions from @a regions_list_p and
 * for each read call @a bwt_generate_cal_list for generate CALs. Finally, all reads with a lot of
 * CALs are stored in @a write_list_p because those reads aren't mapped and will be store in bam file,
 * the rest of reads are stored with their CALs in @a sw_list_p.      
 */
void cal_seeker_server(cal_seeker_input_t* input);

//====================================================================================
// apply_caling
//====================================================================================

void apply_caling(cal_seeker_input_t* input, mapping_batch_t *batch);

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

#endif  // CAL_SEEKER_H
