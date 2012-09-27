#ifndef CAL_SEEKER_H
#define CAL_SEEKER_H

#include "buffers.h"

#define MAX_CALS 200

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct cal_seeker_input {
  size_t batch_size;
  list_t *write_list;
  list_t *regions_list;
  list_t *sw_list;
  list_t *pair_list;
  cal_optarg_t *cal_optarg;
} cal_seeker_input_t;

void cal_seeker_input_init(list_t *regions_list, cal_optarg_t *cal_optarg, 
			   list_t* write_list, unsigned int write_size, 
			   list_t *sw_list, list_t *pair_list,
			   cal_seeker_input_t *input);

void cal_seeker_server(cal_seeker_input_t* input);

//====================================================================================
// apply_caling
//====================================================================================

void apply_caling(cal_seeker_input_t* input, aligner_batch_t *batch);

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

#endif  // CAL_SEEKER_H
