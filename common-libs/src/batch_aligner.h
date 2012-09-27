#ifndef BATCH_ALIGNER_H
#define BATCH_ALIGNER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "commons/commons.h"
#include "commons/system_utils.h"
#include "containers/list.h"

#include "buffers.h"
#include "timing.h"
#include "bwt_server.h"
#include "region_seeker.h"
#include "cal_seeker.h"
#include "pair_server.h"
#include "sw_server.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct batch_aligner_input {
  list_t *read_list;
  list_t *write_list;

  bwt_server_input_t *bwt_input;
  region_seeker_input_t *region_input;
  cal_seeker_input_t *cal_input;
  pair_server_input_t *pair_input;
  sw_server_input_t *sw_input;
} batch_aligner_input_t;

//------------------------------------------------------------------------------------

void batch_aligner_input_init(list_t *read_list, list_t *write_list, 
			      bwt_server_input_t *bwt_input, 
			      region_seeker_input_t *region_input,
			      cal_seeker_input_t *cal_input,
			      pair_server_input_t *pair_input,
			      sw_server_input_t *sw_input,
			      batch_aligner_input_t *input);

//====================================================================================
//  main batch aligner function
//====================================================================================

void batch_aligner(batch_aligner_input_t *input);

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif // BATCH_ALIGNER_H
