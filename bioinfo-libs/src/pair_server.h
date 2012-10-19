#ifndef PAIR_SERVER_H
#define PAIR_SERVER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "commons/commons.h"
#include "commons/system_utils.h"
#include "containers/list.h"

#include "buffers.h"
#include "timing.h"
#include "sw_server.h"
//#include "hashtable.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct pair_server_input {
  pair_mng_t *pair_mng;

  list_t* pair_list;
  list_t* sw_list;
  list_t* write_list;
} pair_server_input_t;

//------------------------------------------------------------------------------------

void pair_server_input_init(pair_mng_t *pair_mng, list_t* pair_list, list_t *sw_list,
			    list_t *write_list, pair_server_input_t* input);

//====================================================================================

void pair_server(pair_server_input_t* input);

//------------------------------------------------------------------------------------

void apply_pair(pair_server_input_t* input, aligner_batch_t *batch);
void prepare_alignments(pair_server_input_t* input, aligner_batch_t *batch);

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif // PAIR_SERVER_H
