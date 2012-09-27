#ifndef REGION_SEEKER_H
#define REGION_SEEKER_H

#include <stdio.h>

#include "containers/list.h"

#include "buffers.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct region_seeker_input{
  list_t *unmppaed_reads_list_p;
  list_t* region_list_p;
  cal_optarg_t *cal_optarg_p;
  bwt_optarg_t *bwt_optarg_p;
  bwt_index_t *bwt_index_p;
}region_seeker_input_t;

void region_seeker_input_init(list_t *unmppaed_reads_list_p, cal_optarg_t *cal_optarg_p, bwt_optarg_t *bwt_optarg_p, bwt_index_t *bwt_index_p, list_t* region_list_p, region_seeker_input_t *input_p);

//--------------------------------------------------------------------------------------

void region_seeker_server(region_seeker_input_t *input);

//====================================================================================

void apply_seeding(region_seeker_input_t* input, aligner_batch_t *batch);

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

#endif // REGION_SEEKER_H
