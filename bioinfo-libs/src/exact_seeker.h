
#ifndef EXACT_SEEKER_H
#define EXACT_SEEKER_H

#include "buffers.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct exact_seeker_input {
  unsigned int num_seeds;
  unsigned int seed_size;
  unsigned int write_size;

  list_t* kl_list_p;
  list_t* seed_list_p;
  list_t* write_list_p;
  bwt_index_t* index_p;
  genome_t* genome_p;
} exact_seeker_input_t;

//------------------------------------------------------------------------------------

void exact_seeker_input_init(list_t* kl_list_p, list_t* split_read_list_p, list_t* write_list_p, unsigned int num_seeds, unsigned int seed_size, unsigned int write_size, bwt_index_t* index_p, genome_t* genome_p, exact_seeker_input_t* input_p);

//------------------------------------------------------------------------------------
// main function
//------------------------------------------------------------------------------------

void exact_seeker(exact_seeker_input_t* input_p);

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif  // EXACT_SEEKER_H
