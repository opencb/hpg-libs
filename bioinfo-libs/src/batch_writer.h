#ifndef BATCH_WRITER_H
#define BATCH_WRITER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "cprops/hashtable.h"

#include "commons/commons.h"
#include "commons/system_utils.h"
#include "containers/list.h"

#include "bioformats/fastq/fastq_file.h"
#include "bioformats/fastq/fastq_batch.h"
#include "bioformats/bam-sam/bam_file.h"

#include "buffers.h"
#include "timing.h"
#include "sw_server.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct batch_writer_input {
  char* match_filename;
  char* mismatch_filename;

  char* splice_exact_filename;
  char* splice_extend_filename;
  
  char* header_filename;

  list_t* list_p;
} batch_writer_input_t;

//------------------------------------------------------------------------------------

void batch_writer_input_init(char* match_filename, char* splice_exact_filename, 
			     char* splice_extend_filename, list_t* list_p, 
			     char* header_filename, batch_writer_input_t* input);

//====================================================================================

void batch_writer(batch_writer_input_t* input_p);
void batch_writer2(batch_writer_input_t* input_p);

//void batch_writer_splice(batch_writer_splice_input_t* input_p);
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif // BATCH_WRITER_H
