#ifndef BATCH_WRITER_H
#define BATCH_WRITER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "commons.h"
#include "list.h"
#include "system_utils.h"
#include "buffers.h"
#include "timing.h"
#include "fastq_file.h"
#include "fastq_batch.h"
#include "bam_file.h"
//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct batch_writer_input {
  char* match_filename;
  char* mismatch_filename;
  char* splice_exact_filename;
  char* splice_extend_filename;
  list_t* list_p;
} batch_writer_input_t;

//------------------------------------------------------------------------------------

void batch_writer_input_init(char* match_filename, char* splice_exact_filename, char* splice_extend_filename, list_t* list_p, batch_writer_input_t* input);
//====================================================================================

void batch_writer(batch_writer_input_t* input_p);
//void batch_writer_splice(batch_writer_splice_input_t* input_p);
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif    // BATCH_WRITER_H
