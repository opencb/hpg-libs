
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fastq_read.h"


//--------------------------------------------
//--------------------------------------------

fastq_read_t *fastq_read_new(char *id, char *sequence, char *quality) {

  fastq_read_t *fq_read = (fastq_read_t*) calloc(1, sizeof(fastq_read_t));

  fq_read->id = strdup(id);
  fq_read->sequence = strdup(sequence);
  fq_read->quality = strdup(quality);

  return fq_read;
}

//--------------------------------------------

void fastq_read_free(fastq_read_t *fq_read) {
  if (fq_read == NULL) return;
  
  if (fq_read->id != NULL) free(fq_read->id);
  if (fq_read->sequence != NULL) free(fq_read->sequence);
  if (fq_read->quality != NULL) free(fq_read->quality);

  free(fq_read);
}


//--------------------------------------------
//--------------------------------------------

