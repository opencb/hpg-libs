#include <string.h>
#include <stdlib.h>

#include "fastq_read.h"

//--------------------------------------------
//--------------------------------------------

fastq_read_t *fastq_read_new(char *id, char *sequence, char *quality) {

  fastq_read_t *fq_read = (fastq_read_t*) calloc(1, sizeof(fastq_read_t));
  /*fq_read->id = strdup(id);
  fq_read->quality = strdup(quality);
  fq_read->sequence = strdup(sequence);*/
  
  size_t id_length = strlen(id); // + 1;
  size_t seq_length = strlen(sequence); // + 1;
  
  fq_read->id = strndup(id, id_length);
  fq_read->sequence = strndup(sequence, seq_length);
  fq_read->quality = strndup(quality, seq_length);
  fq_read->length = seq_length;

//  fq_read->id = (char *)malloc(sizeof(char)*id_length);
//  fq_read->sequence = (char *)malloc(sizeof(char)*seq_length);
//  fq_read->quality = (char *)malloc(sizeof(char)*seq_length);
//
//  memcpy(fq_read->id, id, id_length);
//  memcpy(fq_read->sequence, sequence, seq_length);
//  memcpy(fq_read->quality, quality, seq_length);
  
  return fq_read;
}

void fastq_read_free(fastq_read_t *fq_read) {
  if (fq_read == NULL) return;
  
  if (fq_read->id != NULL) free(fq_read->id);
  if (fq_read->sequence != NULL) free(fq_read->sequence);
  if (fq_read->quality != NULL) free(fq_read->quality);

  free(fq_read);
}


void fastq_read_print(fastq_read_t *read) {
	printf("%s\n", read->id);
	printf("%s\n", read->sequence);
	printf("+\n");
	printf("%s\n", read->quality);
}

//float fastq_quality_average(fastq_read_t *fq_read_t) {
//	int qual = 0;
//	for(size_t i=0; i<fq_read_t->length; i++) {
//		qual += fq_read_t->quality[i];
//	}
//	return (float)qual/fq_read_t->length;
//}

