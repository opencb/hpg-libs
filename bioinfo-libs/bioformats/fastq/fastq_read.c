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

fastq_read_stats_t *fastq_read_stats_new() {
	fastq_read_stats_t *read_stats = (fastq_read_stats_t *)malloc(size_of(fastq_read_stats_t));
	return read_stats;
}

void fastq_read_stats_init(int length, float qual_average, int Ns, int N_out_qual, fastq_read_stats_t *read_stats) {
	if(read_stats == NULL) {
		return;
	}
	read_stats->length = length;
	read_stats->quality_average = qual_average;
	read_stats->Ns = Ns;
	read_stats->N_out_quality = N_out_qual;
}

void fastq_read_stats_free(fastq_read_stats_t *read_stats) {
	if(read_stats != NULL) {
		free(read_stats);
	}
}

int fastq_read_stats(fastq_read_t *fq_read, fastq_read_stats_t *read_stats) {
	for(int i=0; i<fq_read->length; i++) {

	}
	return 0;
}

int fastq_reads_stats(array_list_t *fq_reads, array_list_t *reads_stats) {

	return 0;
}


