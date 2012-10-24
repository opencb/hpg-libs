#include "fastq_stats.h"

fastq_read_stats_t *fastq_stats_new() {
	fastq_read_stats_t *read_stats = (fastq_read_stats_t *)malloc(sizeof(fastq_read_stats_t));

	read_stats->length = 0;
	read_stats->quality_average = 0.0;
	read_stats->Ns = 0;
	read_stats->N_out_quality = 0;

	read_stats->num_A = 0;
	read_stats->num_T = 0;
	read_stats->num_C = 0;
	read_stats->num_G = 0;
	read_stats->num_N = 0;

	return read_stats;
}

//void fastq_stats_init(int length, float qual_average, int Ns, int N_out_qual, fastq_read_stats_t *read_stats) {
//	if(read_stats == NULL) {
//		return;
//	}
//	read_stats->length = length;
//	read_stats->quality_average = qual_average;
//	read_stats->Ns = Ns;
//	read_stats->N_out_quality = N_out_qual;
//}

void fastq_stats_free(fastq_read_stats_t *read_stats) {
	if(read_stats != NULL) {
		free(read_stats);
	}
}

size_t fastq_stats(array_list_t *fq_reads, array_list_t *reads_stats) {
	array_list_new(fq_reads->capacity, fq_reads->realloc_factor, fq_reads->mode);
	if(fq_reads != NULL) {

	}
	return reads_stats->size;
}

//int fastq_stats(fastq_read_t *fq_read, fastq_read_stats_t *read_stats) {
//	for(int i=0; i<fq_read->length; i++) {
//
//	}
//	return 0;
//}
