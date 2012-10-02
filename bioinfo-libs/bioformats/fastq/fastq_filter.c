/*
 * fastq_filter.c
 *
 *  Created on: Sep 27, 2012
 *      Author: imedina
 */

#include "fastq_read.h"

filter_options_t *fastq_filter_new(int min_length, int max_lentgh, float min_quality, float max_quality, int max_Ns, int max_N_out_quality) {
	filter_options_t *options = (filter_options_t *)malloc(size_of(filter_options_t));

	options->min_length = min_length;
	options->max_length = max_lentgh;
	options->min_quality = min_quality;
	options->max_quality = max_quality;
	options->max_Ns = max_Ns;
	options->max_N_out_quality = max_N_out_quality;

	return options;
}

void fastq_filter_free(filter_options_t *options) {
	if(options != NULL) {
		free(options);
	}
}

array_list_t *fastq_filter(array_list_t *reads, array_list_t *passed, array_list_t *failed, filter_options_t *options) {
	fastq_read_t *read;
	for(size_t i=0; i<reads->size; i++) {
		read = array_list_get(i, reads);
		if(read->length < options->min_length || read->length > options->max_length) {
			array_list_insert(read, passed);
		}else {
			// get read stats

		}
	}
	return passed;
}
