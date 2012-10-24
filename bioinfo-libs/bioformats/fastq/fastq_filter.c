#include <stdlib.h>

#include "fastq_filter.h"
#include "fastq_read.h"


fastq_filter_options_t *fastq_filter_options_new(int min_length, int max_lentgh, float min_quality, float max_quality, int max_Ns, int max_N_out_quality) {
	fastq_filter_options_t *options = (fastq_filter_options_t *)malloc(sizeof(fastq_filter_options_t));

	options->min_length = min_length;
	options->max_length = max_lentgh;
	options->min_quality = min_quality;
	options->max_quality = max_quality;
	options->max_Ns = max_Ns;
	options->max_N_out_quality = max_N_out_quality;

	return options;
}

void fastq_filter_options_free(fastq_filter_options_t *options) {
	if(options != NULL) {
		free(options);
	}
}

array_list_t *fastq_filter(array_list_t *reads, array_list_t *passed, array_list_t *failed, fastq_filter_options_t *options) {
	fastq_read_t *read;
	float qual_avg;
//	#pragma omp parallel for schedule(dynamic, 1000000) private(qual_avg)
	for(size_t i=0; i<reads->size; i++) {
//		printf("%i\n", omp_get_thread_num( ));
		read = array_list_get(i, reads);
		if(read->length >= options->min_length && read->length <= options->max_length) {
			qual_avg = fastq_quality_average(read);
			if(qual_avg >= options->min_quality && qual_avg <= options->max_quality) {
//				printf("qual avg: %f\n", qual_avg);
				array_list_insert(read, passed);
			}

		}else {
			// get read stats
			array_list_insert(read, failed);
		}
	}
	return passed;
}
