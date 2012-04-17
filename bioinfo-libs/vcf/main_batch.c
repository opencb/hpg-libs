#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "vcf_file_structure.h"
#include "vcf_file.h"
#include "vcf_read.h"
#include "vcf_reader.h"

#include "list.h"

int main (int argc, char *argv[])
{
	size_t max_batches = 20;
	size_t batch_size = 2000;
	list_t *read_list = (list_t*) malloc (sizeof(list_t));
	list_init("batches", 1, max_batches, read_list);
	
	int ret_code;
	double start, stop, total;
	vcf_file_t* file;
	
#pragma omp parallel sections private(start, stop, total) lastprivate(file)
{
	#pragma omp section
	{
		dprintf("Thread %d reads the VCF file\n", omp_get_thread_num());
		// Reading
		start = omp_get_wtime();
		
		file = vcf_open(argv[1]);
		ret_code = vcf_read_batches(read_list, batch_size, file, 1);
		
		stop = omp_get_wtime();
		total = (stop - start);
		
		if (ret_code) dprintf("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code);
		bprintf("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
		bprintf("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
		
		// Writing to a new file
		if (argc == 3) 
		{
			start = omp_get_wtime();
		
			ret_code = vcf_write(file, argv[2]);
			
			stop = omp_get_wtime();
			total = (stop - start);
			
			if (ret_code) dprintf("[%dW] Error code = %d\n", omp_get_thread_num(), ret_code);
			bprintf("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), total);
			bprintf("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
		}
		
		list_decr_writers(read_list);
		
		vcf_close(file);
	
	}
	#pragma omp section
	{
		dprintf("OMP num threads = %d\n", omp_get_num_threads());
		dprintf("Thread %d prints info\n", omp_get_thread_num());
		
		start = omp_get_wtime();
		
		int i = 0;
		list_item_t* item = NULL;
		while ( (item = list_remove_item(read_list)) != NULL ) {
// 			if (i % 200 == 0) 
// 			{
				int debug = 1;
				dprintf("Batch %d reached by thread %d - %zu/%zu records \n", i, omp_get_thread_num(), 
					((vcf_batch_t*) item->data_p)->length, ((vcf_batch_t*) item->data_p)->max_length);
// 			}
			
// 			vcf_batch_print(stdout, item->data_p);
			vcf_batch_free(item->data_p);
			list_item_free(item);
			i++;
		}
		
		stop = omp_get_wtime();
		total = (stop - start);
		
		bprintf("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
		bprintf("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
	}
}
	
	free(read_list);
	
	return 0;
}
