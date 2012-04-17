#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "vcf_file_structure.h"
#include "vcf_file.h"
#include "vcf_read.h"
#include "vcf_reader.h"

int main (int argc, char *argv[])
{
	double start, stop, total;
	
	// Reading
	start = omp_get_wtime();
	
	vcf_file_t* file = vcf_open(argv[1]);
	int ret_code = vcf_read(file);
	
	stop = omp_get_wtime();
	total = (stop - start);
	
	if (ret_code) dprintf("[R] Error code = %d\n", ret_code);
	bprintf("[R] Time elapsed = %f s\n", total);
	bprintf("[R] Time elapsed = %e ms\n", total*1000);

	// Writing to a new file
	if (argc == 3) 
	{
		start = omp_get_wtime();
	
		ret_code = vcf_write(file, argv[2]);
		
		stop = omp_get_wtime();
		total = (stop - start);
		
		if (ret_code) dprintf("[W] Error code = %d\n", ret_code);
		bprintf("[W] Time elapsed = %f s\n", total);
		bprintf("[W] Time elapsed = %e ms\n", total*1000);
	}
	
	vcf_close(file);
	
	return 0;
}
