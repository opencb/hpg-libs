#include <stdio.h>

#include "vcf_file_structure.h"

int main (int argc, char *argv[])
{
	printf("Size of vcf_file = %zu\n", sizeof(vcf_file_t));
	printf("Size of vcf_header_entry = %zu\n", sizeof(vcf_header_entry_t));
	printf("Size of vcf_record = %zu\n", sizeof(vcf_record_t));
	return 0;
}
