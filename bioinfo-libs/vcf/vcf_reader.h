#ifndef VCF_RAGEL_H
#define VCF_RAGEL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vcf_file_structure.h"
#include "vcf_file.h"
#include "vcf_read.h"
#include "vcf_batch.h"

#include "list.h"

enum Field { CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE };

// int vcf_ragel_read(vcf_file_t *file);
int vcf_ragel_read(list_t *batches_list, size_t batch_size, vcf_file_t *file, int read_samples);

#endif
