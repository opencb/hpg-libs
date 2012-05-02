#ifndef VCF_READER_H
#define VCF_READER_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <file_utils.h>
#include <list.h>
#include <log.h>

#include "vcf_util.h"
#include "vcf_file_structure.h"
#include "vcf_file.h"
#include "vcf_read.h"
#include "vcf_batch.h"

extern int mmap_vcf;


enum VCF_Field { CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE };

int vcf_ragel_read(list_t *batches_list, size_t batch_size, vcf_file_t *file, int read_samples);

#endif
