#ifndef VCF_READER_H
#define VCF_READER_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/list.h>

#include "vcf_util.h"
#include "vcf_file_structure.h"
#include "vcf_file.h"
#include "vcf_read.h"
#include "vcf_batch.h"

#define CHUNK 0x80000

extern int mmap_vcf;


typedef struct {
    vcf_record_t *current_record;
    vcf_header_entry_t *current_header_entry;
    vcf_batch_t *current_batch;
    
    size_t num_samples;
    size_t num_records;
    
    int store_samples;
    int self_contained; // TODO this could not be neccessary, need to check
} vcf_reader_status;


int execute_vcf_ragel_machine(char *p, char *pe, list_t *batches_list, size_t batch_size, vcf_file_t *file, vcf_reader_status *status);

vcf_reader_status *vcf_reader_status_new(size_t batch_size, int store_samples, int self_contained);

void vcf_reader_status_free(vcf_reader_status *status);


int vcf_read_and_parse(list_t *batches_list, size_t batch_size, vcf_file_t *file, int read_samples);

int vcf_gzip_read_and_parse(list_t *batches_list, size_t batch_size, vcf_file_t *file, int read_samples);


int vcf_light_read(list_t *batches_list, size_t batch_size, vcf_file_t *file);

int vcf_gzip_light_read(list_t *batches_list, size_t batch_size, vcf_file_t *file);


int vcf_light_multiread(list_t **batches_list, size_t batch_size, vcf_file_t **files, size_t num_files);


#endif
