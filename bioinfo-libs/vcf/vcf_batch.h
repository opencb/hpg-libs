#ifndef VCF_BATCH_H
#define VCF_BATCH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "vcf_file_structure.h"
#include "util.h"

//====================================================================================
//  vcf_batch.h
//
//  vcf_batch_t structures and prototypes
//====================================================================================

/**
 * Struct which represents a batch of VCF records whose fields have been 
 * already loaded.
 */
typedef list_t vcf_batch_t;
// typedef struct vcf_batch {
//        size_t num_records;	// Number of records read
//        size_t size;             // Max buffer size (can be greater than num_records)
//        vcf_record_t **records;	// Records read
// } vcf_batch_t;

vcf_batch_t* vcf_batch_new(size_t size);

void vcf_batch_free(vcf_batch_t *vcf_batch);

void add_record_to_batch(vcf_record_t *record, vcf_batch_t *vcf_batch);

int batch_is_empty(vcf_batch_t *vcf_batch);

int batch_is_full(vcf_batch_t *vcf_batch);

void vcf_batch_print(FILE *fd, vcf_batch_t *vcf_batch);

#endif
