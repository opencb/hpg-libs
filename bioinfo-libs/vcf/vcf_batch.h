#ifndef VCF_BATCH_H
#define VCF_BATCH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <containers/array_list.h>
#include <commons/log.h>

#include "vcf_file_structure.h"

//====================================================================================
//  vcf_batch.h
//
//  vcf_batch_t structures and prototypes
//====================================================================================

/**
 * List equivalent to a batch of VCF records.
 */
typedef struct vcf_batch {
    array_list_t *records;
    char *text;
} vcf_batch_t;


vcf_batch_t* vcf_batch_new(size_t size);

void vcf_batch_free(vcf_batch_t *vcf_batch);

void add_record_to_vcf_batch(vcf_record_t *record, vcf_batch_t *vcf_batch);


int vcf_batch_is_empty(vcf_batch_t *vcf_batch);

int vcf_batch_is_full(vcf_batch_t *vcf_batch);


void vcf_batch_print(FILE *fd, vcf_batch_t *vcf_batch);

#endif
