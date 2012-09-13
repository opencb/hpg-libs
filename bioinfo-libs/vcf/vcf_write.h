#ifndef VCF_WRITE_H
#define VCF_WRITE_H

#include <stdio.h>
#include <string.h>

#include <containers/array_list.h>

#include "vcf_file_structure.h"
#include "vcf_batch.h"

//====================================================================================
//  vcf_write.h
//
//  vcf writing functions prototypes
//====================================================================================

int vcf_write_to_file(vcf_file_t *vcf_file, FILE *fd);

/* ************ Header management functions **********************/

void write_file_format(vcf_file_t *vcf_file, FILE *fd);

void write_header_entry(vcf_header_entry_t *entry, FILE *fd);

/* ************ Record management functions **********************/

void write_delimiter(vcf_file_t *vcf_file,FILE *fd);

void write_batch(vcf_batch_t *vcf_batch, FILE *fd);

void write_record(vcf_record_t* vcf_record, FILE *fd);


#endif 
