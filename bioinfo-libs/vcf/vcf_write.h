#ifndef VCF_WRITE_H
#define VCF_WRITE_H

#include <stdio.h>
#include <string.h>

#include <containers/array_list.h>

#include "vcf_file_structure.h"
// #include "vcf_batch.h"

//====================================================================================
//  vcf_write.h
//
//  vcf writing functions prototypes
//====================================================================================

int write_vcf_file(vcf_file_t *vcf_file, FILE *fd);

/* ************ Header management functions **********************/

void write_vcf_header(vcf_file_t *vcf_file, FILE *fd);

void write_vcf_fileformat(vcf_file_t *vcf_file, FILE *fd);

void write_vcf_header_entry(vcf_header_entry_t *entry, FILE *fd);

void write_vcf_delimiter(vcf_file_t *vcf_file,FILE *fd);

/* ************ Record management functions **********************/

void write_vcf_batch(vcf_batch_t *vcf_batch, FILE *fd);

void write_vcf_record(vcf_record_t* vcf_record, FILE *fd);


#endif 
