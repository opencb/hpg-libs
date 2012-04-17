#ifndef _VCF_READ_H
#define _VCF_READ_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "vcf_file_structure.h"
#include "util.h"
#include "list.h"

//====================================================================================
//  vcf_read.h
//
//  vcf reading functions prototypes
//====================================================================================



/* ************ Header management functions **********************/

void set_file_format(char *fileformat, vcf_file_t *vcf_file);

vcf_header_entry_t* create_header_entry();

void set_header_entry_name(char *name, vcf_header_entry_t *entry);

void add_header_entry_key(char *key, vcf_header_entry_t *entry);

void add_header_entry_value(char *value, vcf_header_entry_t *entry);

/* ************ Record management functions **********************/

vcf_record_t* create_record();

void set_record_chromosome(char* chromosome, vcf_record_t* vcf_record);

void set_record_position(long position, vcf_record_t* vcf_record);

void set_record_id(char* id, vcf_record_t* vcf_record);

void set_record_reference(char* reference, vcf_record_t* vcf_record);

void set_record_alternate(char* alternate, vcf_record_t* vcf_record);

void set_record_quality(float quality, vcf_record_t* vcf_record);

void set_record_filter(char* filter, vcf_record_t* vcf_record);

void set_record_info(char* info, vcf_record_t* vcf_record);

void set_record_format(char* format, vcf_record_t* vcf_record);

void add_record_sample(char* sample, vcf_record_t* vcf_record, size_t *sample_idx);

#endif
