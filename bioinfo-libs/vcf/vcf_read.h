#ifndef VCF_READ_H
#define VCF_READ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <commons/log.h>
#include <containers/array_list.h>
#include <containers/list.h>

#include "vcf_file_structure.h"

//====================================================================================
//  vcf_read.h
//
//  vcf reading functions prototypes
//====================================================================================



/* ************ Header management functions **********************/

void set_file_format(char *fileformat, int length, vcf_file_t *vcf_file);

void set_header_entry_name(char *name, int length, vcf_header_entry_t *entry);

void add_header_entry_key(char *key, int length, vcf_header_entry_t *entry);

void add_header_entry_value(char *value, int length, vcf_header_entry_t *entry);

/* ************ Record management functions **********************/

void set_record_chromosome(char* chromosome, int length, vcf_record_t* vcf_record);

void set_record_position(long position, vcf_record_t* vcf_record);

void set_record_id(char* id, int length, vcf_record_t* vcf_record);

void set_record_reference(char* reference, int length, vcf_record_t* vcf_record);

void set_record_alternate(char* alternate, int length, vcf_record_t* vcf_record);

void set_record_quality(float quality, vcf_record_t* vcf_record);

void set_record_filter(char* filter, int length, vcf_record_t* vcf_record);

void set_record_info(char* info, int length, vcf_record_t* vcf_record);

void set_record_format(char* format, int length, vcf_record_t* vcf_record);

void add_record_sample(char* sample, int length, vcf_record_t* vcf_record);

#endif
