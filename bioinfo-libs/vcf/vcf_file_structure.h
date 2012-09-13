#ifndef VCF_FILE_STRUCTURE_H
#define VCF_FILE_STRUCTURE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/types.h>

#include <commons/log.h>
#include <containers/array_list.h>
#include <containers/list.h>

/**
 * VCF file structure. The physical file is defined by its file descriptor, its 
 * filename and the mode it has been open.
 *
 * It contains a header with several entries, and a body with several records.
 */
typedef struct vcf_file {
    char* filename;
    char* mode;

    FILE *fd;   /**< Should the file be loaded using IO functions, the file descriptor is set */

    char *data;   /**< Should the file be loaded using mmap, its contents are set */
    size_t data_len;

    char* format;
    int format_len;
    
    array_list_t *header_entries;
    size_t num_header_entries;
    array_list_t *samples_names;
    size_t num_samples;
    array_list_t *records;
    size_t num_records;
} vcf_file_t;

/**
 * Entry in the VCF document header.
 */
typedef struct vcf_header_entry {
    char *name;
    int name_len;
    
    // List of keys of the fields describing the entry (if num_field_keys = 0, then entry = "name=values[0]")
    list_t *keys;
    size_t num_keys;
    // List of values of the fields describing the entry
    list_t *values;
    size_t num_values;
} vcf_header_entry_t;


/**
 * Entry in the VCF document body.
 */
typedef struct vcf_record {
    char* chromosome;
    unsigned long position;
    char* id;
    char* reference;
    char* alternate;
    float quality;
    char* filter;
    char* info;
    char* format;

    int chromosome_len;
    int id_len;
    int reference_len;
    int alternate_len;
    int filter_len;
    int info_len;
    int format_len;
    
    array_list_t *samples;
} vcf_record_t;



/*
 * Creation and destruction of header entries and records
 */

vcf_header_entry_t* create_header_entry();

void vcf_header_entry_free(vcf_header_entry_t *vcf_header_entry);

vcf_record_t* create_record();

void vcf_record_free(vcf_record_t *vcf_record);

/*
 * Addition of header and record entries
 */

int add_header_entry(vcf_header_entry_t *header_entry, vcf_file_t *vcf_file);

int add_sample_name(char *name, int length, vcf_file_t *vcf_file);

int add_record(vcf_record_t* record, vcf_file_t *vcf_file);



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
