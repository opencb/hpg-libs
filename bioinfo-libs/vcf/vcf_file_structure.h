#ifndef VCF_FILE_STRUCTURE_H
#define VCF_FILE_STRUCTURE_H

#include <sys/types.h>

#include <containers/array_list.h>
#include <containers/list.h>

/**
 * Entry in the VCF document header.
 */
typedef struct vcf_header_entry {
    char *name;
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

    array_list_t *samples;
} vcf_record_t;

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
    array_list_t *header_entries;
    size_t num_header_entries;
    array_list_t *samples_names;
    size_t num_samples;
    array_list_t *records;
    size_t num_records;
} vcf_file_t;


#endif
