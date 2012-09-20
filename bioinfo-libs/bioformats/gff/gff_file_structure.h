#ifndef GFF_FILE_STRUCTURE_H
#define GFF_FILE_STRUCTURE_H

#include <sys/types.h>

#include <cprops/linked_list.h>

#include <commons/file_utils.h>

/**
 * Entry in the GFF document header.
 */
typedef struct gff_header_entry {
    char *text;
} gff_header_entry_t;


/**
 * Entry in the GFF document body.
 */
typedef struct gff_record {
    char *sequence;
    char *source;
    char *feature;
    char *attribute;
    unsigned long start;
    unsigned long end;
    int16_t score;
    int8_t frame;
    char strand;
} gff_record_t;


/**
 * GFF file structure. The physical file is defined by its file descriptor, its 
 * filename and the mode it has been open.
 *
 * It contains a header with several entries, and a body with several records.
 * For last version of the spec see: http://www.sanger.ac.uk/resources/software/gff/spec.html
 */
typedef struct gff_file {
    char* filename;
    char* mode;
    char *data;
    size_t data_len;
    
    cp_list *header_entries;
    cp_list *records;
} gff_file_t;

#endif
