#ifndef VCF_FILE_H
#define VCF_FILE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>

#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/array_list.h>
#include <containers/list.h>

#include "vcf_util.h"
#include "vcf_file_structure.h"
#include "vcf_reader.h"
#include "vcf_write.h"

#define INIT_RECORD_SIZE	100

extern int mmap_vcf;

//====================================================================================
//  vcf_file.h
//
//  vcf_file_t structures and prototypes
//====================================================================================


//====================================================================================
//  Function prototypes
//====================================================================================

/*
 * Physical file management
 */

/**
 * Open a file stream in a given file mode (r/w/a) and initialize the associated 
 * vcf_file_t structure.
 * 
 * @param filename The name of the file to open
 * @param mode Open mode (read/write/append)
 */
vcf_file_t *vcf_open(char *filename);

/**
 * Close the file stream associated to a vcf_file_t type.
 * 
 * @param vcf_file vcf_file_t whose file stream is about to being closed
 */
void vcf_close(vcf_file_t *vcf_file);


/*
 * Creation and destruction of header entries and records
 */

vcf_header_entry_t* create_header_entry();

void vcf_header_entry_free(vcf_header_entry_t *vcf_header_entry);

vcf_record_t* create_record();

void vcf_record_free(vcf_record_t *vcf_record);

/*
 * File reading
 */

/**
 * Fill the fields of the vcf_file_t given as argument reading data from a file.
 * 
 * TODO breaks because a batches list isn't provided! 
 * a solution would be creating a temp batches list, then copying its contents to the vcf_file_t
 * 
 * @param vcf_file The vcf_file_t whose fields will be set
 */
int vcf_read(vcf_file_t *vcf_file);

int vcf_read_batches(list_t *batches_list, size_t batch_size, vcf_file_t *vcf_file, int read_samples);

int vcf_parse_batches(list_t *batches_list, size_t batch_size, vcf_file_t *vcf_file, int read_samples);

/**
 * Write the contents of the vcf_file_t given as argument to the given path.
 * 
 * @param vcf_file The vcf_file_t whose information will be written
 * @param filename The path of the file to write the information to
 */
int vcf_write(vcf_file_t *vcf_file, char *filename);

/*
 * Header and record entries management
 */

int add_header_entry(vcf_header_entry_t *header_entry, vcf_file_t *vcf_file);

int add_sample_name(char *name, vcf_file_t *vcf_file);

int add_record(vcf_record_t* record, vcf_file_t *vcf_file);

#endif
