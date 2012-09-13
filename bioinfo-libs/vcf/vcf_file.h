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

#include "vcf_file_structure.h"
#include "vcf_reader.h"
#include "vcf_util.h"
#include "vcf_write.h"

// #define INIT_RECORD_SIZE	100

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
 * Open a file stream and initialize the associated vcf_file_t structure.
 * 
 * @param filename The name of the file to open
 * @param max_simultaneous_batches Maximum number of batches that can be loaded in memory (0 if unlimited)
 */
vcf_file_t *vcf_open(char *filename, size_t max_simultaneous_batches);

/**
 * Close the file stream associated to a vcf_file_t type.
 * 
 * @param vcf_file vcf_file_t whose file stream is about to being closed
 */
void vcf_close(vcf_file_t *vcf_file);


/*
 * File reading
 */

// /**
//  * Fill the fields of the vcf_file_t given as argument reading data from a file.
//  * 
//  * TODO breaks because a batches list isn't provided! 
//  * a solution would be creating a temp batches list, then copying its contents to the vcf_file_t
//  * 
//  * @param vcf_file The vcf_file_t whose fields will be set
//  */
// int vcf_read(vcf_file_t *vcf_file);

int vcf_read_batches(list_t *batches_list, size_t batch_size, vcf_file_t *vcf_file);

int vcf_read_batches_in_bytes(list_t *batches_list, size_t batch_size, vcf_file_t *vcf_file);

int vcf_multiread_batches(list_t **batches_list, size_t batch_size, vcf_file_t **vcf_files, int num_files);

int vcf_parse_batches(list_t *batches_list, size_t batch_size, vcf_file_t *vcf_file, int read_samples);

int vcf_parse_batches_in_bytes(list_t *batches_list, size_t batch_size, vcf_file_t *vcf_file, int read_samples);

/**
 * Write the contents of the vcf_file_t given as argument to the given path.
 * 
 * @param vcf_file The vcf_file_t whose information will be written
 * @param filename The path of the file to write the information to
 */
int vcf_write(vcf_file_t *vcf_file, char *filename);


#endif
