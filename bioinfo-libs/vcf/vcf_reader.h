#ifndef VCF_READER_H
#define VCF_READER_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/list.h>

#include "vcf_util.h"
#include "vcf_file_structure.h"

#define CHUNK 0x80000

extern int mmap_vcf;

/**
 * @file vcf_reader.h
 * @author Cristina Yenyxe Gonzalez Garcia
 * @brief Functions for reading VCF files
 * @details This file includes functions for reading VCF files. Files are read in blocks called 
 * "batches" whose size can be specified in lines or bytes. Files can be stored plainly or compressed 
 * in GZIP format.
 * Since reading involves also parsing contents, this tasks are provided in different schemes. The
 * vcf_read_and_parse functions allow to read and parse in just one step. If these two steps must be performed 
 * separately, the vcf_light_read functions, then run_vcf_parser, must be invoked.
 */



/**
 * @brief Current status of the VCF parser
 * 
 * The status of the parser is defined by the current header entry or record it is processing. 
 * If the current line is a record, it must also keep track of the batch where it will be inserted.
 **/
typedef struct {
    vcf_record_t *current_record;
    vcf_header_entry_t *current_header_entry;
    vcf_batch_t *current_batch;
    
    size_t num_samples;
    size_t num_records;
    
    int store_samples;
} vcf_reader_status;


/**
 * @brief Allocates memory for a vcf_reader_status structure
 * @param batch_lines The initial number of lines in a batch
 * @param store_samples Whether to parse the values of the sequenced samples 
 * @return A new vcf_reader_status structure
 * 
 * Creates a new vcf_reader_status structure. The batch where the records will be inserted during parsing 
 * will have the initial length specified as argument, although it will be resized when neccessary. The 
 * user can also decide whether to store the value of each sample or ignore them for the sake of speed and 
 * memory usage.
 **/
vcf_reader_status *vcf_reader_status_new(size_t batch_lines, int store_samples);

/**
 * @brief Frees memory associated to a vcf_reader_status structure.
 * @param status The structure to be freed
 **/

void vcf_reader_status_free(vcf_reader_status *status);


/**
 * @brief Read and parse the given number of lines from a VCF file
 *
 * @param batch_lines The number of lines to read and parse
 * @param file The file the data will be read from
 * @param read_samples Whether to parse the values of the sequenced samples
 * @return int
 **/
int vcf_read_and_parse(size_t batch_lines, vcf_file_t *file, int read_samples);

/**
 * @brief Read and parse the given number of bytes from a VCF file
 *
 * @param batch_bytes The number of bytes to read and parse
 * @param file The file the data will be read from
 * @param read_samples Whether to parse the values of the sequenced samples
 * @return int
 **/
int vcf_read_and_parse_bytes(size_t batch_bytes, vcf_file_t *file, int read_samples);

/**
 * @brief ...
 *
 * @param batch_lines The number of lines to read and parse
 * @param file The file the data will be read from
 * @param read_samples Whether to parse the values of the sequenced samples
 * @return int
 **/
int vcf_gzip_read_and_parse(size_t batch_lines, vcf_file_t *file, int read_samples);

/**
 * @brief ...
 *
 * @param batch_bytes The number of bytes to read and parse
 * @param file The file the data will be read from
 * @param read_samples Whether to parse the values of the sequenced samples
 * @return int
 **/
int vcf_gzip_read_and_parse_bytes(size_t batch_bytes, vcf_file_t *file, int read_samples);


/**
 * @brief ...
 *
 * @param text_list [out] List of blocks where the current text block will be queued
 * @param batch_lines The number of lines to read
 * @param file The file the data will be read from
 * @return int
 **/
int vcf_light_read(list_t *text_list, size_t batch_lines, vcf_file_t *file);

/**
 * @brief ...
 *
 * @param text_list [out] List of blocks where the current text block will be queued
 * @param batch_bytes The number of bytes to read
 * @param file The file the data will be read from
 * @return int
 **/
int vcf_light_read_bytes(list_t *text_list, size_t batch_bytes, vcf_file_t *file);

/**
 * @brief ...
 *
 * @param text_list [out] List of blocks where the current text block will be queued
 * @param batch_lines The number of lines to read
 * @param file The file the data will be read from
 * @return int
 **/
int vcf_gzip_light_read(list_t *text_list, size_t batch_lines, vcf_file_t *file);

/**
 * @brief ...
 *
 * @param text_list [out] List of blocks where the current text block will be queued
 * @param batch_bytes The number of bytes to read
 * @param file The file the data will be read from
 * @return int
 **/
int vcf_gzip_light_read_bytes(list_t *text_list, size_t batch_bytes, vcf_file_t *file);

/** @cond PRIVATE */
int vcf_light_multiread(list_t **batches_list, size_t batch_size, vcf_file_t **files, size_t num_files);
/** @endcond */

/**
 * 
 * @param p Pointer to the beginning of the input buffer
 * @param pe Pointer to the end of the input buffer
 * @param batch_size The number of lines of a batch, or zero when its size is specified in bytes
 * @param file VCF file whose contents will be read
 * @param status Structure that stores the current status of the parser
 */
int run_vcf_parser(char *p, char *pe, size_t batch_size, vcf_file_t *file, vcf_reader_status *status);

/** @cond PRIVATE */
size_t consume_input(int c, char **data, size_t max_len, int position_in_data);
/** @endcond */

#endif
