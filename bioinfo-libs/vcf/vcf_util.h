#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <commons/log.h>


/**
 * @file vcf_util.h
 * @author Cristina Yenyxe Gonzalez Garcia
 * @brief Functions for getting diverse information from a VCF file
 * @details This file includes functions for getting information from a VCF file that can make 
 * easier information retrieval, such as a certain value in the FORMAT field, the value of the 
 * alleles of a sample, and so on.
 */



/**
 * Flag which defines whether VCF files will be memory-mapped instead of using the I/O API.
 * @see http://www.kernel.org/doc/man-pages/online/pages/man2/mmap.2.html
 */
int mmap_vcf;

size_t count_regions(char *regions_string);

char *get_field_value_in_info(const char *field, char *info);

int get_field_position_in_format(const char *field, char *format);

int get_alleles(char *sample, int genotype_position, int *allele1, int *allele2);

#endif
