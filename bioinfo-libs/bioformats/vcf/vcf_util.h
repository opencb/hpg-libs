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
 * @brief Getting diverse information from a VCF file
 * @details This file includes functions for getting information from a VCF file that can make 
 * easier information retrieval, such as a certain value in the FORMAT field, the value of the 
 * alleles of a sample, and so on.
 */



/**
 * Flag which defines whether VCF files will be memory-mapped instead of using the I/O API.
 * @see http://www.kernel.org/doc/man-pages/online/pages/man2/mmap.2.html
 */
int mmap_vcf;


/**
 * @brief Returns the number of comma-separated regions represented in a string
 *
 * @param regions_string The string that contains the regions
 * @return size_t The number of regions
 * 
 * Given a string containing comma-separated regions, returns the number of them. A region could 
 * be represented as chromosome:position:ref_allele:alt_allele.
 **/
size_t count_regions(char *regions_string);

/**
 * @brief Given the acronym of a field of the INFO column, retrieves its value
 *
 * @param field The acronym of the field
 * @param info The contents of the INFO column
 * @return The value associated to the acronym
 * 
 * The INFO column in a VCF file is composed of pairs (acronym, value). The structure is 
 * the following: acronym1=value1,acronym2=value2,...,acronymN=valueN. Given one of the 
 * acronyms, this function returns its associated value.
 **/
char *get_field_value_in_info(const char *field, char *info);

/**
 * @brief Given the acronym of a field in the FORMAT column, retrives its position
 *
 * @param field The acronym of the field
 * @param format The contents of the FORMAT column
 * @return The position associated to the acronym
 * 
 * The FORMAT column in a VCF file is of the form acronym1:acronym2:...:acronymN. Given one 
 * of the acronyms, this function returns its position.
 **/
int get_field_position_in_format(const char *field, char *format);

/**
 * @brief Given the value of a sample and the position of its genotype, returns the value of its alleles
 *
 * @param sample The contents of a sample column
 * @param genotype_position Position of the genotype field
 * @param allele1 [out] Value of the first allele
 * @param allele2 [out] Value of the second allele
 * @return 0 if both alleles are present, 1 if the first one is missing, 2 if the second one is missing,
 * 3 if both are missing
 * 
 * A sample in a VCF file is described in the FORMAT column, so its values are of the form 
 * value1:value2:value3. Given the position where the genotype field is, this function returns
 * the value of its alleles.
 * 
 * @see get_field_position_in_format()
 **/

int get_alleles(char *sample, int genotype_position, int *allele1, int *allele2);

#endif
