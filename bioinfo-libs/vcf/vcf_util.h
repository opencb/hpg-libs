#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <log.h>

// static int debug = 0;
// static int benchmark = 1;
// 
// #define dprintf(...) { if (debug) { fprintf(stderr, __VA_ARGS__); } }
// #define bprintf(...) { if (benchmark) { fprintf(stderr, __VA_ARGS__); } }

/**
 * Flag which defines whether VCF files will be loaded using mmap.
 */
int mmap_vcf;

size_t count_regions(char *regions_string);

char *get_field_value_in_info(const char *field, char *info);

int get_field_position_in_format(const char *field, char *format);

int get_alleles(char *sample, int genotype_position, int *allele1, int *allele2);

#endif
