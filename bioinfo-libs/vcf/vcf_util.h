#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <sys/types.h>

static int debug = 0;
static int benchmark = 1;

#define dprintf(...) { if (debug) { fprintf(stderr, __VA_ARGS__); } }
#define bprintf(...) { if (benchmark) { fprintf(stderr, __VA_ARGS__); } }

/**
 * Flag which defines whether VCF files will be loaded using mmap.
 */
int mmap_vcf;

size_t count_regions(char *regions_string);

#endif
