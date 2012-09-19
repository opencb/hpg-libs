#ifndef REGION_TABLE_UTILS_H
#define REGION_TABLE_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <limits.h>
#include <omp.h>

#include <bioformats/features/region/region.h>
#include <bioformats/gff/gff_file.h>
#include <bioformats/gff/gff_file_structure.h>
#include <commons/log.h>
#include <containers/region_table.h>

// region_table_t *parse_regions(char *input_regions, int as_positions);
region_table_t *parse_regions(char *input_regions, int as_positions, const char *url, const char *species, const char *version);

// region_table_t *parse_regions_from_gff_file(char *filename);
region_table_t *parse_regions_from_gff_file(char *filename, const char *url, const char *species, const char *version);

#endif
