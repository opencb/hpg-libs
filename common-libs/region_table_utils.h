#ifndef REGION_TABLE_UTILS_H
#define REGION_TABLE_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <limits.h>
#include <omp.h>

#include <gff_file.h>
#include <gff_file_structure.h>
#include <region.h>
#include <region_table.h>

region_table_t *parse_regions(char *input_regions);

region_table_t *parse_regions_from_gff_file(char *filename);

#endif
