#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <limits.h>

#include "region.h"
#include "region_table.h"

static int debug = 0;
static int benchmark = 1;

#define dprintf(...) { if (debug) { fprintf(stderr, __VA_ARGS__); } }
#define bprintf(...) { if (benchmark) { fprintf(stderr, __VA_ARGS__); } }

size_t count_regions(char *regions_string);

region_table_t *parse_regions(char *input_regions);

#endif
