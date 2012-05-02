#include "vcf_util.h"

size_t count_regions(char *regions_string)
{
    size_t num_regions = 0;
    char *aux = regions_string;
    while (*aux) {
        if (*aux++ == ',') ++num_regions;
    }
    return ++num_regions;
}
