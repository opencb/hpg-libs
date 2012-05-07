#include "vcf_util.h"

size_t count_regions(char *regions_string) {
    size_t num_regions = 0;
    char *aux = regions_string;
    while (*aux) {
        if (*aux++ == ',') ++num_regions;
    }
    return ++num_regions;
}

int get_alleles(char* sample, int* allele1, int* allele2) {
    char *aux_buffer = NULL;
    char *genotype = strtok_r(sample, "/|:", &aux_buffer);
}
