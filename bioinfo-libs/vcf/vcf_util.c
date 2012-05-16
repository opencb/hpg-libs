#include "vcf_util.h"

size_t count_regions(char *regions_string) {
    size_t num_regions = 0;
    char *aux = regions_string;
    while (*aux) {
        if (*aux++ == ',') ++num_regions;
    }
    return ++num_regions;
}

int get_genotype_position_in_format(char *format) {
    int gt_pos = 0, cur_pos = 0;
    char *save_strtok, *token;
    token = strtok_r(format, ":", &save_strtok);
    while (token != NULL && strcmp(token, "GT")) {
        token = strtok_r(NULL, ":", &save_strtok);
        gt_pos++;
    }
    
    return (token == NULL) ? -1 : gt_pos;
}

int get_alleles(char* sample, int genotype_position, int* allele1, int* allele2) {
    char *aux_buffer, *allele, *genotype;
    int ret_code = 0, cur_pos = -1;
    while (cur_pos < genotype_position) {
        genotype = strtok_r(sample, ":", &aux_buffer);
        sample = NULL;
        cur_pos++;
    }
    
    LOG_DEBUG_F("genotype = %s\n", genotype);
    
    allele = strtok_r(genotype, "/|", &aux_buffer);
    if (strcmp(".", allele) == 0) {
        *allele1 = -1;
        ret_code += 1;
    } else {
        *allele1 = atoi(allele);
    }
    
    LOG_DEBUG_F("allele 1 = %s\n", allele);
    
    allele = strtok_r(NULL, "/|", &aux_buffer);
    if (strcmp(".", allele) == 0) {
        *allele2 = -1;
        ret_code += 2;
    } else {
        *allele2 = atoi(allele);
    }
    
    LOG_DEBUG_F("allele2 = %s\n", allele);
    
    return ret_code;
}
