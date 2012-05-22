#include "vcf_util.h"

size_t count_regions(char *regions_string) {
    size_t num_regions = 0;
    char *aux = regions_string;
    while (*aux) {
        if (*aux++ == ',') ++num_regions;
    }
    return ++num_regions;
}

char *get_field_value_in_info(const char *field, char *info) {
    char *save_strtok, *token;
    char *value;
    
    // Search for field in info (has to begin with '<field>=')
    token = strtok_r(info, ";", &save_strtok);
    while (token != NULL && !starts_with(token, field) && strlen(token) > strlen(field)+1 && token[strlen(field)] == '=') {
        token = strtok_r(NULL, ";", &save_strtok);
    }
    
    if (token == NULL) {
        return NULL;  // Field not found
    } else {
    value = strtok_r(token, ";", &save_strtok);
    }
    
    // Search for the field value
    value = strtok_r(token, "=", &save_strtok);
    if (!strcmp(value, field)) {
//         free(value);
        value = strtok_r(NULL, "=", &save_strtok);
    }
    
    return value;
}

int get_field_position_in_format(const char *field, char *format) {
    int field_pos = 0, cur_pos = 0;
    char *save_strtok, *token;
    token = strtok_r(format, ":", &save_strtok);
    while (token != NULL && strcmp(token, field)) {
        token = strtok_r(NULL, ":", &save_strtok);
        field_pos++;
    }
    
    return (token == NULL) ? -1 : field_pos;
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
