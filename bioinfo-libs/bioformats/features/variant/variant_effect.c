/*
 * Copyright (c) 2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2013 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of bioinfo-libs.
 *
 * bioinfo-libs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * bioinfo-libs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with bioinfo-libs. If not, see <http://www.gnu.org/licenses/>.
 */

#include "variant_effect.h"


int invoke_effect_ws(const char *url, vcf_record_t **records, int num_records, char *excludes) {
    int variants_len = 512, current_index = 0;
    char *variants = (char*) calloc (variants_len, sizeof(char));
    const char *output_format = "txt";
    
    int chr_len, reference_len, alternate_len;
    int new_len_range;

    for (int i = 0; i < num_records; i++) {
        vcf_record_t *record = records[i];
        chr_len = record->chromosome_len;
        reference_len = record->reference_len;
        alternate_len = record->alternate_len;
        
        int num_alternates;
        char *alternates_aux = strndup(record->alternate, alternate_len);
        char **alternates = split(alternates_aux, ",", &num_alternates);

        // If a position has many alternates, each pair reference-alternate will be concatenated
        new_len_range = (current_index + chr_len + reference_len + alternate_len + 32) * num_alternates;

        // Reallocate memory if next record won't fit
        if (variants_len < (current_index + new_len_range + 1)) {
            char *aux = (char*) realloc(variants, (variants_len + new_len_range + 1) * sizeof(char));
            if (aux) { 
                variants = aux; 
                variants_len += new_len_range;
            } else {
                LOG_FATAL("Not enough memory for composing URL request\n");
            }
        }

        for (int j = 0; j < num_alternates; j++) {
            // Append region info to buffer
    //         printf("record chromosome = %.*s\n", record->chromosome_len, record->chromosome);
            strncat(variants, record->chromosome, chr_len);
            strncat(variants, ":", 1);
            current_index += chr_len + 1;
            sprintf(variants + current_index, "%lu:", record->position);
            strncat(variants, record->reference, reference_len);
            strncat(variants, ":", 1);
            strcat(variants, alternates[j]);
            strncat(variants, ",", 1);
            current_index = strlen(variants);
        }

        for (int j = 0; j < num_alternates; j++) {
            free(alternates[j]);
        }
        free(alternates);
        free(alternates_aux);
    }
    
    LOG_DEBUG_F("variants = %.*s\n", 100, variants);
//     LOG_DEBUG_F("excludes = %s\n", excludes);
    
    char *params[3] = { "of", "variants", "exclude" };
    char *params_values[3] = { output_format, variants, excludes };
    
    CURLcode ret_code = http_post(url, params, params_values, 3, save_effect_response, NULL);
    
    free(variants);
    
    return ret_code;
}

int invoke_snp_phenotype_ws(const char *url, vcf_record_t **records, int num_records) {
    CURLcode ret_code = CURLE_OK;

    const char *output_format = "txt";
    int variants_len = 512, current_index = 0;
    char *variants = (char*) calloc (variants_len, sizeof(char));
    
    int id_len, new_len_range;

    for (int i = 0; i < num_records; i++) {
        vcf_record_t *record = records[i];
        if (!strcmp(".", record->id)) {
            continue;
        }
        
        id_len = record->id_len;
        new_len_range = current_index + id_len + 32;
        
//         LOG_DEBUG_F("%s:%lu:%s:%s\n", record->chromosome, record->position, record->reference, record->alternate);
        
        // Reallocate memory if next record won't fit
        if (variants_len < (current_index + new_len_range + 1)) {
            char *aux = (char*) realloc(variants, (variants_len + new_len_range + 1) * sizeof(char));
            if (aux) { 
                variants = aux; 
                variants_len += new_len_range;
            }
        }
        
        // Append region info to buffer
        strncat(variants, record->id, id_len);
        strncat(variants, ",", 1);
        current_index += id_len + 2;
    }
    
    LOG_DEBUG_F("snps = %s\n", variants);
    
    if (current_index > 0) {
        char *params[2] = { "of", "snps" };
        char *params_values[2] = { output_format, variants };
        ret_code = http_post(url, params, params_values, 2, save_snp_phenotype_response, NULL);
    }
    
    free(variants);
    
    return ret_code;
}

int invoke_mutation_phenotype_ws(const char *url, vcf_record_t **records, int num_records) {
    CURLcode ret_code = CURLE_OK;

    const char *output_format = "txt";
    int variants_len = 512, current_index = 0;
    char *variants = (char*) calloc (variants_len, sizeof(char));
    
    int chr_len, reference_len, alternate_len;
    int new_len_range;

    for (int i = 0; i < num_records; i++) {
        vcf_record_t *record = records[i];
        if (strcmp(".", record->id)) {
            continue;
        }
        
        chr_len = record->chromosome_len;
        reference_len = record->reference_len;
        alternate_len = record->alternate_len;
        new_len_range = current_index + chr_len + reference_len + alternate_len + 32;
        
//         LOG_DEBUG_F("mutation phenotype of %s:%lu:%s:%s\n", record->chromosome, record->position, record->reference, record->alternate);
        
        // Reallocate memory if next record won't fit
        if (variants_len < (current_index + new_len_range + 1)) {
            char *aux = (char*) realloc(variants, (variants_len + new_len_range + 1) * sizeof(char));
            if (aux) { 
                variants = aux; 
                variants_len += new_len_range;
            }
        }
        
        // Append region info to buffer
        sprintf(variants + current_index, "%.*s:%lu:%.*s:%.*s,", chr_len, record->chromosome, record->position, 
                reference_len, record->reference, alternate_len, record->alternate);
        current_index = strlen(variants);
    }
    
    LOG_DEBUG_F("mutations = %s\n", variants);
    
    if (current_index > 0) {
        char *params[2] = { "of", "variants" };
        char *params_values[2] = { output_format, variants };
        ret_code = http_post(url, params, params_values, 2, save_mutation_phenotype_response, NULL);
    }
    
    free(variants);
    
    return ret_code;
}



static size_t save_effect_response(char *contents, size_t size, size_t nmemb, void *userdata) {
    int tid = omp_get_thread_num();
    
    strncat(effect_line[tid], contents, size * nmemb);
    
    char *buffer = realloc (effect_line[tid], max_line_size[tid] + size * nmemb);
    if (buffer) {
        effect_line[tid] = buffer;
        max_line_size[tid] += size * nmemb;
    } else {
        LOG_FATAL("Error while allocating memory for effect web service response");
    }
    
    return size * nmemb;
}

static size_t save_snp_phenotype_response(char *contents, size_t size, size_t nmemb, void *userdata) {
    int tid = omp_get_thread_num();
    
    strncat(snp_line[tid], contents, size * nmemb);
    
    char *buffer = realloc (snp_line[tid], snp_max_line_size[tid] + size * nmemb);
    if (buffer) {
        snp_line[tid] = buffer;
        snp_max_line_size[tid] += size * nmemb;
    } else {
        LOG_FATAL("Error while allocating memory for SNP phenotype web service response");
    }
    
    return size * nmemb;
}

static size_t save_mutation_phenotype_response(char *contents, size_t size, size_t nmemb, void *userdata) {
    int tid = omp_get_thread_num();
    
    strncat(mutation_line[tid], contents, size * nmemb);
    
    char *buffer = realloc (mutation_line[tid], mutation_max_line_size[tid] + size * nmemb);
    if (buffer) {
        mutation_line[tid] = buffer;
        mutation_max_line_size[tid] += size * nmemb;
    } else {
        LOG_FATAL("Error while allocating memory for mutation phenotype web service response");
    }
    
    return size * nmemb;
}


int initialize_ws_buffers(int num_threads) {
    // Create a buffer for each thread
    effect_line = (char**) calloc (num_threads, sizeof(char*));
    max_line_size = (int*) calloc (num_threads, sizeof(int));
    
    snp_line = (char**) calloc (num_threads, sizeof(char*));
    snp_max_line_size = (int*) calloc (num_threads, sizeof(int));
    
    mutation_line = (char**) calloc (num_threads, sizeof(char*));
    mutation_max_line_size = (int*) calloc (num_threads, sizeof(int));
    
    for (int i = 0; i < num_threads; i++) {
        max_line_size[i] = snp_max_line_size[i] = mutation_max_line_size[i] = CURL_MAX_WRITE_SIZE;
        
        effect_line[i] = (char*) calloc (max_line_size[i], sizeof(char));
        snp_line[i] = (char*) calloc (snp_max_line_size[i], sizeof(char));
        mutation_line[i] = (char*) calloc (mutation_max_line_size[i], sizeof(char));
    }
         
    return 0;
}

int free_ws_buffers(int num_threads) {
    // Free line buffers
    for (int i = 0; i < num_threads; i++) {
        free(effect_line[i]);
        free(snp_line[i]);
        free(mutation_line[i]);
    }
    
    free(max_line_size);
    free(effect_line);
        
    free(snp_max_line_size);
    free(snp_line);
        
    free(mutation_max_line_size);
    free(mutation_line);
    
    return 0;
}
