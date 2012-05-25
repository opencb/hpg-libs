#include "vcf_filters.h"


filter_t *create_coverage_filter(int min_coverage) {
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    filter->type = COVERAGE;
    filter->filter_func = coverage_filter;
    filter->free_func = free_coverage_filter;
    filter->priority = 4;
    
    coverage_filter_args *filter_args = (coverage_filter_args*) malloc (sizeof(coverage_filter_args));
    filter_args->min_coverage = min_coverage;
    filter->args = filter_args;
    
    return filter;
}

void free_coverage_filter(filter_t *filter) {
    free(filter->args);
    free(filter);
}

filter_t* create_num_alleles_filter(int num_alleles) {
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    filter->type = NUM_ALLELES;
    filter->filter_func = num_alleles_filter;
    filter->free_func = free_num_alleles_filter;
    filter->priority = 4;
    
    num_alleles_filter_args *filter_args = (num_alleles_filter_args*) malloc (sizeof(num_alleles_filter_args));
    filter_args->num_alleles = num_alleles;
    filter->args = filter_args;
    
    return filter;
}

void free_num_alleles_filter(filter_t* filter) {
    free(filter->args);
    free(filter);
}


filter_t *create_quality_filter(int min_quality) {
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    filter->type = QUALITY;
    filter->filter_func = quality_filter;
    filter->free_func = free_quality_filter;
    filter->priority = 4;
    
    quality_filter_args *filter_args = (quality_filter_args*) malloc (sizeof(quality_filter_args));
    filter_args->min_quality = min_quality;
    filter->args = filter_args;
    
    return filter;
}

void free_quality_filter(filter_t *filter) {
    free(filter->args);
    free(filter);
}

filter_t *create_region_filter(char *region_descriptor, int use_region_file) {
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    filter->type = REGION;
    filter->filter_func = region_filter;
    filter->free_func = free_region_filter;
    filter->priority = 2;

    region_filter_args *filter_args = (region_filter_args*) malloc (sizeof(region_filter_args));
    if (use_region_file) {
        filter_args->regions = parse_regions_from_gff_file(region_descriptor);
    } else {
        filter_args->regions = parse_regions(region_descriptor, 0);
    }
    filter->args = filter_args;

    return filter;
}

filter_t *create_region_exact_filter(char *region_descriptor, int use_region_file) {
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    filter->type = REGION;
    filter->filter_func = region_filter;
    filter->free_func = free_region_filter;
    filter->priority = 2;

    region_filter_args *filter_args = (region_filter_args*) malloc (sizeof(region_filter_args));
    if (use_region_file) {
        filter_args->regions = parse_regions_from_gff_file(region_descriptor);
    } else {
        filter_args->regions = parse_regions(region_descriptor, 1);
    }
    filter->args = filter_args;

    return filter;
}

void free_region_filter(filter_t *filter) {
    region_table_t *regions = ((region_filter_args*) filter->args)->regions;
    // Free ordering array
    char **ordering = regions->ordering;
    for (int i = 0; i < regions->max_chromosomes; i++) {
        free(ordering[i]);
    }
    free(ordering);

    // Free hashtable
    cp_hashtable *table = regions->storage;
    cp_hashtable_destroy(table);

    // Free pointers to args and to the filter itself
    free(regions);
    free(filter->args);
    free(filter);
}


int filter_compare(const void *filter1, const void *filter2) {
    return ((filter_t*) filter1)->priority - ((filter_t*) filter2)->priority;
}

filter_t *create_snp_filter(char *include_snps) {
    filter_t *filter =  (filter_t*) malloc (sizeof(filter_t));
    filter->type = SNP;
    filter->filter_func = snp_filter;
    filter->free_func = free_snp_filter;
    filter->priority = 5;

    snp_filter_args *filter_args = (snp_filter_args*) malloc (sizeof(snp_filter_args));
    filter_args->include_snps = 1;	// Default: Include SNPs

    if (include_snps != NULL) {
        if (strcmp("include", include_snps) == 0) {
            filter_args->include_snps = 1;
        } else if (strcmp("exclude", include_snps) == 0) {
            filter_args->include_snps = 0;
        }
    }

    filter->args = filter_args;

    return filter;
}

void free_snp_filter(filter_t *filter) {
    free(filter->args);
    free(filter);
}



filter_chain *add_to_filter_chain(filter_t *filter, filter_chain *chain) {
    filter_chain *result = chain;

    if (result == NULL) {
        result = cp_heap_create((cp_compare_fn) filter_compare);
    }
    cp_heap_push(result, filter);

    return result;
}

filter_t **sort_filter_chain(filter_chain *chain, int *num_filters) {
    *num_filters = cp_heap_count(chain);
    filter_t **filters = (filter_t**) malloc (cp_heap_count(chain) * sizeof(filter_t*));

    // Pop all filters from the heap and make an ordered list from them
    filter_t *filter = NULL;
    for (int i = 0; (filter = cp_heap_pop(chain)) != NULL; i++) {
        LOG_DEBUG_F("Filter %d type = %d\n", i, filter->type);
        filters[i] = filter;
    }

    return filters;
}

list_t *run_filter_chain(list_t *input_records, list_t *failed, filter_t **filters, int num_filters) {
    list_t *passed = input_records;
    list_t *aux_passed;

    LOG_DEBUG_F("Applying filter chain of %d filters\n", num_filters);

    // Apply each filter with the arguments provided
    for (int i = 0; i < num_filters; i++) {
        filter_t *filter = filters[i];
        aux_passed = filter->filter_func(passed, failed, filter->args);
    // 		free(passed);
        passed = aux_passed;
    }

    return passed;
}

void free_filter_chain(filter_chain* chain) {
    if (chain) { cp_heap_destroy(chain); }
}



list_t* coverage_filter(list_t* input_records, list_t* failed, void* f_args) {
    list_t *passed = (list_t*) malloc (sizeof(list_t));
    list_init("passed", 1, input_records->max_length, passed);

    int min_coverage = ((coverage_filter_args*)f_args)->min_coverage;

    LOG_DEBUG_F("coverage_filter (min coverage = %d) over %zu records\n", min_coverage, input_records->length);
    vcf_record_t *record;
    char *aux_buffer = (char*) calloc (128, sizeof(char));
    for (list_item_t *item = input_records->first_p; item != NULL; item = item->next_p) {
        record = item->data_p;
        list_item_t *new_item = list_item_new(item->id, item->type, record);
        
        if (strlen(record->info) > strlen(aux_buffer)) {
            aux_buffer = realloc (aux_buffer, strlen(record->info)+1);
            memset(aux_buffer, 0, (strlen(record->info)+1) * sizeof(char));
        }
        
        strncpy(aux_buffer, record->info, strlen(record->info));
        
        char *record_coverage = get_field_value_in_info("DP", aux_buffer);
        if (record_coverage != NULL && is_numeric(record_coverage)) {
            if (atoi(record_coverage) >= min_coverage) {
                list_insert_item(new_item, passed);
            } else {
                list_insert_item(new_item, failed);
            }
        } else {
            list_insert_item(new_item, failed);
        }
        
    }

    free(aux_buffer);
    return passed;
}

list_t* num_alleles_filter(list_t* input_records, list_t* failed, void* args) {
    list_t *input_stats = (list_t*) malloc (sizeof(list_t));
    list_init("stats", 1, input_records->max_length, input_stats);
    file_stats_t *file_stats = new_file_stats();
    
    list_t *passed = (list_t*) malloc (sizeof(list_t));
    list_init("passed", 1, input_records->max_length, passed);

    int num_alleles = ((num_alleles_filter_args*)args)->num_alleles;

    // TODO candidate for parallelization
    get_variants_stats(input_records->first_p, input_records->length, input_stats, file_stats);
    
    list_item_t *stats_item = NULL;
    list_item_t *record_item = input_records->first_p;
    variant_stats_t *variant_stats;
    // The stats returned by get_variants_stats are related to the records in the same
    // position of the input_records list, so when a variant_stats_t fulfills the condition,
    // it means the related vcf_record_t passes the filter
    while (record_item != NULL) {
        stats_item = list_remove_item(input_stats);
        variant_stats = stats_item->data_p;
        
        list_item_t *new_item = list_item_new(record_item->id, record_item->type, record_item->data_p);
        if (variant_stats->num_alleles == num_alleles) {
            list_insert_item(new_item, passed);
        } else {
            list_insert_item(new_item, failed);
        }
        
        free_variant_stats(variant_stats);
        list_item_free(stats_item);
        
        record_item = record_item->next_p;
    }
    
    list_decr_writers(input_stats);
    
    return passed;
}


list_t* quality_filter(list_t* input_records, list_t* failed, void* f_args) {
    list_t *passed = (list_t*) malloc (sizeof(list_t));
    list_init("passed", 1, input_records->max_length, passed);

    int min_quality = ((quality_filter_args*)f_args)->min_quality;

    LOG_DEBUG_F("quality_filter (min quality = %d) over %zu records\n", min_quality, input_records->length);
    vcf_record_t *record;
    for (list_item_t *item = input_records->first_p; item != NULL; item = item->next_p) {
        record = item->data_p;
        list_item_t *new_item = list_item_new(item->id, item->type, record);
        if (record->quality >= min_quality) {
            list_insert_item(new_item, passed);
        } else {
            list_insert_item(new_item, failed);
        }
    }

    return passed;
}

list_t *region_filter(list_t *input_records, list_t *failed, void *f_args) {
    char *field;
    list_t *passed = (list_t*) malloc (sizeof(list_t));
    list_init("passed", 1, INT_MAX, passed);

    region_filter_args *args = (region_filter_args*) f_args;
    region_table_t *regions = args->regions;

    LOG_DEBUG_F("region_filter over %zu records\n", input_records->length);

    int i = 0;
    region_t *region = (region_t*) malloc (sizeof(region_t));
    for (list_item_t *item = input_records->first_p; item != NULL; item = item->next_p) {
        vcf_record_t *record = item->data_p;
        list_item_t *new_item = list_item_new(item->id, item->type, record);
        
        LOG_DEBUG_F("record = %s, %ld\n", record->chromosome, record->position);
        
        region->chromosome = record->chromosome;
        region->start_position = record->position;
        region->end_position = record->position;
        
        if (find_region(region, regions)) {
            // Add to the list of records that pass all checks for at least one region
            list_insert_item(new_item, passed);
            LOG_DEBUG_F("%s, %ld passed\n", record->chromosome, record->position);
        } else {
            // Add to the list of records that fail all checks for all regions
            list_insert_item(new_item, failed);
        }
        
        i++;
    }

    free(region);

    return passed;
}

list_t *snp_filter(list_t *input_records, list_t *failed, void *f_args) {
    list_t *passed = (list_t*) malloc (sizeof(list_t));
    list_init("passed", 1, input_records->max_length, passed);

    int include_snps = ((snp_filter_args*)f_args)->include_snps;

    LOG_DEBUG_F("snp_filter (preserve SNPs = %d) over %zu records\n", include_snps, input_records->length);
    for (list_item_t *item = input_records->first_p; item != NULL; item = item->next_p) {
        vcf_record_t *record = item->data_p;
        list_item_t *new_item = list_item_new(item->id, item->type, record);
        // TODO check 'id' field is not empty (modifications to the parser needed!)
        if (strcmp(".", record->id) == 0) {
            if (include_snps) {
                list_insert_item(new_item, failed);
            } else {
                list_insert_item(new_item, passed);
            }
        } else {
            if (include_snps) {
                list_insert_item(new_item, passed);
            } else {
                list_insert_item(new_item, failed);
            }
        }
    }

    return passed;
}
