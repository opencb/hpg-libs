#include "vcf_filters.h"

static void annotate_failed_record(char *filter_name, size_t filter_name_len, vcf_record_t *record);

//====================================================================================
//  Filtering functions
//====================================================================================

array_list_t* coverage_filter(array_list_t* input_records, array_list_t* failed, char *filter_name, void* f_args) {
    assert(input_records);
    assert(failed);
    
    array_list_t *passed = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
    size_t filter_name_len = strlen(filter_name);

    int min_coverage = ((coverage_filter_args*)f_args)->min_coverage;

    LOG_DEBUG_F("coverage_filter (min coverage = %d) over %zu records\n", min_coverage, input_records->size);
    char *aux_buffer = (char*) calloc (128, sizeof(char));
    vcf_record_t *record;
    for (int i = 0; i < input_records->size; i++) {
        record = input_records->items[i];
        
        if (record->info_len > strlen(aux_buffer)) {
            aux_buffer = realloc (aux_buffer, record->info_len+1);
            memset(aux_buffer, 0, (record->info_len+1) * sizeof(char));
        }
        
        strncpy(aux_buffer, record->info, record->info_len);
        
        char *record_coverage = get_field_value_in_info("DP", aux_buffer);
        if (record_coverage != NULL && is_numeric(record_coverage)) {
            if (atoi(record_coverage) >= min_coverage) {
                array_list_insert(record, passed);
            } else {
                annotate_failed_record(filter_name, filter_name_len, record);
                array_list_insert(record, failed);
            }
        } else {
            annotate_failed_record(filter_name, filter_name_len, record);
            array_list_insert(record, failed);
        }
        
    }

    free(aux_buffer);
    return passed;
}

array_list_t* maf_filter(array_list_t* input_records, array_list_t* failed, char *filter_name, void* args) {
    assert(input_records);
    assert(failed);
    
    list_t *input_stats = (list_t*) malloc (sizeof(list_t));
    list_init("stats", 1, input_records->size + 1, input_stats);
    file_stats_t *file_stats = file_stats_new();
    
    array_list_t *passed = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
    size_t filter_name_len = strlen(filter_name);

    float max_maf = ((maf_filter_args*) args)->max_maf;
    float record_maf = 1.0;

    // TODO candidate for parallelization
    get_variants_stats((vcf_record_t**) input_records->items, input_records->size, input_stats, file_stats);
    
    list_item_t *stats_item = NULL;
    variant_stats_t *variant_stats;
    // The stats returned by get_variants_stats are related to a record in the same
    // position of the input_records list, so when a variant_stats_t fulfills the condition,
    // it means the related vcf_record_t passes the filter
    vcf_record_t *record;
    for (int i = 0; i < input_records->size; i++) {
        record = input_records->items[i];
        stats_item = list_remove_item(input_stats);
        variant_stats = stats_item->data_p;
        
        record_maf = 1.0;
        for (int j = 0; j < variant_stats->num_alleles; j++) {
            record_maf = fmin(record_maf, variant_stats->alleles_freq[j]);
        }
        
        if (record_maf <= max_maf) {
            array_list_insert(record, passed);
        } else {
            annotate_failed_record(filter_name, filter_name_len, record);
            array_list_insert(record, failed);
        }
        
        variant_stats_free(variant_stats);
        list_item_free(stats_item);
    }
    
    list_decr_writers(input_stats);
    file_stats_free(file_stats);
    
    return passed;
}

array_list_t* missing_values_filter(array_list_t* input_records, array_list_t* failed, char* filter_name, void* args) {
    assert(input_records);
    assert(failed);
    
    list_t *input_stats = (list_t*) malloc (sizeof(list_t));
    list_init("stats", 1, input_records->size + 1, input_stats);
    file_stats_t *file_stats = file_stats_new();
    
    array_list_t *passed = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
    size_t filter_name_len = strlen(filter_name);

    float max_missing = ((missing_values_filter_args*) args)->max_missing;
    float record_missing;
    float allele_count;

    // TODO candidate for parallelization
    get_variants_stats((vcf_record_t**) input_records->items, input_records->size, input_stats, file_stats);
    
    list_item_t *stats_item = NULL;
    variant_stats_t *variant_stats;
    // The stats returned by get_variants_stats are related to a record in the same
    // position of the input_records list, so when a variant_stats_t fulfills the condition,
    // it means the related vcf_record_t passes the filter
    vcf_record_t *record;
    for (int i = 0; i < input_records->size; i++) {
        record = input_records->items[i];
        stats_item = list_remove_item(input_stats);
        variant_stats = stats_item->data_p;
        allele_count = 0;
        
        for (int j = 0; j < variant_stats->num_alleles; j++) {
            allele_count += variant_stats->alleles_count[j];
        }
        record_missing = variant_stats->missing_alleles / (allele_count + variant_stats->missing_alleles);
        
        if (record_missing <= max_missing) {
            array_list_insert(record, passed);
        } else {
            annotate_failed_record(filter_name, filter_name_len, record);
            array_list_insert(record, failed);
        }
        
        variant_stats_free(variant_stats);
        list_item_free(stats_item);
    }
    
    list_decr_writers(input_stats);
    file_stats_free(file_stats);
    
    return passed;
}


array_list_t* num_alleles_filter(array_list_t* input_records, array_list_t* failed, char *filter_name, void* args) {
    assert(input_records);
    assert(failed);
    
    list_t *input_stats = (list_t*) malloc (sizeof(list_t));
    list_init("stats", 1, input_records->size + 1, input_stats);
    file_stats_t *file_stats = file_stats_new();
    
    array_list_t *passed = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
    size_t filter_name_len = strlen(filter_name);

    int num_alleles = ((num_alleles_filter_args*)args)->num_alleles;

    // TODO candidate for parallelization
    get_variants_stats((vcf_record_t**) input_records->items, input_records->size, input_stats, file_stats);
    
    list_item_t *stats_item = NULL;
    variant_stats_t *variant_stats;
    // The stats returned by get_variants_stats are related to a record in the same
    // position of the input_records list, so when a variant_stats_t fulfills the condition,
    // it means the related vcf_record_t passes the filter
    vcf_record_t *record;
    for (int i = 0; i < input_records->size; i++) {
        record = input_records->items[i];
        stats_item = list_remove_item(input_stats);
        variant_stats = stats_item->data_p;
        
        if (variant_stats->num_alleles == num_alleles) {
            array_list_insert(record, passed);
        } else {
            annotate_failed_record(filter_name, filter_name_len, record);
            array_list_insert(record, failed);
        }
        
        variant_stats_free(variant_stats);
        list_item_free(stats_item);
    }
    
    list_decr_writers(input_stats);
    file_stats_free(file_stats);
    
    return passed;
}


array_list_t* quality_filter(array_list_t* input_records, array_list_t* failed, char *filter_name, void* f_args) {
    assert(input_records);
    assert(failed);
    
    array_list_t *passed = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
    size_t filter_name_len = strlen(filter_name);

    int min_quality = ((quality_filter_args*)f_args)->min_quality;

    LOG_DEBUG_F("quality_filter (min quality = %d) over %zu records\n", min_quality, input_records->size);
    vcf_record_t *record;
    for (int i = 0; i < input_records->size; i++) {
        record = input_records->items[i];
        if (record->quality >= min_quality) {
            array_list_insert(record, passed);
        } else {
            annotate_failed_record(filter_name, filter_name_len, record);
            array_list_insert(record, failed);
        }
    }

    return passed;
}

array_list_t *region_filter(array_list_t *input_records, array_list_t *failed, char *filter_name, void *f_args) {
    assert(input_records);
    assert(failed);
    
    char *field;
    array_list_t *passed = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
    size_t filter_name_len = strlen(filter_name);

    region_filter_args *args = (region_filter_args*) f_args;
    region_table_t *regions = args->regions;

    LOG_DEBUG_F("region_filter over %zu records\n", input_records->size);

    int i = 0;
    vcf_record_t *record;
    region_t *region = (region_t*) malloc (sizeof(region_t));
    for (int i = 0; i < input_records->size; i++) {
        record = input_records->items[i];
        
//         LOG_DEBUG_F("record = %s, %ld\n", record->chromosome, record->position);
        
        region->chromosome = strndup(record->chromosome, record->chromosome_len);
        region->start_position = record->position;
        region->end_position = record->position;
        
        if (find_region(region, regions)) {
            // Add to the list of records that pass all checks for at least one region
            array_list_insert(record, passed);
//             LOG_DEBUG_F("%.*s, %ld passed\n", record->chromosome_len, record->chromosome, record->position);
        } else {
            // Add to the list of records that fail all checks for all regions
            annotate_failed_record(filter_name, filter_name_len, record);
            array_list_insert(record, failed);
        }
        
        free(region->chromosome);
    }

    free(region);

    return passed;
}

array_list_t *snp_filter(array_list_t *input_records, array_list_t *failed, char *filter_name, void *f_args) {
    assert(input_records);
    assert(failed);
    
    array_list_t *passed = array_list_new(input_records->size + 1, 1, COLLECTION_MODE_ASYNCHRONIZED);
    size_t filter_name_len = strlen(filter_name);

    int include_snps = ((snp_filter_args*)f_args)->include_snps;

    LOG_DEBUG_F("snp_filter (preserve SNPs = %d) over %zu records\n", include_snps, input_records->size);
    vcf_record_t *record;
    for (int i = 0; i < input_records->size; i++) {
        record = input_records->items[i];
        if (record->id_len == 1 && strncmp(".", record->id, 1) == 0) {
            if (include_snps) {
                annotate_failed_record(filter_name, filter_name_len, record);
                array_list_insert(record, failed);
            } else {
                array_list_insert(record, passed);
            }
        } else {
            if (include_snps) {
                array_list_insert(record, passed);
            } else {
                annotate_failed_record(filter_name, filter_name_len, record);
                array_list_insert(record, failed);
            }
        }
    }

    return passed;
}


//====================================================================================
//  Filter management (creation, comparison...) functions
//====================================================================================

filter_t *coverage_filter_new(int min_coverage) {
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    sprintf(filter->name, "cov%d", min_coverage);
    sprintf(filter->description, "Coverage >= %d", min_coverage);
    
    filter->type = COVERAGE;
    filter->filter_func = coverage_filter;
    filter->free_func = coverage_filter_free;
    filter->priority = 4;
    
    coverage_filter_args *filter_args = (coverage_filter_args*) malloc (sizeof(coverage_filter_args));
    filter_args->min_coverage = min_coverage;
    filter->args = filter_args;
    
    return filter;
}

void coverage_filter_free(filter_t *filter) {
    assert(filter);
    free(filter->args);
    free(filter);
}


filter_t *maf_filter_new(float max_maf) {
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    sprintf(filter->name, "maf%.2f", max_maf);
    sprintf(filter->description, "MAF <= %.2f", max_maf);
    
    filter->type = MAF;
    filter->filter_func = maf_filter;
    filter->free_func = maf_filter_free;
    filter->priority = 3;
    
    maf_filter_args *filter_args = (maf_filter_args*) malloc (sizeof(maf_filter_args));
    filter_args->max_maf = max_maf;
    filter->args = filter_args;
    
    return filter;
}

void maf_filter_free(filter_t *filter) {
    assert(filter);
    free(filter->args);
    free(filter);
}


filter_t* missing_values_filter_new(float max_missing) {
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    sprintf(filter->name, "missign%.2f", max_missing);
    sprintf(filter->description, "Missing values <= %.2f", max_missing);
    
    filter->type = MISSING_VALUES;
    filter->filter_func = missing_values_filter;
    filter->free_func = missing_values_filter_free;
    filter->priority = 5;
    
    missing_values_filter_args *filter_args = (missing_values_filter_args*) malloc (sizeof(missing_values_filter_args));
    filter_args->max_missing = max_missing;
    filter->args = filter_args;
    
    return filter;
}

void missing_values_filter_free(filter_t* filter) {
    assert(filter);
    free(filter->args);
    free(filter);
}


filter_t* num_alleles_filter_new(int num_alleles) {
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    sprintf(filter->name, "%dalleles", num_alleles);
    sprintf(filter->description, "Number of alleles = %d", num_alleles);
    
    filter->type = NUM_ALLELES;
    filter->filter_func = num_alleles_filter;
    filter->free_func = num_alleles_filter_free;
    filter->priority = 4;
    
    num_alleles_filter_args *filter_args = (num_alleles_filter_args*) malloc (sizeof(num_alleles_filter_args));
    filter_args->num_alleles = num_alleles;
    filter->args = filter_args;
    
    return filter;
}

void num_alleles_filter_free(filter_t* filter) {
    assert(filter);
    free(filter->args);
    free(filter);
}


filter_t *quality_filter_new(int min_quality) {
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    sprintf(filter->name, "q%d", min_quality);
    sprintf(filter->description, "Quality >= %d", min_quality);
    
    filter->type = QUALITY;
    filter->filter_func = quality_filter;
    filter->free_func = quality_filter_free;
    filter->priority = 4;
    
    quality_filter_args *filter_args = (quality_filter_args*) malloc (sizeof(quality_filter_args));
    filter_args->min_quality = min_quality;
    filter->args = filter_args;
    
    return filter;
}

void quality_filter_free(filter_t *filter) {
    assert(filter);
    free(filter->args);
    free(filter);
}


filter_t *region_filter_new(char *region_descriptor, int use_region_file, const char *url, const char *species, const char *version) {
    assert(region_descriptor);
    assert(url);
    assert(species);
    assert(version);
    
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    filter->type = REGION;
    filter->filter_func = region_filter;
    filter->free_func = region_filter_free;
    filter->priority = 2;

    region_filter_args *filter_args = (region_filter_args*) malloc (sizeof(region_filter_args));
    if (use_region_file) {
        snprintf(filter->name, 11, "RegionFile");
        snprintf(filter->description, 64, "Regions read from '%s'", region_descriptor);
        filter_args->regions = parse_regions_from_gff_file(region_descriptor, url, species, version);
    } else {
        snprintf(filter->name, 11, "RegionList");
        snprintf(filter->description, 64, "Regions (could be more) %s", region_descriptor);
        filter_args->regions = parse_regions(region_descriptor, 0, url, species, version);
    }
    filter->args = filter_args;

    return filter;
}

filter_t *region_exact_filter_new(char *region_descriptor, int use_region_file, const char *url, const char *species, const char *version) {
    assert(region_descriptor);
    assert(url);
    assert(species);
    assert(version);
    
    filter_t *filter = (filter_t*) malloc (sizeof(filter_t));
    filter->type = REGION;
    filter->filter_func = region_filter;
    filter->free_func = region_filter_free;
    filter->priority = 2;

    region_filter_args *filter_args = (region_filter_args*) malloc (sizeof(region_filter_args));
    if (use_region_file) {
        snprintf(filter->name, 11, "RegionFile");
        snprintf(filter->description, 64, "Regions read from '%s'", region_descriptor);
        filter_args->regions = parse_regions_from_gff_file(region_descriptor, url, species, version);
    } else {
        snprintf(filter->name, 11, "RegionList");
        snprintf(filter->description, 64, "Regions (could be more) %s", region_descriptor);
        filter_args->regions = parse_regions(region_descriptor, 1, url, species, version);
    }
    filter->args = filter_args;

    return filter;
}

void region_filter_free(filter_t *filter) {
    assert(filter);
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


filter_t *snp_filter_new(int include_snps) {
    filter_t *filter =  (filter_t*) malloc (sizeof(filter_t));
    if (include_snps) {
        snprintf(filter->name, 7, "SNPyes");
        snprintf(filter->description, 12, "To be a SNP");
    } else {
        snprintf(filter->name, 6, "SNPno");
        snprintf(filter->description, 16, "Not to be a SNP");
    }
    filter->type = SNP;
    filter->filter_func = snp_filter;
    filter->free_func = snp_filter_free;
    filter->priority = 5;

    snp_filter_args *filter_args = (snp_filter_args*) malloc (sizeof(snp_filter_args));
    filter_args->include_snps = include_snps;
    filter->args = filter_args;

    return filter;
}

void snp_filter_free(filter_t *filter) {
    assert(filter);
    free(filter->args);
    free(filter);
}


int filter_compare(const void *filter1, const void *filter2) {
    assert(filter1);
    assert(filter2);
    return ((filter_t*) filter1)->priority - ((filter_t*) filter2)->priority;
}



//====================================================================================
//  Filter chain functions
//====================================================================================

filter_chain *add_to_filter_chain(filter_t *filter, filter_chain *chain) {
    assert(filter);
    
    filter_chain *result = chain;
    if (result == NULL) {
        result = cp_heap_create((cp_compare_fn) filter_compare);
    }
    assert(result);
    cp_heap_push(result, filter);

    return result;
}

filter_t **sort_filter_chain(filter_chain *chain, int *num_filters) {
    assert(chain);
    assert(num_filters);
    
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

array_list_t *run_filter_chain(array_list_t *input_records, array_list_t *failed, filter_t **filters, int num_filters) {
    assert(input_records);
    assert(failed);
    assert(filters);
    
    array_list_t *passed = input_records;
    array_list_t *aux_passed;

    LOG_DEBUG_F("Applying filter chain of %d filters\n", num_filters);

    // Apply each filter with the arguments provided
    for (int i = 0; i < num_filters; i++) {
        filter_t *filter = filters[i];
        aux_passed = filter->filter_func(passed, failed, filter->name, filter->args);
    // 		free(passed);
        passed = aux_passed;
    }

    // Mark records that passed all filters
    vcf_record_t *record;
    for (int i = 0; i < passed->size; i++) {
        record = passed->items[i];
        // If the filter didn't fail filters from another run
        if (!strncmp(record->filter, ".", 1)) {
            set_vcf_record_filter("PASS", 4, record);
        }
    }
    
    return passed;
}

void free_filter_chain(filter_chain* chain) {
    assert(chain);
    filter_t *filter;
    while ((filter = cp_heap_pop(chain) != NULL)) {
        filter->free_func(filter);
    }
    cp_heap_destroy(chain);
}

void free_filters(filter_t **filters, int num_filters) {
    assert(filters);
    filter_t *filter;
    for (int i = 0; i < num_filters; i++) {
        filter = filters[i];
        filter->free_func(filter);
    }
    free(filters);
}



//====================================================================================
//  Other functions
//====================================================================================

vcf_header_entry_t **get_filters_as_vcf_headers(filter_t **filters, int num_filters) {
    vcf_header_entry_t **headers = malloc (num_filters * sizeof(vcf_header_entry_t*));
    vcf_header_entry_t *entry;
    char *value;
    size_t value_len;
    
    for (int i = 0; i < num_filters; i++) {
        entry = vcf_header_entry_new();
        entry->name = strndup("FILTER", 6);
        entry->name_len = 6;
        
        value_len = strlen(filters[i]->name) + 4;
        value = malloc(value_len * sizeof(char));
        sprintf(value, "ID=%s", filters[i]->name);
        value[value_len-1] = '\0';
        LOG_DEBUG_F("ID of the defined filter = %s\n", value);
        add_vcf_header_entry_value(value, value_len-1, entry);
        free(value);
        
        value_len = strlen(filters[i]->description) + 15;
        value = malloc(value_len * sizeof(char));
        sprintf(value, "Description=\"%s\"", filters[i]->description);
        value[value_len-1] = '\0';
        LOG_DEBUG_F("Description of the defined filter = %s\n", value);
        add_vcf_header_entry_value(value, value_len-1, entry);
        free(value);
        
        headers[i] = entry;
    }
    
    return headers;
}

static void annotate_failed_record(char *filter_name, size_t filter_name_len, vcf_record_t *record) {
    if (!strncmp(record->filter, ".", record->filter_len) || !strncmp(record->filter, "PASS", record->filter_len)) {
        record->filter = filter_name;
        record->filter_len = filter_name_len;
    } else {
        char *aux = calloc(record->filter_len + filter_name_len + 2, sizeof(char));
        if (aux) {
            strncat(aux, record->filter, record->filter_len);
            strncat(aux, ",", 1);
            strncat(aux, filter_name, filter_name_len);
            record->filter = aux;
            record->filter_len = record->filter_len + filter_name_len + 2;
        } else {
            LOG_FATAL_F("Memory allocation error while annotating record %.*%s:%ld\n", 
                        record->chromosome_len, record->chromosome, record->position);
        }
    }
}
