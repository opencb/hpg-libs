#ifndef VCF_FILTERS_H 
#define VCF_FILTERS_H

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <cprops/heap.h>

#include <list.h>
#include <log.h>
#include <region_table_utils.h>

#include "vcf_file_structure.h"
#include "vcf_stats.h"
#include "vcf_util.h"

//====================================================================================
//  vcf_filter.h
//
//  vcf_filter structures and prototypes
//====================================================================================

enum filter_type { COVERAGE, NUM_ALLELES, QUALITY, REGION, SNP  };

typedef struct {
    int min_coverage;
} coverage_filter_args;

typedef struct {
    int num_alleles;
} num_alleles_filter_args;

typedef struct {
    int min_quality;
} quality_filter_args;

typedef struct {
    region_table_t *regions;
} region_filter_args;

typedef struct {
    int include_snps;	// 1 = preserve SNPs, 0 = remove SNPs
} snp_filter_args;

/**
 * A filter selects a subcollection of records which fulfill some condition.
 * It is mandatory to provide the list of records to filter and a list to store in 
 * the records that failed the filter's test.
 * 
 * The args variable stores filter-dependant arguments.
 * 
 * If more than one filter are applied, they must be ordered by priority (max = 0).
 */
typedef struct filter {
	unsigned int priority;
	enum filter_type type;
	list_t* (*filter_func) (list_t *input_records, list_t *failed, void *args);
	void (*free_func) (struct filter *f);
	void *args;
} filter_t;

typedef cp_heap filter_chain;


//====================================================================================
//  Filter management (creation, comparison...) functions prototypes
//====================================================================================

filter_t *create_coverage_filter(int min_coverage);

void free_coverage_filter(filter_t *filter);

filter_t *create_coverage_filter(int min_coverage);

void free_coverage_filter(filter_t *filter);

filter_t *create_num_alleles_filter(int num_alleles);

void free_num_alleles_filter(filter_t *filter);

filter_t *create_quality_filter(int min_quality);

void free_quality_filter(filter_t *filter);

filter_t *create_region_filter(char *region_descriptor, int use_region_file, const char *url, const char *species, const char *version);

filter_t *create_region_exact_filter(char *region_descriptor, int use_region_file, const char *url, const char *species, const char *version);

void free_region_filter(filter_t *filter);

filter_t *create_snp_filter(char *include_snps);

void free_snp_filter(filter_t *filter);


int filter_compare(const void *filter1, const void *filter2);


//====================================================================================
//  Filter chain functions prototypes
//====================================================================================

// filter_chain *create_filter_chain(int writers);

/**
 * Add a filter to the given filter chain. If the chain is NULL, the filter is added 
 * after creating that chain.
 * 
 * @param filter
 * 	Filter to add to the filter chain
 * @param chain
 * 	Filter chain the filter is inserted in
 * 
 * @return
 * 	The new state of the filter chain
 */
filter_chain *add_to_filter_chain(filter_t *filter, filter_chain *chain);

/**
 * Given a chain of several filters, creates a list sorted by priority.
 * 
 * @param chain
 * 	Filter chain to order
 * 
 * @return
 * 	Sorted list of filters
 */
filter_t **sort_filter_chain(filter_chain *chain, int *num_filters);

/**
 * Free memory allocated to store a filter chain.
 * 
 * @param chain
 * Chain of filters to apply
 */
void free_filter_chain(filter_chain *chain);

/**
 * Applies a collection of filters to a list of records.
 * 
 * @param input_records
 * 	List of records to filter
 * @param failed
 * 	Records that failed the filter's test
 * @param filters
 * 	Filters to apply
 * @param num_filters
 * 	Number of filters to apply
 * 
 * @return Records that passed the filters' tests
 */
list_t *run_filter_chain(list_t *input_records, list_t *failed, filter_t **filters, int num_filters);


//====================================================================================
//  Filtering functions prototypes
//====================================================================================

/**
 * Given a list of records, check which ones have a coverage greater or equals than 
 * the one specified.
 * 
 * @param records List of records to filter
 * @param failed Records that failed the filter's test
 * 
 * @return Records that passed the filter's test
 */
list_t *coverage_filter(list_t *input_records, list_t *failed, void *args);

/**
 * Given a list of records, check which ones have a num_alleles equals to 
 * the one specified.
 * 
 * @param records List of records to filter
 * @param failed Records that failed the filter's test
 * 
 * @return Records that passed the filter's test
 */
list_t *num_alleles_filter(list_t *input_records, list_t *failed, void *args);

/**
 * Given a list of records, check which ones have a quality greater or equals than 
 * the one specified.
 * 
 * @param records List of records to filter
 * @param failed Records that failed the filter's test
 * 
 * @return Records that passed the filter's test
 */
list_t *quality_filter(list_t *input_records, list_t *failed, void *args);

/**
 * Given a list of records, check which ones are positioned in certain genome region.
 * A region is defined by a pair of fields: the chromosome and a position or range 
 * of positions.
 * 
 * @param input_records List of records to filter
 * @param failed Records that failed the filter's test
 * @param regions List of regions where the records must be placed
 * @param num_regions Number of regions to check
 * 
 * @return Records that passed the filter's test
 */
list_t *region_filter(list_t *input_records, list_t *failed, void *args);

/**
 * Given a list of records, check which ones represent a SNP.
 * 
 * @param records List of records to filter
 * @param failed Records that failed the filter's test
 * 
 * @return Records that passed the filter's test
 */
list_t *snp_filter(list_t *input_records, list_t *failed, void *args);


#endif

