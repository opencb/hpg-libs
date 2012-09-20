#ifndef VCF_FILTERS_H 
#define VCF_FILTERS_H

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <cprops/heap.h>

#include <containers/array_list.h>
#include <containers/list.h>
#include <containers/region_table_utils.h>
#include <commons/log.h>

#include "vcf_file_structure.h"
#include "vcf_stats.h"
#include "vcf_util.h"

/**
 * @file vcf_filters.h
 * @author Cristina Yenyxe Gonzalez Garcia
 * @brief Filters for VCF files
 * @details This file includes functions for filtering VCF files. All filters receive as input a list of 
 * records, and return two lists: one with the records that passed the filters and another with the ones 
 * that were rejected.
 */



/**
 * @brief The type of the filter to apply
 **/
enum filter_type { COVERAGE, NUM_ALLELES, QUALITY, REGION, SNP  };

/**
 * @brief Arguments for the filter by coverage
 * 
 * The only argument of a filter by coverage is the minimum coverage of a record, as specified on its 
 * INFO field.
 **/
typedef struct {
    int min_coverage;
} coverage_filter_args;

/**
 * @brief Arguments for the filter by number of alleles
 * 
 * The only argument of a filter by number of alleles is precisely that number.
 **/
typedef struct {
    int num_alleles;
} num_alleles_filter_args;

/**
 * @brief Arguments for the filter by quality
 * 
 * The only argument of a filter by quality is the minimum quality of a record, as specified on its 
 * QUAL field.
 **/
typedef struct {
    int min_quality;
} quality_filter_args;

/**
 * @brief Arguments for the filter by region
 * 
 * The argument of a filter by region is a set of one or more regions of the form 
 * chromosome:position:ref_allele:alt_allele.
 **/
typedef struct {
    region_table_t *regions;
} region_filter_args;

/**
 * @brief Arguments for the filter by SNP
 * 
 * The only argument of a filter by SNP specifies whether to include (1) or exclude (0) a SNP.
 **/
typedef struct {
    int include_snps;	// 1 = preserve SNPs, 0 = remove SNPs
} snp_filter_args;

/**
 * @brief A filter selects a subcollection of records which fulfill some condition.
 * 
 * A filter selects a subcollection of records which fulfill some condition.
 * It is mandatory to provide the list of records to filter and a list to store in 
 * the records that failed the filter's test.
 * 
 * If more than one filter is applied, they must be ordered by priority (max = 0).
 */
typedef struct filter {
	unsigned int priority;  /**< Sorting criteria when several filters are applied */
	enum filter_type type;  /**< Filtering criteria */
	array_list_t* (*filter_func) (array_list_t *input_records, array_list_t *failed, void *args);  /**< Filtering function itself */
	void (*free_func) (struct filter *f);   /**< Filter deallocation function */
	void *args;             /**< Filter-dependant arguments */
} filter_t;

typedef cp_heap filter_chain;


//====================================================================================
//  Filtering functions prototypes
//====================================================================================

/**
 * Given a list of records, check which ones have a coverage greater or equals than 
 * the one specified.
 * 
 * @param records List of records to filter
 * @param failed Records that failed the filter's test
 * @return Records that passed the filter's test
 */
array_list_t *coverage_filter(array_list_t *input_records, array_list_t *failed, void *args);

/**
 * Given a list of records, check which ones have a num_alleles equals to 
 * the one specified.
 * 
 * @param records List of records to filter
 * @param failed Records that failed the filter's test
 * @return Records that passed the filter's test
 */
array_list_t *num_alleles_filter(array_list_t *input_records, array_list_t *failed, void *args);

/**
 * Given a list of records, check which ones have a quality greater or equals than 
 * the one specified.
 * 
 * @param records List of records to filter
 * @param failed Records that failed the filter's test
 * 
 * @return Records that passed the filter's test
 */
array_list_t *quality_filter(array_list_t *input_records, array_list_t *failed, void *args);

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
array_list_t *region_filter(array_list_t *input_records, array_list_t *failed, void *args);

/**
 * Given a list of records, check which ones represent a SNP.
 * 
 * @param records List of records to filter
 * @param failed Records that failed the filter's test
 * 
 * @return Records that passed the filter's test
 */
array_list_t *snp_filter(array_list_t *input_records, array_list_t *failed, void *args);


//====================================================================================
//  Filter management (creation, comparison...) functions prototypes
//====================================================================================

/**
 * @brief Creates a new filter by minimum coverage
 *
 * @param min_coverage Minimum coverage for the records to pass the filter
 * @return The new filter
 **/
filter_t *coverage_filter_new(int min_coverage);

/**
 * @brief Deallocates memory of a filter by minimum coverage
 *
 * @param filter The filter to deallocate
 **/
void coverage_filter_free(filter_t *filter);

/**
 * @brief Creates a new filter by number of alleles
 *
 * @param num_alleles Number of alleles of the records that pass the filter
 * @return The new filter
 **/
filter_t *num_alleles_filter_new(int num_alleles);

/**
 * @brief Deallocates memory of a filter by number of alleles
 *
 * @param filter The filter to deallocate
 **/
void num_alleles_filter_free(filter_t *filter);

/**
 * @brief Creates a new filter by minimum quality
 *
 * @param min_quality Minimum quality for the records to pass the filter
 * @return The new filter
 **/
filter_t *quality_filter_new(int min_quality);

/**
 * @brief Deallocates memory of a filter by minimum quality
 *
 * @param filter The filter to deallocate
 **/
void quality_filter_free(filter_t *filter);

/**
 * @brief ...
 *
 * @param region_descriptor ...
 * @param use_region_file ...
 * @param url ...
 * @param species ...
 * @param version ...
 * @return The new filter
 **/
filter_t *region_filter_new(char *region_descriptor, int use_region_file, const char *url, const char *species, const char *version);

/**
 * @brief ...
 *
 * @param region_descriptor ...
 * @param use_region_file ...
 * @param url ...
 * @param species ...
 * @param version ...
 * @return The new filter
 **/
filter_t *region_exact_filter_new(char *region_descriptor, int use_region_file, const char *url, const char *species, const char *version);

/**
 * @brief Deallocates memory of a filter by region
 *
 * @param filter The filter to deallocate
 **/
void region_filter_free(filter_t *filter);

/**
 * @brief Creates a new filter by SNP
 *
 * @param include_snps Whether to include or exclude a SNP.
 * @return The new filter
 **/
filter_t *snp_filter_new(char *include_snps);

/**
 * @brief Deallocates memory of a filter by SNP
 *
 * @param filter The filter to deallocate
 **/
void snp_filter_free(filter_t *filter);


/**
 * @brief Compares the priority of two filters
 *
 * @param filter1 First filter to compare
 * @param filter2 Second filter to compare
 * @return 0 if both filters have the same priority, less than zero if the 1st filter 
 * has less priority, and more than zero if it has more priority
 **/
int filter_compare(const void *filter1, const void *filter2);


//====================================================================================
//  Filter chain functions prototypes
//====================================================================================

/**
 * @brief Adds a filter to the given filter chain
 * @param filter Filter to add to the filter chain
 * @param chain Filter chain the filter is inserted in
 * @return The new state of the filter chain
 * 
 * Adds a filter to the given filter chain. If the chain is NULL, the filter is added 
 * after creating that chain.
 */
filter_chain *add_to_filter_chain(filter_t *filter, filter_chain *chain);

/**
 * @brief Sorts a chain of several filters by priority
 * @param chain Filter chain to sort
 * @return Sorted list of filters
 * 
 * Given a chain of several filters, creates a list sorted by priority.
 */
filter_t **sort_filter_chain(filter_chain *chain, int *num_filters);

/**
 * @brief Frees memory allocated to store a filter chain
 * @param chain Chain of filters to free
 */
void free_filter_chain(filter_chain *chain);

/**
 * Applies a collection of filters to a list of records.
 * 
 * @param input_records List of records to filter
 * @param failed Records that failed the filter's test
 * @param filters Filters to apply
 * @param num_filters Number of filters to apply
 * @return Records that passed the filters' tests
 */
array_list_t *run_filter_chain(array_list_t *input_records, array_list_t *failed, filter_t **filters, int num_filters);

#endif

