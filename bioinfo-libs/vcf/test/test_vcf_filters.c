#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "vcf_batch.h"
#include "vcf_file.h"
#include "vcf_file_structure.h"
#include "vcf_filters.h"

#include "list.h" 

/*
 * How to get the constants defined below
 * cat CEU.exon.2010_03.genotypes__head400.vcf | grep "#" | wc -l
 * Result: 11 lines, so there are 389 records
 * tail -n 389 CEU.exon.2010_03.genotypes__head400.vcf | cut -f 3 | grep "rs" | wc -l
 * tail -n 389 CEU.exon.2010_03.genotypes__head400.vcf | cut -f 3 | grep "\." | wc -l
 */
#define MAX_RECORDS	    389
#define SNPS_IN_FILE	266


Suite *create_test_suite();
list_t *read_test_datasuite(vcf_file_t *file);
void free_test_datasuite(list_t *datasuite, vcf_file_t *file);


list_t *datasuite, *quality_datasuite;
list_t *passed, *failed;

filter_t *coverage_f, *quality_f, *region_f, *snp_f;
filter_chain *chain;



/* ******************************
 *       Unchecked fixtures     *
 * ******************************/

void setup_snp(void)
{
	printf("Begin SNP filter testing\n");
	snp_f = create_snp_filter(NULL);
}

void teardown_snp(void)
{
	printf("Finished SNP filter testing\n");
	snp_f->free_func(snp_f);
}

void setup_region(void)
{
	printf("Begin region filter testing\n");
}

void teardown_region(void)
{
	printf("Finished region filter testing\n");
}

void setup_quality(void)
{
    printf("Begin quality filter testing\n");
    quality_f = create_quality_filter(30);
}

void teardown_quality(void)
{
    printf("Finished quality filter testing\n");
    quality_f->free_func(quality_f);
}

void setup_coverage(void)
{
    printf("Begin coverage filter testing\n");
    coverage_f = create_coverage_filter(50);
}

void teardown_coverage(void)
{
    printf("Finished coverage filter testing\n");
    coverage_f->free_func(coverage_f);
}

void setup_snp_region(void)
{
	printf("Begin SNP + region filter chain testing\n");
	
	snp_f = create_snp_filter(NULL);
	char input[] = "1:1000000-5000000";
	region_f = create_region_filter(input, 0);
	
	chain = add_to_filter_chain(region_f, chain);
	chain = add_to_filter_chain(snp_f, chain);
}

void teardown_snp_region(void)
{
	printf("Finished SNP + region filter chain testing\n");
	
	snp_f->free_func(snp_f);
	region_f->free_func(region_f);
    free_filter_chain(chain);
}

/* ******************************
 *       Checked fixtures       *
 * ******************************/

void create_passed_failed(void)
{
	failed = (list_t*) malloc (sizeof(list_t));
	list_init("failed", 1, MAX_RECORDS, failed);
}

void free_passed_failed(void)
{
//     list_free_deep(passed, vcf_record_free);
//     list_free_deep(failed, vcf_record_free);
	list_item_t* item = NULL;
	
	while ( (item = list_remove_item_async(passed)) != NULL ) {
		list_item_free(item);
	}
	
	item = NULL;
	while ( (item = list_remove_item_async(failed)) != NULL ) {
		list_item_free(item);
	}
	free(failed);
}


/* ******************************
 *           Unit tests         *
 * ******************************/

START_TEST (snp_include)
{
	((snp_filter_args*) snp_f->args)->include_snps = 1;
	passed = snp_f->filter_func(datasuite, failed, snp_f->args);
	
	// size(accepted) = SNPS_IN_FILE
	fail_unless(passed->length == SNPS_IN_FILE, "The number of SNP recognized is not correct");
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->length + failed->length == datasuite->length,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	list_item_t *item = NULL;
	// no accepted ID = '.'
	for (item = passed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		fail_if(strcmp(".", record->id) == 0, "A known SNP must have an ID");
	}
	
	// no rejected ID != '.'
	for (item = failed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		fail_unless(strcmp(".", record->id) == 0, "An unknown SNP can't have an ID defined");
	}
}
END_TEST


START_TEST (snp_exclude)
{
	((snp_filter_args*) snp_f->args)->include_snps = 0;
	passed = snp_f->filter_func(datasuite, failed, snp_f->args);
	
	// size(failed) = SNPS_IN_FILE
	fail_unless(failed->length == SNPS_IN_FILE, "The number of SNP recognized is not correct");
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->length + failed->length == datasuite->length,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	list_item_t *item = NULL;
	// no accepted ID != '.'
	for (item = passed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		fail_unless(strcmp(".", record->id) == 0, 
			"An unknown SNP can't have an ID defined");
	}
	
	// no rejected ID = '.'
	for (item = failed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		fail_if(strcmp(".", record->id) == 0, 
			"A known SNP must have an ID");
	}
}
END_TEST


START_TEST (region_chrom_1)
{	
	// create filter for just one chromosome
	char input[] = "1";
	region_f = create_region_filter(input, 0);
	passed = region_f->filter_func(datasuite, failed, region_f->args);
	
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->length + failed->length == datasuite->length,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	list_item_t *item = NULL;
	// no accepted chromosome != 1
	for (item = passed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		fail_unless(strcmp("1", record->chromosome) == 0, 
			"The record must be in chromosome 1");
	}
	
	// no rejected chromosome == 1
	for (item = failed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		fail_if(strcmp("1", record->chromosome) == 0, 
			"The record must not be in chromosome 1");
	}
	
	region_f->free_func(region_f);
}
END_TEST


START_TEST (region_chrom_1_2)
{	
	// create filter for both chromosomes in test file
	char input[] = "1,2";
	region_f = create_region_filter(input, 0);
	passed = region_f->filter_func(datasuite, failed, region_f->args);
	
	// all records pass
	fail_if(passed->length < datasuite->length, "All records must pass the test");
	fail_if(failed->length > 0, "There must not be rejected records");
	
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->length + failed->length == datasuite->length,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	list_item_t *item = NULL;
	// all records must be in chrom 1/2
	for (item = passed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		fail_unless(strcmp("1", record->chromosome) == 0 || strcmp("2", record->chromosome) == 0, 
			"The record must be in chromosome 1 or 2");
	}
	
	region_f->free_func(region_f);
}
END_TEST


START_TEST (region_chrom_start)
{
	char input[] = "1:10000000,2:20000000";
	region_f = create_region_filter(input, 0);
	passed = region_f->filter_func(datasuite, failed, region_f->args);
	
	// rejected = 18
	fail_unless(failed->length == 18, "18 records must be discarded");
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->length + failed->length == datasuite->length,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	list_item_t *item = NULL;
	// all accepted records must be in chrom 1 after position 10M, or in chrom 2 after position 20M
	int chrom1_ok, chrom2_ok;
	for (item = passed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		chrom1_ok = strcmp("1", record->chromosome) == 0 && record->position >= 10000000;
		chrom2_ok = strcmp("2", record->chromosome) == 0 && record->position >= 20000000;
		fail_unless(chrom1_ok || chrom2_ok, 
			"The record must be in chr1 pos >= 10M, or in chr2 pos >= 20M");
	}
	
	// all rejected records must be neither in chrom 1 after position 10M, or in chrom 2 after position 20M
	for (item = failed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		chrom1_ok = strcmp("1", record->chromosome) == 0 && record->position >= 10000000;
		chrom2_ok = strcmp("2", record->chromosome) == 0 && record->position >= 20000000;
		fail_if(chrom1_ok || chrom2_ok, 
			"The record must be neither in chr1 pos >= 10M or chr2 pos >= 20M");
	}
	
	region_f->free_func(region_f);
}
END_TEST


START_TEST (region_chrom_start_end)
{
	char input[] = "1:1000000-6000000,2:24000000-38000000";
	region_f = create_region_filter(input, 0);
	passed = region_f->filter_func(datasuite, failed, region_f->args);
	
	// accepted = 14
	fail_unless(passed->length == 22, "22 records must be accepted");
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->length + failed->length == datasuite->length,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	list_item_t *item = NULL;
	// all accepted records must be in chrom 1 after position (1M, 6M), or in chrom 2 after position (24M, 38M)
	int chrom1_ok, chrom2_ok;
	for (item = passed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		chrom1_ok = strcmp("1", record->chromosome) == 0 && record->position >= 1000000 && record->position <= 6000000;
		chrom2_ok = strcmp("2", record->chromosome) == 0 && record->position >= 24000000 && record->position <= 38000000;
		fail_unless(chrom1_ok || chrom2_ok, 
			"The record must be in chr1 pos in (1M, 6M), or in chr2 pos in (24M, 38M)");
	}
	
	// all rejected records must be neither in chrom 1 in position (1M, 6M), or in chrom 2 in position (24M, 38M)
	for (item = failed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		chrom1_ok = strcmp("1", record->chromosome) == 0 && record->position >= 1000000 && record->position <= 6000000;
		chrom2_ok = strcmp("2", record->chromosome) == 0 && record->position >= 24000000 && record->position <= 38000000;
		fail_if(chrom1_ok || chrom2_ok, 
			"The record must be neither in chr1 pos in (1M, 6M), or in chr2 pos in (24M, 38M)");
	}
	
	region_f->free_func(region_f);
}
END_TEST


START_TEST (quality_limit_bound)
{
    ((quality_filter_args*) quality_f->args)->min_quality = 100;
    passed = quality_f->filter_func(quality_datasuite, failed, quality_f->args);
    
    fail_unless(passed->length == 31, "Q100: The number of qualities found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->length + failed->length == quality_datasuite->length,
            "Q100: The sum of the number of accepted and rejected records must be the same as the input records");
    
    list_item_t *item = NULL;
    // no accepted ID < min_qual
    for (item = passed->first_p; item != NULL; item = item->next_p)
    {
        vcf_record_t *record = item->data_p;
        fail_if(record->quality < 100, "Q100: An accepted record can't have less quality than specified");
    }
    
    // no rejected ID >= min_qual
    for (item = failed->first_p; item != NULL; item = item->next_p)
    {
        vcf_record_t *record = item->data_p;
        fail_unless(record->quality < 100, "Q100: A rejected record can't have greater or equal quality than specified");
    }
}
END_TEST

START_TEST (quality_limit_over_bound)
{
    ((quality_filter_args*) quality_f->args)->min_quality = 101;
    passed = quality_f->filter_func(quality_datasuite, failed, quality_f->args);
    
    fail_unless(passed->length == 2, "Q101: The number of qualities found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->length + failed->length == quality_datasuite->length,
            "Q101: The sum of the number of accepted and rejected records must be the same as the input records");
    
    list_item_t *item = NULL;
    // no accepted ID < min_qual
    for (item = passed->first_p; item != NULL; item = item->next_p)
    {
        vcf_record_t *record = item->data_p;
        fail_if(record->quality < 101, "Q101: An accepted record can't have less quality than specified");
    }
    
    // no rejected ID >= min_qual
    for (item = failed->first_p; item != NULL; item = item->next_p)
    {
        vcf_record_t *record = item->data_p;
        fail_unless(record->quality < 101, "Q101: A rejected record can't have greater or equal quality than specified");
    }
}
END_TEST

START_TEST (quality_all_excluded)
{
    ((quality_filter_args*) quality_f->args)->min_quality = 200;
    passed = quality_f->filter_func(quality_datasuite, failed, quality_f->args);
    
    fail_unless(passed->length == 0, "Q200: The number of qualities found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->length + failed->length == quality_datasuite->length,
            "Q200: The sum of the number of accepted and rejected records must be the same as the input records");
    
    list_item_t *item = NULL;
    // no rejected ID >= min_qual
    for (item = failed->first_p; item != NULL; item = item->next_p)
    {
        vcf_record_t *record = item->data_p;
        fail_unless(record->quality < 200, "Q200: A rejected record can't have greater or equal quality than specified");
    }
}
END_TEST

START_TEST (quality_all_included)
{
    ((quality_filter_args*) quality_f->args)->min_quality = 60;
    passed = quality_f->filter_func(quality_datasuite, failed, quality_f->args);
    
    fail_unless(passed->length == 32, "Q60: The number of qualities found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->length + failed->length == quality_datasuite->length,
            "Q60: The sum of the number of accepted and rejected records must be the same as the input records");
    
    list_item_t *item = NULL;
    // no accepted ID < min_qual
    for (item = passed->first_p; item != NULL; item = item->next_p)
    {
        vcf_record_t *record = item->data_p;
        fail_if(record->quality < 60, "Q60: An accepted record can't have less quality than specified");
    }
}
END_TEST


START_TEST (snpinclude_regionchromstartend_chain)
{
	int num_filters;
	filter_t **filters = sort_filter_chain(chain, &num_filters);
	
	fail_unless(filters[0]->type == SNP, "The first filter to apply must be a SNP filter");
	fail_unless(filters[1]->type == REGION, "The second filter to apply must be a region filter");
	
	passed = run_filter_chain(datasuite, failed, filters, num_filters);
	
	// accepted = 5
	fail_unless(passed->length == 5, "5 records must be accepted");
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->length + failed->length == datasuite->length,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	list_item_t *item = NULL;
	// all accepted records must be a SNP in chr1 pos in (1M, 5M)
	int snp_ok, region_ok;
	for (item = passed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		snp_ok = strcmp(".", record->id) != 0;
		region_ok = strcmp("1", record->chromosome) == 0 && record->position >= 1000000 && record->position <= 5000000;
		fail_unless(region_ok && snp_ok, 
			"The record must be a SNP in chr1 pos in (1M, 5M)");
	}
	
	// all rejected records must not be a SNP in chr1 pos in (1M, 5M)
	for (item = failed->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		snp_ok = strcmp(".", record->id) != 0;
		region_ok = strcmp("1", record->chromosome) == 0 && record->position >= 1000000 && record->position <= 5000000;
		fail_if(region_ok && snp_ok, 
			"The record must not be a SNP in chr1 pos in (1M, 5M)");
	}
}
END_TEST



/* ******************************
 * 	Main entry point	*
 * ******************************/

int main (int argc, char *argv)
{
	vcf_file_t *file = vcf_open("CEU.exon.2010_03.genotypes__head400.vcf");
    vcf_file_t *quality_file = vcf_open("qualities.vcf");
	list_t *batches = read_test_datasuite(file);
	datasuite = batches->first_p->data_p;
    batches = read_test_datasuite(quality_file);
    quality_datasuite = batches->first_p->data_p;
	
	Suite *fs = create_test_suite();
	SRunner *fs_runner = srunner_create(fs);
	srunner_run_all(fs_runner, CK_NORMAL);
	int number_failed = srunner_ntests_failed (fs_runner);
	srunner_free (fs_runner);
	
	free_test_datasuite(datasuite, file);	// TODO exceeds check timeout
    free_test_datasuite(quality_datasuite, quality_file);
	vcf_close(file);
    vcf_close(quality_file);
	
	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


Suite *create_test_suite()
{
	// SNP filter (include and exclude)
	TCase *tc_snp = tcase_create("SNP filters");
	tcase_add_unchecked_fixture(tc_snp, setup_snp, teardown_snp);
	tcase_add_checked_fixture(tc_snp, create_passed_failed, free_passed_failed);
	tcase_add_test(tc_snp, snp_include);
	tcase_add_test(tc_snp, snp_exclude);
	
	// Region filter (chromosome, chrom+start position, chrom+start+end position)
	TCase *tc_region = tcase_create("Region filters");
	tcase_add_unchecked_fixture(tc_region, setup_region, teardown_region);
	tcase_add_checked_fixture(tc_region, create_passed_failed, free_passed_failed);
	tcase_add_test(tc_region, region_chrom_1);
	tcase_add_test(tc_region, region_chrom_1_2);
	tcase_add_test(tc_region, region_chrom_start);
	tcase_add_test(tc_region, region_chrom_start_end);
	
    // Quality filter
    TCase *tc_quality = tcase_create("Quality filters");
    tcase_add_unchecked_fixture(tc_quality, setup_quality, teardown_quality);
    tcase_add_checked_fixture(tc_quality, create_passed_failed, free_passed_failed);
    tcase_add_test(tc_quality, quality_limit_bound);
    tcase_add_test(tc_quality, quality_limit_over_bound);
    tcase_add_test(tc_quality, quality_all_included);
    tcase_add_test(tc_quality, quality_all_excluded);
    
	// Chains of filter (SNP+region...)
	TCase *tc_filterchain = tcase_create("Filter chains");
	tcase_add_unchecked_fixture(tc_filterchain, setup_snp_region, teardown_snp_region);
	tcase_add_checked_fixture(tc_filterchain, create_passed_failed, free_passed_failed);
	tcase_add_test(tc_filterchain, snpinclude_regionchromstartend_chain);
	
	// Add test cases to a test suite
	Suite *fs = suite_create("VCF filters");
	suite_add_tcase(fs, tc_snp);
	suite_add_tcase(fs, tc_region);
    suite_add_tcase(fs, tc_quality);
	suite_add_tcase(fs, tc_filterchain);
	
	return fs;
}

list_t *read_test_datasuite(vcf_file_t *file)
{
	list_t *batches = (list_t*) malloc (sizeof(list_t));
	list_init("batches", 1, 2, batches);
	
	int read = vcf_read_batches(batches, MAX_RECORDS, file, 0);
	if (read != 0)
	{
		fprintf(stderr, "Error reading file\n");
		return batches;
	}
	
	printf("Read %zu/%zu batches\n", batches->length, batches->max_length);
	
	printf("Batch contains %zu records\n", ((vcf_batch_t*) batches->first_p->data_p)->length);
	list_decr_writers(batches);
	
	return batches;
}

void free_test_datasuite(list_t *datasuite, vcf_file_t *file)
{
	// Free batches
	list_item_t* item = NULL;
	while ( (item = list_remove_item_async(datasuite)) != NULL ) {
		vcf_batch_free(item->data_p);
		list_item_free(item);
	}
}

