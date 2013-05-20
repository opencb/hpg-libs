#include <stdlib.h>
#include <check.h>

#include <bioformats/features/region/region_table.h>
#include <commons/file_utils.h>
#include <containers/array_list.h>


Suite *create_test_suite();


static region_table_t *table;

static region_t *reg_2, *reg_3, *reg_4_1, *reg_4_2, *reg_X_1, *reg_X_2, *reg_Y;


/* ******************************
 * 	Unchecked fixtures	*
 * ******************************/

void setup_region_table(void) {
    delete_files_by_extension("/tmp/", ".db");
    table = create_region_table("http://ws.bioinfo.cipf.es/", "hsa", "latest");
}

void teardown_region_table(void) {
    free_region_table(table);
}

void setup_regions(void) {
    reg_2 = region_new(strdup("2"), 10000, 20000, "+", NULL);
    reg_3 = region_new(strdup("3"), 30000, 40000, "-", NULL);
    
    reg_4_1 = region_new(strdup("4"), 60000, 80000, NULL, "regulatory");
    reg_4_2 = region_new(strdup("4"), 200000, 250000, NULL, NULL);
    
    reg_X_1 = region_new(strdup("X"), 1000000, 1500000, NULL, NULL);
    reg_X_2 = region_new(strdup("X"), 1000000, 1800000, NULL, NULL);   // Different end than the previous

    reg_Y = region_new(strdup("Y"), 2000000, 3000000, NULL, NULL);
}

void teardown_regions(void) { }


/* ******************************
 * 	     Unit tests	        *
 * ******************************/

/* Data structure initialization */

START_TEST(table_structure) {
    fail_unless(table->max_chromosomes == 25, "The number of chromosomes for HSA is 25");
    for (int i = 0; i < 25; i++) {
        fail_if(count_regions_in_chromosome(table->ordering[i], table) > 0, "There must be no elements in the table after its creation");
    }
}
END_TEST


/* Data structure manipulation */

START_TEST(insert_region_and_chromosome) {
    fail_if(insert_region(reg_2, table), "Region reg_2 could not be inserted");
    
    array_list_t *list = get_chromosome("2", table);
    region_t *region = array_list_get(0, list);
    fail_if(strcmp("2", region->chromosome), "The inserted region must be in chromosome 2");
    fail_if(10000 != region->start_position || 20000 != region->end_position, "The inserted region must be in interval 10000-20000");
    
    fail_if(!count_regions_in_chromosome("2", table), "There must be one region in chromosome 2");
    fail_if(!find_exact_region(reg_2, table), "Region 2:10000-20000 must be inserted");
}
END_TEST

START_TEST(insert_several_regions) {
    // Insert regions in different chromosomes, in the same chromosome and in the same start position
    region_t *regions[] = { reg_2, reg_3, reg_4_1, reg_4_2, reg_X_1, reg_X_2, reg_Y };
    int num_regions = 7;
    fail_if(insert_regions(regions, num_regions, table), "Regions must be successfully inserted");

    // Check number of elements inserted in each chromosome
    fail_unless(count_regions_in_chromosome("2", table) == 1, "There must be 1 element(s) in chr2");
    fail_unless(count_regions_in_chromosome("3", table) == 1, "There must be 1 element(s) in chr3");
    fail_unless(count_regions_in_chromosome("4", table) == 2, "There must be 2 element(s) in chr4");
    fail_unless(count_regions_in_chromosome("X", table) == 2, "There must be 2 element(s) in chrX");
    fail_unless(count_regions_in_chromosome("Y", table) == 1, "There must be 1 element(s) in chrY");
    
    fail_unless(get_chromosome("2", table)->size == 1, "There must be 1 element(s) in chr2");
    fail_unless(get_chromosome("3", table)->size == 1, "There must be 1 element(s) in chr3");
    fail_unless(get_chromosome("4", table)->size == 2, "There must be 2 element(s) in chr4");
    fail_unless(get_chromosome("X", table)->size == 2, "There must be 2 element(s) in chrX");
    fail_unless(get_chromosome("Y", table)->size == 1, "There must be 1 element(s) in chrY");
}
END_TEST

START_TEST(insert_several_regions_one_by_one) {
    // Insert regions in different chromosomes, in the same chromosome and in the same start position
    fail_if(insert_region(reg_2, table), "Insertion of region in chr2 must be successfully performed");
    fail_if(insert_region(reg_3, table), "Insertion of region in chr3 must be successfully performed");
    fail_if(insert_region(reg_4_1, table), "Insertion of region in chr4, position 60K must be successfully performed");
    fail_if(insert_region(reg_4_2, table), "Insertion of region in chr4, position 200-250K must be successfully performed");
    fail_if(insert_region(reg_X_1, table), "Insertion of region in chrX, position 1M-1.5M must be successfully performed");
    fail_if(insert_region(reg_X_2, table), "Insertion of region in chrX, position 1M-1.8M must be successfully performed");
    fail_if(insert_region(reg_Y, table), "Insertion of region in chrY must be successfully performed");

    // Check number of elements inserted in each chromosome
    fail_unless(count_regions_in_chromosome("2", table) == 1, "There must be 1 element(s) in chr2");
    fail_unless(count_regions_in_chromosome("3", table) == 1, "There must be 1 element(s) in chr3");
    fail_unless(count_regions_in_chromosome("4", table) == 2, "There must be 2 element(s) in chr4");
    fail_unless(count_regions_in_chromosome("X", table) == 2, "There must be 2 element(s) in chrX");
    fail_unless(count_regions_in_chromosome("Y", table) == 1, "There must be 1 element(s) in chrY");
    
    fail_unless(get_chromosome("2", table)->size == 1, "There must be 1 element(s) in chr2");
    fail_unless(get_chromosome("3", table)->size == 1, "There must be 1 element(s) in chr3");
    fail_unless(get_chromosome("4", table)->size == 2, "There must be 2 element(s) in chr4");
    fail_unless(get_chromosome("X", table)->size == 2, "There must be 2 element(s) in chrX");
    fail_unless(get_chromosome("Y", table)->size == 1, "There must be 1 element(s) in chrY");
}
END_TEST

START_TEST(remove_regions_and_chromosomes) {
    // Insert regions
    region_t *regions[] = { reg_2, reg_3, reg_4_1, reg_4_2, reg_X_1, reg_X_2, reg_Y };
    int num_regions = 7;
    fail_if(insert_regions(regions, num_regions, table), "Regions must be successfully inserted");
    finish_region_table_loading(table);
    
    // Remove region by exact coordinates
    fail_if(remove_exact_region(reg_3, table), "Deletion of region in chr3 must be successfully performed");
    fail_if(count_regions_in_chromosome("3", table), "There must be 0 element(s) in chr3");
    
    // Remove all features inside region limits
    fail_if(remove_region(reg_X_2, table), "Deletion of regions in chrX must be successfully performed");
    fail_if(count_regions_in_chromosome("X", table), "There must be 0 element(s) in chrX");
    
    // Remove all features inside region limits
    fail_if(remove_region(reg_4_2, table), "Deletion of regions in chr4 must be successfully performed");
    fail_if(count_regions_in_chromosome("4", table) != 1, "There must be 1 element(s) in chr4");
}
END_TEST


/* Data structure search */

START_TEST(search_region) {
    // Insert regions
    fail_if(insert_region(reg_2, table), "Insertion of region in chr2 must be successfully performed");
    fail_if(insert_region(reg_3, table), "Insertion of region in chr3 must be successfully performed");
    fail_if(insert_region(reg_4_1, table), "Insertion of region in chr4, position 60-80K must be successfully performed");
    fail_if(insert_region(reg_4_2, table), "Insertion of region in chr4, position 200-250K must be successfully performed");
    fail_if(insert_region(reg_X_1, table), "Insertion of region in chrX, position 1M-1.5M must be successfully performed");
    fail_if(insert_region(reg_X_2, table), "Insertion of region in chrX, position 1M-1.8M must be successfully performed");
    fail_if(insert_region(reg_Y, table), "Insertion of region in chrY must be successfully performed");

    // Region in the limits of chr3 tree
    region_t t_reg_3 = { .chromosome = "3", .start_position = 30000, .end_position = 40000 };
    fail_if(!find_region(&t_reg_3, table), "Region 3:30000-40000 must be found");
    
    // Region contained in chrY tree
    region_t t_reg_Y = { .chromosome = "Y", .start_position = 2100000, .end_position = 2200000 };
    fail_if(!find_region(&t_reg_Y, table), "Region Y:2100000-2200000 must be found");
    
    // Region before the ones in chr4
    region_t t_reg_4_1 = { .chromosome = "4", .start_position = 30000, .end_position = 50000 };
    fail_if(find_region(&t_reg_4_1, table), "Region 4:30000-50000 must not be found");
    
    // Regions between the ones in chr4
    region_t t_reg_4_2 = { .chromosome = "4", .start_position = 30000, .end_position = 65000 };
    fail_if(!find_region(&t_reg_4_2, table), "Region 4:30000-65000 must be found");
    
    region_t t_reg_4_3 = { .chromosome = "4", .start_position = 60000, .end_position = 70000 };
    fail_if(!find_region(&t_reg_4_3, table), "Region 4:60000-70000 must be found");
    
    region_t t_reg_4_4 = { .chromosome = "4", .start_position = 50000, .end_position = 90000 };
    fail_if(!find_region(&t_reg_4_4, table), "Region 4:50000-90000 must be found");
    
    region_t t_reg_4_5 = { .chromosome = "4", .start_position = 80000, .end_position = 90000 };
    fail_if(!find_region(&t_reg_4_5, table), "Region 4:80000-90000 must be found");
    
    // Region after the ones in chr4
    region_t t_reg_4_6 = { .chromosome = "4", .start_position = 80001, .end_position = 90000 };
    fail_if(find_region(&t_reg_4_6, table), "Region 4:80001-90000 must not be found");
}
END_TEST


/* ******************************
 * 	Main entry point	*
 * ******************************/


int main(int argc, char *argv[]) {
    Suite *fs = create_test_suite();
    SRunner *fs_runner = srunner_create(fs);
    srunner_run_all(fs_runner, CK_NORMAL);
    int number_failed = srunner_ntests_failed(fs_runner);
    srunner_free (fs_runner);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Suite *create_test_suite() {
    // Creation of the table for storing regions
    TCase *tc_table = tcase_create("Table creation");
    tcase_add_unchecked_fixture(tc_table, setup_region_table, teardown_region_table);
    tcase_add_test(tc_table, table_structure);

    // Manipulation of the data structure
    TCase *tc_manipulation = tcase_create("Data structure manipulation");
    tcase_add_checked_fixture(tc_manipulation, setup_region_table, teardown_region_table);
    tcase_add_checked_fixture(tc_manipulation, setup_regions, teardown_regions);
    tcase_add_test(tc_manipulation, insert_region_and_chromosome);
    tcase_add_test(tc_manipulation, insert_several_regions);
    tcase_add_test(tc_manipulation, insert_several_regions_one_by_one);
    tcase_add_test(tc_manipulation, remove_regions_and_chromosomes);

    // Region searching
    TCase *tc_searching = tcase_create("Searching in region data structure");
    tcase_add_checked_fixture(tc_searching, setup_region_table, teardown_region_table);
    tcase_add_checked_fixture(tc_searching, setup_regions, teardown_regions);
    tcase_add_test(tc_searching, search_region);

    // Add test cases to a test suite
    Suite *fs = suite_create("Region searching table");
    suite_add_tcase(fs, tc_table);
    suite_add_tcase(fs, tc_manipulation);
    suite_add_tcase(fs, tc_searching);

    return fs;
}
