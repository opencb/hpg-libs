#include <stdlib.h>
#include <check.h>
#include <string.h>

#include "string_utils.h"
#include "array_list.h"

Suite *create_test_suite();

static array_list_t *array_list_p;

int compare_items2(const void *i1, const void *i2) {
//	void *i1 = (void*)item1;
//	void *i2 = (void*)item2;
	return strcmp((char*)i1, (char*)i2);
//	return item1 == item2;
}

/* ******************************
 * 	Unchecked fixtures	*
 * ******************************/

void setup_region_table(void)
{
//	array_list_p = array_list_new();
	array_list_p =array_list_new(3, 1.4, COLLECTION_MODE_SYNCHRONIZED);
}

void teardown_region_table(void)
{
	array_list_free(array_list_p, free);
}

void setup_regions(void)
{

}

void teardown_regions(void) { array_list_free(array_list_p, NULL);}

/* ******************************
 * 	     Unit tests	        *
 * ******************************/

/* Data structure manipulation */

START_TEST(array_list_insert_test)
{
	printf("array_list_insert_test:\n");
	void *item1 = strdup("item1");
	void *item2 = strdup("item2");
	void *item3 = strdup("item3");

	fail_if(array_list_insert(item1, array_list_p) != 1, "Insertion must be successfully performed");
	fail_unless(array_list_size(array_list_p) == 1, "There must be one element in the chromosome table");
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);
	
	fail_if(array_list_insert(item2, array_list_p) != 1, "Insertion must be successfully performed");
	fail_unless(array_list_size(array_list_p) == 2, "There must be one element in the chromosome table");
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);

	fail_if(array_list_insert(item3, array_list_p) != 1, "Insertion must be successfully performed");
	fail_unless(array_list_size(array_list_p) == 3, "There must be one element in the chromosome table");
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);

	fail_if(array_list_insert(item3, array_list_p) != 1, "Insertion must be successfully performed");
	fail_unless(array_list_size(array_list_p) == 4, "There must be one element in the chromosome table");
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);

	array_list_set(3, item1, array_list_p);
	fail_if(array_list_insert(item3, array_list_p) != 1, "Insertion must be successfully performed");
	fail_unless(array_list_size(array_list_p) == 5, "There must be one element in the chromosome table");
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);

	array_list_print(array_list_p);
	printf("\n\n");
}
END_TEST

START_TEST(array_list_insert_all_test)
{
	printf("array_list_insert_all_test:\n");
	void *item0 = strdup("item1");
	void *item1 = strdup("item2");
	void *item2 = strdup("item3");

	void **items = (void **) malloc(3*sizeof(void *));
	items[0] = item0;
	items[1] = item1;
	items[2] = item2;

	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);

	array_list_insert_all(items, 3, array_list_p);
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);
	array_list_print(array_list_p);
	printf("\n");

	array_list_insert_all(items, 3, array_list_p);
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);
	array_list_print(array_list_p);
	printf("\n\n");
}
END_TEST

START_TEST(array_list_remove_test)
{
	printf("array_list_remove_test:\n");
	void *item1 = strdup("item1");
	void *item2 = strdup("item2");
	void *item3 = strdup("item2");
	void *item4 = strdup("item3");
	void *item5 = strdup("item4");

	fail_if(array_list_insert(item1, array_list_p) != 1, "Insertion must be successfully performed");
	fail_unless(array_list_size(array_list_p) == 1, "There must be one element in the chromosome table");
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);

	fail_if(array_list_insert(item2, array_list_p) != 1, "Insertion must be successfully performed");
	fail_unless(array_list_size(array_list_p) == 2, "There must be one element in the chromosome table");
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);

	fail_if(array_list_insert(item3, array_list_p) != 1, "Insertion must be successfully performed");
	fail_unless(array_list_size(array_list_p) == 3, "There must be one element in the chromosome table");
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);

	fail_if(array_list_insert(item4, array_list_p) != 1, "Insertion must be successfully performed");
	fail_unless(array_list_size(array_list_p) == 4, "There must be one element in the chromosome table");
	printf("Capacity: %lu, size: %lu\n", array_list_p->capacity, array_list_p->size);


	size_t index = array_list_index_of(item4, array_list_p);
	fail_unless(index == 3, "Text 'item3' must be in position 3");
	printf("Index: %lu\n", index);

	array_list_p->compare_fn = compare_items2;
	index = array_list_index_of(item3, array_list_p);
	fail_unless(index == 1, "Text 'item2' must be in position 1");
	printf("Index: %lu\n", index);

	index = array_list_index_of(item5, array_list_p);
	fail_unless(index == ULONG_MAX, "Text 'item4' must not be in the collection");
	printf("Index: %lu\n", index);


	array_list_print(array_list_p);
	printf("\n");

	array_list_remove(item1, array_list_p);
	array_list_print(array_list_p);
	printf("\n");
}
END_TEST



/* ******************************
 * 	Main entry point	*
 * ******************************/


int main(int argc, char *argv[])
{
	Suite *fs = create_test_suite();
	SRunner *fs_runner = srunner_create(fs);
	srunner_run_all(fs_runner, CK_NORMAL);
	int number_failed = srunner_ntests_failed(fs_runner);
	srunner_free (fs_runner);
	
	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Suite *create_test_suite()
{
	// Creation of the table for storing regions
	// Manipulation of the data structure
	TCase *tc_manipulation = tcase_create("Data structure manipulation");
	tcase_add_checked_fixture(tc_manipulation, setup_region_table, teardown_region_table);
//	tcase_add_checked_fixture(tc_manipulation, setup_regions, teardown_regions);
	tcase_add_test(tc_manipulation, array_list_insert_test);
	tcase_add_test(tc_manipulation, array_list_insert_all_test);
	tcase_add_test(tc_manipulation, array_list_remove_test);
	
	// Add test cases to a test suite
	Suite *fs = suite_create("Region searching table");
	suite_add_tcase(fs, tc_manipulation);
	
	return fs;
}
