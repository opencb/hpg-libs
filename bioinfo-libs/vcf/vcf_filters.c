#include "vcf_filters.h"


filter_t *create_snp_filter(char *include_snps)
{
	filter_t *snp_f =  (filter_t*) malloc (sizeof(filter_t));
	snp_f->type = SNP;
	snp_f->filter_func = snp_filter;
	snp_f->free_func = free_snp_filter;
	snp_f->priority = 5;
	
	snp_filter_args *snp_args = (snp_filter_args*) malloc (sizeof(snp_filter_args));
	snp_args->include_snps = 1;	// Default: Include SNPs
	
	if (include_snps != NULL)
	{
		if (strcmp("include", include_snps) == 0)
		{
			snp_args->include_snps = 1;
		} else if (strcmp("exclude", include_snps) == 0)
		{
			snp_args->include_snps = 0;
		}
	}
	
	snp_f->args = snp_args;
	
	return snp_f;
}

void free_snp_filter(filter_t *filter)
{
	free(filter->args);
	free(filter);
}

filter_t *create_region_filter(char *region_descriptor, int use_region_file)
{
	filter_t *region_f = (filter_t*) malloc (sizeof(filter_t));
	region_f->type = REGION;
	region_f->filter_func = region_filter;
	region_f->free_func = free_region_filter;
	region_f->priority = 2;
	
	region_filter_args *reg_args = (region_filter_args*) malloc (sizeof(region_filter_args));
	if (use_region_file)
	{
		// TODO read regions from file
		reg_args->regions = parse_regions(region_descriptor);
	} else
	{
		reg_args->regions = parse_regions(region_descriptor);
	}
	region_f->args = reg_args;
	
	return region_f;
}

void free_region_filter(filter_t *filter)
{
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

int filter_compare(const void *filter1, const void *filter2)
{
	return ((filter_t*) filter1)->priority - ((filter_t*) filter2)->priority;
}



filter_chain *add_to_filter_chain(filter_t *filter, filter_chain *chain)
{
	filter_chain *result = chain;
	
	if (result == NULL)
	{
		result = cp_heap_create((cp_compare_fn) filter_compare);
	}
	cp_heap_push(result, filter);
	
	return result;
}

filter_t **sort_filter_chain(filter_chain *chain, int *num_filters)
{
	*num_filters = cp_heap_count(chain);
	filter_t **filters = (filter_t**) malloc (cp_heap_count(chain) * sizeof(filter_t*));
	
	// Pop all filters from the heap and make an ordered list from them
	filter_t *filter = NULL;
	for (int i = 0; (filter = cp_heap_pop(chain)) != NULL; i++)
	{
		dprintf("Filter %d type = %d\n", i, filter->type);
		filters[i] = filter;
	}
	
	return filters;
}

list_t *run_filter_chain(list_t *input_records, list_t *failed, filter_t **filters, int num_filters)
{
	list_t *passed = input_records;
	list_t *aux_passed;
	
	dprintf("Applying filter chain of %d filters\n", num_filters);
	
	// Apply each filter with the arguments provided
	for (int i = 0; i < num_filters; i++)
	{
		filter_t *filter = filters[i];
		aux_passed = filter->filter_func(passed, failed, filter->args);
// 		free(passed);
		passed = aux_passed;
	}
	
	return passed;
}



list_t *region_filter(list_t *input_records, list_t *failed, void *f_args)
{
	char *field;
	list_t *passed = (list_t*) malloc (sizeof(list_t));
	list_init("passed", 1, INT_MAX, passed);
	
	region_filter_args *args = (region_filter_args*) f_args;
	region_table_t *regions = args->regions;
	
	dprintf("region_filter over %zu records\n", input_records->length);
	
	int i = 0;
	region_t *region = (region_t*) malloc (sizeof(region_t));
	for (list_item_t *item = input_records->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		list_item_t *new_item = list_item_new(item->id, item->type, record);
		
		dprintf("record = %s, %ld\n", record->chromosome, record->position);
		
		region->chromosome = record->chromosome;
		region->start_position = record->position;
		region->end_position = record->position;
		
		if (find_region(region, regions))
		{
			// Add to the list of records that pass all checks for at least one region
			list_insert_item(new_item, passed);
			dprintf("%s, %ld passed\n", record->chromosome, record->position);
		} else {
			// Add to the list of records that fail all checks for all regions
			list_insert_item(new_item, failed);
		}
		
		i++;
	}
	
	free(region);
	
	return passed;
}

list_t *snp_filter(list_t *input_records, list_t *failed, void *f_args)
{	
	list_t *passed = (list_t*) malloc (sizeof(list_t));
	list_init("passed", 1, input_records->max_length, passed);
	
	int include_snps = ((snp_filter_args*)f_args)->include_snps;
	
	dprintf("snp_filter (preserve SNPs = %d) over %zu records\n", include_snps, input_records->length);
	for (list_item_t *item = input_records->first_p; item != NULL; item = item->next_p)
	{
		vcf_record_t *record = item->data_p;
		list_item_t *new_item = list_item_new(item->id, item->type, record);
		// TODO check 'id' field is not empty (modifications to the parser needed!)
		if (strcmp(".", record->id) == 0)
		{
			if (include_snps) 
				list_insert_item(new_item, failed);
			else	
				list_insert_item(new_item, passed);
		} else
		{
			if (include_snps) 
				list_insert_item(new_item, passed);
			else	
				list_insert_item(new_item, failed);
		}
	}
	
	return passed;
}
