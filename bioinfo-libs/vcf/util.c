#include "util.h"

size_t count_regions(char *regions_string)
{
	size_t num_regions = 0;
	char *aux = regions_string;
	while (*aux) { 
		if (*aux++ == ',') ++num_regions;
	}
	
	return ++num_regions;
}

region_table_t *parse_regions(char *input_regions)//, size_t num_regions)
{
	region_table_t *regions_table = create_table(NULL);
	
	char *str_1 = input_regions;
    char *str_2 = (char*) malloc (64 * sizeof(char));
	char *saveptr1, *saveptr2;
	char *token, *subtoken;
	size_t token_len, subtoken_len;
	
	int i = 0;
	while ((token = strtok_r(str_1, ",", &saveptr1)) != NULL)
	{
		region_t *region = (region_t*) malloc (sizeof(region_t));
		token_len = strlen(token);
		
		dprintf("token = %s, len = %zu\n", token, token_len);
		
		strncpy(str_2, token, 63);
		str_2[token_len] = '\0';
		
		// Set chromosome
		subtoken = strtok_r(str_2, ":", &saveptr2);
		subtoken_len = strlen(subtoken);
		region->chromosome = (char*) malloc ((subtoken_len+1) * sizeof(char));
		strncpy(region->chromosome, subtoken, subtoken_len);
		region->chromosome[subtoken_len] = '\0';
		
		dprintf("region %s", region->chromosome);
		
		// Set start position
		subtoken = strtok_r(NULL, "-", &saveptr2);
		region->start_position = (subtoken != NULL) ? atol(subtoken) : 1;
		
		dprintf(":%u", region->start_position);
		
		// Set end position
		subtoken = strtok_r(NULL, "-", &saveptr2);
		region->end_position = (subtoken != NULL) ? atol(subtoken) : UINT_MAX;
		
		dprintf("-%u\n", region->end_position);
		
		insert_region(region, regions_table);
		
        str_1 = NULL;
        
		i++;
	}
	
	free(str_1); 
	free(str_2);
	
	return regions_table;
}
