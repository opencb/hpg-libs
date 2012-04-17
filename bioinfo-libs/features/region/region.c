#include "region.h"


char **get_chromosome_order(const char *chromosome_file, int *num_chromosomes)
{
	char **ordering = NULL;
	int count = 0;
	if (chromosome_file == NULL)
	{
		count = 25;
		// Use default criteria: 1st numerics, then lexicographic order
		ordering = (char**) malloc (count * sizeof(char*));
		for (int i = 0; i < 22; i++)
		{
			ordering[i] = (char*) malloc (16 * sizeof(char));
			sprintf(ordering[i], "%d", i+1);
		}
		ordering[22] = (char*) malloc (16 * sizeof(char));
		ordering[23] = (char*) malloc (16 * sizeof(char));
		ordering[24] = (char*) malloc (16 * sizeof(char));
		strcpy(ordering[22], "X");
		strcpy(ordering[23], "Y");
		strcpy(ordering[24], "MT");
	} else
	{
		// Use list from file
		FILE *file = fopen(chromosome_file, "r");
		if (!file) {
			*num_chromosomes = 0;
			return ordering;
		}
	
		// Count chromosomes lines
		char line[16];
		while ( fgets(line, sizeof line, file) )
		{
			if ( line[0] != '\n' ) { ++count; }
		}
		if (count == 0) { 
			*num_chromosomes = 0;
			return ordering;
		}
		
		ordering = (char**) malloc (count * sizeof(char*));
		
		// Read file
		rewind(file);
		int aux_count = 0;
		while ( fgets(line, sizeof line, file) )
		{
			if ( line[0] != '\n' ) { 
				ordering[aux_count] = (char*) malloc (16 * sizeof(char));
				strcpy(ordering[aux_count], trim(line));
				aux_count++;
			}
		}
		
		// Close file
		fclose(file);
	}
	
	*num_chromosomes = count;
	
	return ordering;
}


int compare_regions(void *region_1, void *region_2, char **chromosome_ordering, int num_chromosomes)
{
	if (region_1 == NULL || region_2 == NULL)
	{
		return INT_MIN;
	}
	
	region_t *reg_1 = (region_t *) region_1;
	region_t *reg_2 = (region_t *) region_2;
	
	// TODO This could be avoided while inserting, becuase regions are classified by chromosome
	int result = compare_chromosomes(reg_1->chromosome, reg_2->chromosome, chromosome_ordering, num_chromosomes);
	if (result != 0)
	{
		return result;
	} else
	{
// 		return compare_position_ranges(reg_1, reg_2);
		return compare_positions(&reg_1->start_position, &reg_2->start_position);
	}
}

static int compare_chromosomes(char *chromosome_1, char *chromosome_2, char **chromosome_ordering, int num_chromosomes)
{
	int chr_1_found = 0, chr_2_found = 0;
	for (int i = 0; i < num_chromosomes; i++)
	{
		if (strcasecmp(chromosome_ordering[i], chromosome_1) == 0)
		{
			if (strcasecmp(chromosome_ordering[i], chromosome_2) == 0)
			{
				return 0;
			}
			return -1;
		}
		else if (strcasecmp(chromosome_ordering[i], chromosome_2) == 0)
		{
			return 1;
		}
	}
	return 0;
}

int compare_positions(uint32_t *position_1, uint32_t *position_2)
{
	return *position_1 - *position_2;
}

int compare_position_ranges(region_t *region_1, region_t *region_2)
{
	int result = region_1->start_position - region_2->start_position;
	if (result == 0)
	{
		result = region_1->end_position - region_2->end_position;
	}
	return result;
}

int region_contains_other(region_t *container, region_t *content)
{
// 	printf("container = %s:%d:%d\t", container->chromosome, container->start_position, container->end_position);
// 	printf("content = %s:%d:%d\n", content->chromosome, content->start_position, content->end_position);
// 	printf("start ok = %d, end ok = %d\n", container->start_position <= content->start_position, container->end_position >= content->end_position);
	
	int result = strcasecmp(container->chromosome, content->chromosome);
	if (result != 0) { return result; }
	
	if (container->start_position > content->start_position)
	{
		return 1;
	}
	
	if (container->end_position < content->end_position)
	{
		return -1;
	}
	
	return 0;
}
