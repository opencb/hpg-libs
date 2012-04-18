#ifndef REGION_TABLE_H
#define REGION_TABLE_H

#include <ctype.h>

#include <cprops/avl.h>
#include <cprops/hashtable.h>
#include <cprops/vector.h>

#include "string_utils.h"
#include "region.h"


/**
 * @struct region_table
 * 
 * Stores a collection of regions grouped by chromosome. It uses two data 
 * structures, one for guaranteing access to each chromosome in O(1) time, 
 * and another one for preserving the order among chromosomes.
 * 
 * @var region_table::storage
 * 	Stores chromosomes and a collection of regions contained in them
 * @var region_table::ordering
 * 	Keeps the order among chromosomes
 * @var region_table::max_chromosomes
 * 	Number of chromosomes registered in the 'ordering' variable, which should be 
 * 	the maximum to keep in 'storage'
 */
typedef struct region_table
{
	int max_chromosomes;
	char **ordering;
	cp_hashtable *storage;
} region_table_t;


/**
 * Create the structure for storing a collection of chromosomes and regions
 * contained in them.
 * 
 * @param chromosome_file 
 * 	(Optional) File the list of ordered chromosomes is read from
 * 
 * @return
 * 	A structure
 */
region_table_t *create_table(const char *chromosome_file);


/* ******************************
 *  		Regions		*
 * ******************************/

/**
 * Insert a region in a chromosome table. If the region's chromosome is not in the 
 * table, it is inserted before the region.
 * 
 * @param region
 * 	Region to insert
 * @param table
 * 	Region table where the region will be inserted in
 * 
 * @return
 * 	0: If the region was successfully inserted
 * 	1: If the region's chromosome could not be inserted
 * 	2: If the region could not be inserted in the chromosome tree
 */
int insert_region(region_t *region, region_table_t *table);

/**
 * Check whether a region is stored in a chromosome table. Its chromosome, start and end 
 * position must be a perfect match of one of the regions already inserted.
 * 
 * @param region
 * 	Region presumably contained in the region table
 * @param table
 * 	Table that could contain the region
 * 
 * @return
 * 	1 if the region was found, 0 otherwise
 */
int contains_region(region_t *region, region_table_t *table);

/**
 * Find a region in a chromosome table. Its chromosome must be a perfect match, but its 
 * start and end positions can be just contained in the range of one of the regions 
 * stored. For example, if region 1:100-200 is stored, 1:150-160 would be a match.
 * 
 * @param region
 * 	Region to find
 * @param table
 * 	Table that could contain the region
 * 
 * @return
 * 	1 if the region matches another one, 0 otherwise
 */
int find_region(region_t *region, region_table_t *table);

/**
 * Remove a region from a chromosome table. In more than one region begin in the same 
 * position as this one, only the exact match is removed.
 * 
 * @param region
 * 	Region to remove
 * @param table
 * 	Table the region will be removed from
 * 
 * @return
 * 	The removed region
 */
region_t *remove_region(region_t *region, region_table_t *table);

/**
 * Free a region from the chromosome table.
 * 
 * @param region
 * 	Region to free
 */
void free_region(region_t *region);


/* ******************************
 *  Chromosomes (region trees)	*
 * ******************************/

/**
 * Insert a chromosome in the table, and initializes the tree of regions it contains.
 * 
 * @param key 
 * 	The name of the chromosome to be inserted
 * @param table
 * 	The regions structure the chromosome is inserted in
 * 
 * @return 
 * 	The created chromosome tree, or NULL if the chromosome could not be inserted
 */
void *insert_chromosome(const char *key, region_table_t *table);

/**
 * Find a chromosome in the table and returns its related tree containing all the inserted regions.
 * 
 * @param key
 * 	The name of the chromosome to find
 * 
 * @return
 * 	The tree with all regions related to the chromosome
 */
cp_avltree *get_chromosome(const char *key, region_table_t *table);

/**
 * Remove a chromosome from the table and returns the tree containing all the inserted regions.
 * 
 * @param key
 * 	The name of the chromosome to remove
 * 
 * @return
 * 	The removed tree, related to the chromosome
 */
cp_avltree *remove_chromosome(const char *key, region_table_t *table);

/**
 * Free a chromosome from the table and its related regions tree.
 * 
 * @param key
 * 	The name of the chromosome to free
 * 
 */
void free_chromosome(const char *key);

/**
 * Free the contents of a chromosome tree, being it inside the chromosome table or not.
 * 
 * @param chromosome_tree
 * 	The tree to be freed
 * 
 */
void free_chromosome_tree(cp_avltree *chromosome_tree);


#endif
