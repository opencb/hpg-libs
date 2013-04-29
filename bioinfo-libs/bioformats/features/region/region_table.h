#ifndef REGION_TABLE_H
#define REGION_TABLE_H

#include <ctype.h>
#include <stdlib.h>

#include <bioformats/db/db_utils.h>

#include <commons/log.h>
#include <commons/string_utils.h>
#include <commons/sqlite/sqlite3.h>

#include <containers/cprops/avl.h>
#include <containers/cprops/hashtable.h>
#include <containers/cprops/vector.h>

#include "region.h"

#define REGIONS_CHUNKSIZE       10000

/**
 * @struct region_table
 * 
 * Stores a collection of regions grouped by chromosome. It uses two data 
 * structures, one for guaranteing access to each chromosome in O(1) time, 
 * and another one for preserving the order among chromosomes.
 */
typedef struct region_table {
	int max_chromosomes;            /**< Number of chromosomes registered in the 'ordering' variable, which should be the maximum to keep in 'storage' */
	char **ordering;                /**< Order among chromosomes */
        sqlite3 *storage;               /**< Set of regions contained in the different chromosomes */
        khash_t(stats_chunks) *chunks;  /**< Hashtable containing groups of regions (chunks) */
        int is_ready;           /**< Flag that notifies that chunks are saved and the storage is indexed */
} region_table_t;


/**
 * Create the structure for storing a collection of chromosomes and regions
 * contained in them.
 * 
 * @param chromosome_file 
 * 	(Optional) File the list of ordered chromosomes is read from
 * @return A structure for a region table
 */
region_table_t *create_region_table(const char *url, const char *species, const char *version);

/**
 * Free the region table given as parameter.
 * 
 * @param regions region table to be freed
 */
void free_region_table(region_table_t *regions);


void finish_region_table_loading(region_table_t *table);


/* ******************************
 *  		Regions		*
 * ******************************/

int insert_regions(region_t **regions, int num_regions, region_table_t *table);

/**
 * Insert a region in a chromosome table. If the region's chromosome is not in the 
 * table, it is inserted before the region.
 * 
 * @param region region to insert
 * @param table region table where the region will be inserted in
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
 * @param region region presumably contained in the region table
 * @param table table that could contain the region
 * 
 * @return 1 if the region was found, 0 otherwise
 */
int contains_region(region_t *region, region_table_t *table);

/**
 * Find a region in a chromosome table. Its chromosome must be a perfect match, but its 
 * start and end positions can be just contained in the range of one of the regions 
 * stored. For example, if region 1:100-200 is stored, 1:150-160 would be a match.
 * 
 * @param region region to find
 * @param table table that could contain the region
 * 
 * @return 1 if the region matches another one, 0 otherwise
 */
int find_region(region_t *region, region_table_t *table);

/**
 * Remove a region from a chromosome table. In more than one region begin in the same 
 * position as this one, only the exact match is removed.
 * 
 * @param region region to remove
 * @param table table the region will be removed from
 * 
 * @return The removed region
 */
region_t *remove_region(region_t *region, region_table_t *table);


/* ******************************
 *  Chromosomes (region trees)	*
 * ******************************/

/**
 * Find a chromosome in the table and returns its related tree containing all the inserted regions.
 * 
 * @param key name of the chromosome to find
 * @param table regions structure the chromosome is get from
 * @return The tree with all regions related to the chromosome
 */
cp_avltree *get_chromosome(const char *key, region_table_t *table);

/**
 * Return the number of regions that have been inserted and belong to a certain chromosome.
 * 
 * @param key name of the chromosome whose regions are counted
 * @param table regions structure where the chromosome was inserted
 */
int count_regions_in_chromosome(const char *key, region_table_t *table);


#endif
