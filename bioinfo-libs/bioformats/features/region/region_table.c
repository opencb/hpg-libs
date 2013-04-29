#include "region_table.h"

region_table_t *create_region_table(const char *url, const char *species, const char *version)
{
    int num_chromosomes;
    region_table_t *table = (region_table_t*) malloc (sizeof(region_table_t));

    table->ordering = get_chromosome_order(url, species, version, &num_chromosomes);
    table->max_chromosomes = num_chromosomes;

    create_regions_db("/tmp/regions.db", REGIONS_CHUNKSIZE, &(table->storage));

    return table;
}

void free_region_table(region_table_t* regions) {
    // Close database
    sqlite3_close_v2(regions->storage);
    
    // Free ordering array
    char **ordering = regions->ordering;
    for (int i = 0; i < regions->max_chromosomes; i++) {
        free(ordering[i]);
    }
    free(ordering);

    free(regions);
}


void finish_region_table_loading(region_table_t *table) {
    // Save chunks from hashtable
    insert_chunk_hash(REGIONS_CHUNKSIZE, table->chunks, table->storage);
    // Index contents
    create_regions_index(table->storage);
    // Mark as ready!
    table->is_ready = 1;
}


/* ******************************
 *  		Regions		*
 * ******************************/

int insert_region(region_t *region, region_table_t *table) {
    insert_regions(&region, 1, table);
}

int insert_regions(region_t **regions, int num_regions, region_table_t *table) {
    int rc;
    sqlite3_stmt *stmt;
    sqlite3* db = table->storage;
    
    char *sql_begin = "BEGIN TRANSACTION";
    rc = exec_sql(sql_begin, db);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not insert regions: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }

    char sql[] = "INSERT INTO regions VALUES (?1, ?2, ?3, ?4, ?5)";
    rc = sqlite3_prepare_v2(db, sql, strlen(sql) + 300, &stmt, NULL);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not insert regions: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }
    
    for (int i = 0; i < num_regions; i++) {
        region_t *region = regions[i];
        sqlite3_bind_text(stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(stmt, 2, region->start_position);
        sqlite3_bind_int64(stmt, 3, region->end_position);
        sqlite3_bind_text(stmt, 4, region->strand, strlen(region->strand), SQLITE_STATIC);
        sqlite3_bind_text(stmt, 5, region->type, strlen(region->type), SQLITE_STATIC);

        if (rc = sqlite3_step(stmt) != SQLITE_DONE) {
            LOG_ERROR_F("Could not insert region %s:%ld-%ld: %s (%d)\n", 
                        region->chromosome, region->start_position, region->end_position, 
                        sqlite3_errmsg(db), sqlite3_errcode(db));
        } else {
            // Update value in chunks hashtable
            update_chunks_hash(region->chromosome, UINT_MAX, REGIONS_CHUNKSIZE, 
                               region->start_position, region->end_position, table->chunks);
        }

        sqlite3_reset(stmt);
    }
    
    char *sql_end = "END TRANSACTION";
    rc = exec_sql(sql_end, db);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not insert regions: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }

    return rc;
}


int contains_region(region_t *region, region_table_t *table) {
    if (!table->is_ready) { // Don't allow queries over an un-indexed DB
        finish_region_table_loading(table);
    }
    
    int count = 0;
    sqlite3* db = table->storage;
    
    char *sql_begin = "BEGIN TRANSACTION";
    exec_sql(sql_begin, db);

    sqlite3_stmt *query_stmt;
    char *sql = "SELECT COUNT(*) FROM regions WHERE chromosome = ?1 AND start = ?2 AND end = ?3";
    sqlite3_prepare_v2(db, sql, strlen(sql), &query_stmt, NULL);
    sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
    sqlite3_bind_int64(query_stmt, 2, region->start_position);
    sqlite3_bind_int64(query_stmt, 3, region->end_position);
    
    if (sqlite3_step(query_stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(query_stmt, 1);
    } else {
        // Retry once
        sqlite3_reset(query_stmt);
        sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(query_stmt, 2, region->start_position);
        sqlite3_bind_int64(query_stmt, 3, region->end_position);
        
        if (sqlite3_step(query_stmt) == SQLITE_ROW) {
            count = sqlite3_column_int(query_stmt, 1);
        } else {
            LOG_ERROR_F("Regions table failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }
    }
    
    sqlite3_finalize(query_stmt);
    
    char *sql_end = "END TRANSACTION";
    exec_sql(sql_end, db);

    return count > 0;
}


int find_region(region_t *region, region_table_t *table) {
    if (!table->is_ready) { // Don't allow queries over an un-indexed DB
        finish_region_table_loading(table);
    }
    
    int count = 0;
    sqlite3* db = table->storage;
    
    char *sql_begin = "BEGIN TRANSACTION";
    exec_sql(sql_begin, db);

    sqlite3_stmt *query_stmt;
    char *sql = "SELECT COUNT(*) FROM regions WHERE chromosome = ?1 AND start <= ?2 AND end >= ?3";
    sqlite3_prepare_v2(db, sql, strlen(sql), &query_stmt, NULL);
    sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
    sqlite3_bind_int64(query_stmt, 2, region->start_position);
    sqlite3_bind_int64(query_stmt, 3, region->end_position);
    
    if (sqlite3_step(query_stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(query_stmt, 1);
    } else {
        // Retry once
        sqlite3_reset(query_stmt);
        sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(query_stmt, 2, region->start_position);
        sqlite3_bind_int64(query_stmt, 3, region->end_position);
        
        if (sqlite3_step(query_stmt) == SQLITE_ROW) {
            count = sqlite3_column_int(query_stmt, 1);
        } else {
            LOG_ERROR_F("Regions table failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }
    }
    
    sqlite3_finalize(query_stmt);
    
    char *sql_end = "END TRANSACTION";
    exec_sql(sql_end, db);

    return count > 0;
}

region_t *remove_region(region_t *region, region_table_t *table) {
/*
	if (!cp_hashtable_contains(table->storage, region->chromosome))
	{
		return NULL;
	}
	
	region_t *removed = NULL;
	cp_avltree *chr = get_chromosome(region->chromosome, table);
	// - If the tree node contains >1 element: remove all elements in that exact range. 
	//   If all of them are erased, remove the node.
	// - If the tree node contains 1 element: remove the node
	cp_vector *region_vector = cp_avltree_get(chr, &region->start_position);
	int regions_in_node = cp_vector_size(region_vector);
	
	if (regions_in_node > 1)
	{
		int i = 0;
		do
		{
			// Remove all elements in that exact range. 
			if (compare_position_ranges(region, cp_vector_element_at(region_vector, i)) == 0)
			{
				removed = cp_vector_remove_element_at(region_vector, i);
				regions_in_node--;
			} else
			{
				i++;
			}
		} while (i < regions_in_node);
		
		// If all elements are erased, remove the node
		regions_in_node = cp_vector_size(region_vector);
		if (regions_in_node == 0)
		{
			removed = cp_avltree_delete(chr, &region->start_position);
		}
	} else
	{
		removed = cp_avltree_delete(chr, &region->start_position);
	}
	
	return removed;
*/
/*
    int rc;
    sqlite3_stmt *stmt;
    sqlite3* db = table->storage;
    
    char *sql_begin = "BEGIN TRANSACTION";
    rc = exec_sql(sql_begin, db);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not remove region: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }

    char sql[] = "INSERT INTO regions VALUES (?1, ?2, ?3, ?4, ?5)";
    rc = sqlite3_prepare_v2(db, sql, strlen(sql) + 300, &stmt, NULL);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not remove region: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }
    
    for (int i = 0; i < num_regions; i++) {
        region_t *region = regions[i];
        sqlite3_bind_text(stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(stmt, 2, region->start_position);
        sqlite3_bind_int64(stmt, 3, region->end_position);
        sqlite3_bind_text(stmt, 4, region->strand, strlen(region->strand), SQLITE_STATIC);
        sqlite3_bind_text(stmt, 5, region->type, strlen(region->type), SQLITE_STATIC);

        if (rc = sqlite3_step(stmt) != SQLITE_DONE) {
            LOG_ERROR_F("Could not insert region %s:%ld-%ld: %s (%d)\n", 
                        region->chromosome, region->start_position, region->end_position, 
                        sqlite3_errmsg(db), sqlite3_errcode(db));
        }

        sqlite3_reset(stmt);
    }
    
    char *sql_end = "END TRANSACTION";
    rc = exec_sql(sql_end, db);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not insert regions: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }

    return rc == SQLITE_OK;
*/
    return NULL;
}



/* ******************************
 *  Chromosomes (region trees)	*
 * ******************************/

cp_avltree *get_chromosome(const char *key, region_table_t *table) {
/*
    return (cp_avltree *) cp_hashtable_get_by_option(table->storage, (void*) key, COLLECTION_MODE_NOSYNC);
*/
    return NULL;
}

int count_regions_in_chromosome(const char *key, region_table_t *table) {
    if (!table->is_ready) { // Don't allow queries over an un-indexed DB
        finish_region_table_loading(table);
    }
    
    int count = -1;
    sqlite3* db = table->storage;
    
    char *sql_begin = "BEGIN TRANSACTION";
    exec_sql(sql_begin, db);

    sqlite3_stmt *query_stmt;
    char *sql = "SELECT COUNT(*) FROM regions WHERE chromosome = ?1";
    sqlite3_prepare_v2(db, sql, strlen(sql), &query_stmt, NULL);
    sqlite3_bind_text(query_stmt, 1, key, strlen(key), SQLITE_STATIC);
    
    if (sqlite3_step(query_stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(query_stmt, 1);
    } else {
        // Retry once
        sqlite3_reset(query_stmt);
        sqlite3_bind_text(query_stmt, 1, key, strlen(key), SQLITE_STATIC);
        if (sqlite3_step(query_stmt) == SQLITE_ROW) {
            count = sqlite3_column_int(query_stmt, 1);
        } else {
            LOG_ERROR_F("Regions table failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }
    }
    
    sqlite3_finalize(query_stmt);
    
    char *sql_end = "END TRANSACTION";
    exec_sql(sql_end, db);

    return count;
}

