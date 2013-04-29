#ifndef DB_UTILS_H
#define DB_UTILS_H

/*
 * db_utils.h
 *
 *  Created on: Apr 22, 2013
 *      Author: jtarraga
 */

//------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#include <bioformats/features/region/region.h>
#include <commons/log.h>
#include <commons/sqlite/sqlite3.h>
#include <containers/array_list.h>
#include <containers/khash.h>

//------------------------------------------------------------------------

KHASH_MAP_INIT_STR(stats_chunks, int)    

//------------------------------------------------------------------------
//-------------------------- STATISTICS DATABASE -------------------------

int create_stats_db(const char *db_name, int chunksize, 
		    int (*record_stats_creation)(sqlite3 *), sqlite3** db);

int create_stats_index(int (*create_custom_index)(sqlite3 *), sqlite3* db);

//------------------------------------------------------------------------

int insert_global_stats(const char *id, const char *title, 
			const char *value, sqlite3 *db);

int prepare_statement_global_stats(sqlite3 *db, sqlite3_stmt **stmt);

int finalize_statement_global_stats(sqlite3_stmt *stmt);

int insert_statement_global_stats(const char *id, const char *title, 
				  const char *value, sqlite3_stmt *stmt, 
				  sqlite3 *db);


//------------------------------------------------------------------------
//--------------------------- REGIONS DATABASE ---------------------------

int create_regions_db(const char *db_name, int chunksize, sqlite3** db);

int create_regions_index(sqlite3* db);

//------------------------------------------------------------------------

int insert_region_query_fields_list(array_list_t *list, sqlite3 *db);


//------------------------------------------------------------------------
//----------------------------- CHUNKS TABLE -----------------------------

int insert_chunk(const char *chromosome, int chunk_id, int start, int end, 
		 int features_count, sqlite3 *db);

int inc_chunk(const char *chr, int chunk_id, int chunk_start, int chunk_end, 
	      sqlite3 *db);

int update_chunks_hash(const char *chr, int chr_length, int chunksize, 
		       int start, int end, khash_t(stats_chunks) *hash);

int insert_chunk_hash(int chunksize, khash_t(stats_chunks) *hash, sqlite3 *db);


#endif // end of DB_UTILS_H
