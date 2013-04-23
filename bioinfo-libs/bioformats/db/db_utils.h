#ifndef DB_UTILS_H
#define DB_UTILS_H

/*
 * db_utils.h
 *
 *  Created on: Apr 22, 2013
 *      Author: jtarraga
 */

//------------------------------------------------------------------------

#include "commons/sqlite/sqlite3.h"


//------------------------------------------------------------------------

int create_stats_db(const char *db_name, int chunksize, 
		    int (*record_stats_creation)(sqlite3 *), sqlite3** db);

//------------------------------------------------------------------------

int insert_global_stats(const char *id, const char *title, 
			const char *value, sqlite3 *db);

int prepare_statement_global_stats(sqlite3 *db, sqlite3_stmt **stmt);

int finalize_statement_global_stats(sqlite3_stmt *stmt);

int insert_statement_global_stats(const char *id, const char *title, 
				  const char *value, sqlite3_stmt *stmt, 
				  sqlite3 *db);

//------------------------------------------------------------------------

int insert_chunk(const char *chromosome, int chunk_id, int start, int end, 
		 int features_count, sqlite3 *db);

int inc_chunk(const char *chr, int chunk_id, int chunk_start, int chunk_end, 
	      sqlite3 *db);

//------------------------------------------------------------------------

int insert_record_stats(const char *chr, int chunk, void *stats,
			int (*insert_format_fn)(void *, sqlite3 *), sqlite3* db);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // end of DB_UTILS_H
