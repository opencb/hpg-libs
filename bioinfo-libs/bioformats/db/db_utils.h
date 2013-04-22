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


int create_stats_db(const char *db_name, int chunksize, 
		    int (*format_stats)(sqlite3 *), sqlite3 **db);

//------------------------------------------------------------------------

int insert_global_stats(const char *id, const char *title, 
			const char *value, sqlite3 *db);

//------------------------------------------------------------------------

int insert_chunk(const char *chr, int chunk, int counter, sqlite3 *db);
int inc_chunk(const char *chr, int chunk, sqlite3 *db);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // end of DB_UTILS_H
