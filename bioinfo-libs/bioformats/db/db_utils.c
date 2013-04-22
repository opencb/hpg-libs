/*
 * db_utils.c
 *
 *  Created on: Apr 22, 2013
 *      Author: jtarraga
 */

#include <stdio.h>
#include <stdlib.h>

#include "commons/log.h"

#include "db_utils.h"

//------------------------------------------------------------------------

inline void exec_sql(char *sql, sqlite3* db) {
  int rc;
  char *error_msg;
  if (rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg)) {
    LOG_FATAL_F("Stats database failed: %s\n", error_msg);
  }
}

//------------------------------------------------------------------------
//             C R E A T I O N        F U N C T I O N S
//------------------------------------------------------------------------

int create_stats_db(const char *db_name, int chunksize, 
		    int (*format_table_fn)(sqlite3 *), sqlite3** db) {

  // create sqlite db
  if (sqlite3_open(db_name, db)) {
    LOG_FATAL_F("Could not open stats database (%s): %s\n", 
		db_name, sqlite3_errmsg(*db));
  }

  int rc;
  char sql[128];
  
  sprintf(sql, "BEGIN TRANSACTION");
  exec_sql(sql, *db);

  // create global stats table and index
  sprintf(sql, "CREATE TABLE global_stats (id VARCHAR(40) PRIMARY KEY, title VARCHAR(128), value VARCHAR(40))");
  exec_sql(sql, *db);

  sprintf(sql, "CREATE INDEX id_idx ON global_stats (id)");
  exec_sql(sql, *db);

  // create chunks table and index
  sprintf(sql, "CREATE TABLE chunks (chr VARCHAR(40), chunk INTEGER, counter INTEGER)");
  exec_sql(sql, *db);

  sprintf(sql, "CREATE INDEX chr_chunk_idx ON chunks (chr, chunk)");
  exec_sql(sql, *db);

  // create format table: for bam, vcf.. files
  rc = format_table_fn(*db);
  
  sprintf(sql, "END TRANSACTION");
  exec_sql(sql, *db);

  return rc;
}

//------------------------------------------------------------------------
//             G L O B A L     S T A T S     T A B L E
//------------------------------------------------------------------------

int insert_global_stats(const char *id, const char *title, 
			const char *value, sqlite3 *db) {

  char sql[strlen(id) + strlen(title) + strlen(value) + 200];

  sprintf(sql, "INSERT INTO global_stats VALUES('%s', '%s', '%s')", 
	  id, title, value);
  exec_sql(sql, db);
}

//------------------------------------------------------------------------
//                     C H U N K       T A B L E
//------------------------------------------------------------------------


int insert_chunk(const char *chr, int chunk, int counter, sqlite3 *db) {

  char sql[strlen(chr) + 200];

  sprintf(sql, "INSERT INTO chunks VALUES('%s', %i, %i)", chr, chunk, counter);
  exec_sql(sql, db);
}

//------------------------------------------------------------------------

int inc_chunk(const char *chr, int chunk, sqlite3 *db) {
  int rc;
  sqlite3_stmt *stmt;
  char sql[strlen(chr) + 200];

  sprintf(sql, "SELECT counter FROM chunks WHERE chr = '%s' AND chunk = %i", 
	  chr, chunk);

  if (sqlite3_prepare_v2(db, sql, -1, &stmt, 0)) {
    LOG_FATAL_F("Stats database failed: %s\n", sqlite3_errmsg(db));
  }

  rc = sqlite3_step(stmt);

  if (rc == SQLITE_ROW) {
    sprintf(sql, "UPDATE chunks SET counter = counter + 1 WHERE chr = '%s' AND chunk = %i", 
	    chr, chunk);
    exec_sql(sql, db);
  } else if (rc == SQLITE_DONE) {
    sprintf(sql, "INSERT INTO chunks VALUES('%s', %i, 1)", chr, chunk);
    exec_sql(sql, db);
  }
}


//------------------------------------------------------------------------
//------------------------------------------------------------------------
