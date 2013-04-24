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

inline int exec_sql(char *sql, sqlite3* db) {
  int rc;
  char *error_msg;
  if (rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg)) {
    LOG_DEBUG_F("Stats database failed (%s): %s\n", sql, error_msg);
  }
  return rc;
}

//------------------------------------------------------------------------
//             C R E A T I O N        F U N C T I O N S
//------------------------------------------------------------------------

int create_stats_db(const char *db_name, int chunksize, 
		    int (*create_custom_fields)(sqlite3 *), sqlite3** db) {

  // create sqlite db
  if (sqlite3_open(db_name, db)) {
    LOG_FATAL_F("Could not open stats database (%s): %s\n", 
		db_name, sqlite3_errmsg(*db));
  }

  int rc;
  char sql[128];
  
  sprintf(sql, "BEGIN TRANSACTION");
  rc = exec_sql(sql, *db);

  // create global stats table and index, and insert the chunksize
  sprintf(sql, "CREATE TABLE global_stats (name TEXT PRIMARY KEY, title TEXT, value TEXT)");
  rc = exec_sql(sql, *db);


  sprintf(sql, "%i", chunksize);
  rc = insert_global_stats("CHUNK_SIZE", "Chunk size", sql, *db);
  rc = insert_global_stats("CHR_PREFIX", "Chromosome prefix", "", *db);

  // create chunks table
  sprintf(sql, "CREATE TABLE chunk (chromosome TEXT, chunk_id INT, start INT, end INT, features_count INT)");
  rc = exec_sql(sql, *db);

  // create record_query_fields table for bam, vcf.. files
  if (create_custom_fields) {
    rc = create_custom_fields(*db);
  }
  
  sprintf(sql, "END TRANSACTION");
  rc = exec_sql(sql, *db);

  return rc;
}

//------------------------------------------------------------------------

int create_stats_index(int (*create_custom_index)(sqlite3 *), sqlite3* db) {
  //  sprintf(sql, "CREATE INDEX id_idx ON global_stats (id)");
  //  rc = exec_sql(sql, *db);

  int rc;
  char sql[128];

  // create chunks index
  sprintf(sql, "CREATE INDEX chunk_chromosome_chunk_id_idx ON chunk (chromosome, chunk_id)");
  rc = exec_sql(sql, db);

  // create custeom index (for the record_query_fields table)
  if (create_custom_index) {
    rc = create_custom_index(db);
  }

  return rc;
}

//------------------------------------------------------------------------
//             G L O B A L     S T A T S     T A B L E
//------------------------------------------------------------------------

int insert_global_stats(const char *name, const char *title, 
			const char *value, sqlite3 *db) {
  int rc;

  // IMPORTANT !!!
  //
  // The sqlite3_exec() interface is a convenience wrapper around 
  // sqlite3_prepare_v2(), sqlite3_step(), and sqlite3_finalize(), 
  // that allows an application to run multiple statements of SQL 
  // without having to use a lot of C code.


  char sql[strlen(name) + strlen(title) + strlen(value) + 100];
  sprintf(sql, "INSERT INTO global_stats VALUES('%s', '%s', '%s')", 
	  name, title, value);
  rc = exec_sql(sql, db);

/*
  sqlite3_stmt* insert_stmt;
  char sql[] = "INSERT INTO global_stats VALUES(?1, ?2, ?3)";

  sqlite3_prepare_v2(db, sql, strlen(sql) + strlen(name) + strlen(title) + strlen(value), &insert_stmt, NULL);
  sqlite3_bind_text(insert_stmt, 1, name, strlen(name), SQLITE_STATIC);
  sqlite3_bind_text(insert_stmt, 2, title, strlen(title), SQLITE_STATIC);
  sqlite3_bind_text(insert_stmt, 3, value, strlen(value), SQLITE_STATIC);

  if (rc = sqlite3_step(insert_stmt) != SQLITE_DONE) {
    LOG_DEBUG_F("Stats databases failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
  } else {
    rc = SQLITE_OK;
  }
  sqlite3_finalize(insert_stmt);
*/

  return rc;

}

//------------------------------------------------------------------------

int prepare_statement_global_stats(sqlite3 *db, sqlite3_stmt **stmt) {
  char sql[] = "INSERT INTO global_stats VALUES(?1, ?2, ?3)"; 
  return sqlite3_prepare_v2(db, sql, strlen(sql) + 1024, stmt, NULL);
}

//------------------------------------------------------------------------

int finalize_statement_global_stats(sqlite3_stmt *stmt) {
  return sqlite3_finalize(stmt);
}

//------------------------------------------------------------------------

int insert_statement_global_stats(const char *name, const char *title, 
				  const char *value, sqlite3_stmt *stmt, 
				  sqlite3 *db) {  
  int rc;

  sqlite3_bind_text(stmt, 1, name, strlen(name), SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 2, title, strlen(title), SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 3, value, strlen(value), SQLITE_TRANSIENT);

  if (rc = sqlite3_step(stmt) != SQLITE_DONE) {
    LOG_DEBUG_F("Stats databases failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
  } else {
    rc = SQLITE_OK;
  }
  sqlite3_reset(stmt);

  return rc;
}

//------------------------------------------------------------------------
//                     C H U N K       T A B L E
//------------------------------------------------------------------------

int insert_chunk(const char *chr, int chunk_id, int start, int end, 
		 int features_count, sqlite3 *db) {

  int rc;
  sqlite3_stmt* stmt;
  char sql[] = "INSERT INTO chunk VALUES(?1, ?2, ?3, ?4, ?5)";

  sqlite3_prepare_v2(db, sql, strlen(sql) + strlen(chr), &stmt, NULL);
  sqlite3_bind_text(stmt, 1, chr, strlen(chr), SQLITE_STATIC);
  sqlite3_bind_int(stmt, 2, chunk_id);
  sqlite3_bind_int(stmt, 3, start);
  sqlite3_bind_int(stmt, 4, end);
  sqlite3_bind_int(stmt, 5, features_count);

  if (rc = sqlite3_step(stmt) != SQLITE_DONE) {
    LOG_DEBUG_F("Stats databases failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
  } else {
    rc = SQLITE_OK;
  }
  sqlite3_finalize(stmt);

  return rc;
}

//------------------------------------------------------------------------

int inc_chunk(const char *chr, int chunk_id, int chunk_start, int chunk_end, sqlite3 *db) {
  int rc;
  sqlite3_stmt *stmt;
  char sql[strlen(chr) + 200];

  strcpy(sql, "SELECT features_count FROM chunk WHERE chromosome = ?1 AND chunk_id = ?2");

  sqlite3_prepare_v2(db, sql, strlen(sql) + strlen(chr), &stmt, NULL);
  sqlite3_bind_text(stmt, 1, chr, strlen(chr), SQLITE_STATIC);
  sqlite3_bind_int(stmt, 2, chunk_id);
  
  rc = sqlite3_step(stmt);

  if (rc == SQLITE_ROW) {
    sprintf(sql, "UPDATE chunk SET features_count = features_count + 1 WHERE chromosome = '%s' AND chunk_id = %i", 
	    chr, chunk_id);
    rc = exec_sql(sql, db);
  } else if (rc == SQLITE_DONE) {
    sprintf(sql, "INSERT INTO chunk VALUES('%s', %i, %i, %i, 1)", chr, chunk_id,
	    chunk_start, chunk_end);
    rc = exec_sql(sql, db);
  }

  sqlite3_finalize(stmt);
}

//------------------------------------------------------------------------
//         R E C O R D     Q U E R Y     F I E L D S    T A B L E
//------------------------------------------------------------------------

int insert_record_query_fields(const char *chr, int chr_length, int chunksize, 
			       int start, int end, void *fields,
			       int (*insert_custom_fields)(void *, sqlite3 *), 
			       sqlite3* db) {
  int rc;

  if (rc = insert_custom_fields(fields, db)) {
    LOG_DEBUG_F("Stats database failed: %s\n", sqlite3_errmsg(db));
  } else {
    /*
    // if custom fields insertion succeeded then
    // update features counter for "touched" chunks
    int chunk_start, chunk_end;
    int chunk_id_start = start / chunksize;
    int chunk_id_end = end / chunksize;

    for (int chunk_id = chunk_id_start; chunk_id <= chunk_id_end; chunk_id++) {
      chunk_start = chunk_id * chunksize;
      chunk_end = chunk_start + chunksize - 1;
      if (chunk_end > chr_length) {
	chunk_end = chr_length - 1;
      }

      rc += inc_chunk(chr, chunk_id, chunk_start, chunk_end, db);
    }
    */
  }

  return rc;
}

//------------------------------------------------------------------------

int insert_statement_record_query_fields(const char *chr, int chr_length, int chunksize, 
					 int start, int end, void *fields,
					 int (*insert_statement_custom_fields)(void *, sqlite3_stmt *, sqlite3 *), 
					 sqlite3_stmt *custom_stmt,
					 sqlite3* db) {
  int rc;


  if (rc = insert_statement_custom_fields(fields, custom_stmt, db)) {
    LOG_DEBUG_F("Stats database failed: %s\n", sqlite3_errmsg(db));
  } else {
    /*
    // if custom fields insertion succeeded then
    // update features counter for "touched" chunks
    int chunk_start, chunk_end;
    int chunk_id_start = start / chunksize;
    int chunk_id_end = end / chunksize;

    for (int chunk_id = chunk_id_start; chunk_id <= chunk_id_end; chunk_id++) {
      chunk_start = chunk_id * chunksize;
      chunk_end = chunk_start + chunksize - 1;
      if (chunk_end > chr_length) {
	chunk_end = chr_length - 1;
      }

      rc += inc_chunk(chr, chunk_id, chunk_start, chunk_end, db);
    }
    */
  }

  return rc;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
