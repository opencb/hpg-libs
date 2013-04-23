/*
 * bam_db.c
 *
 *  Created on: Apr 23, 2013
 *      Author: jtarraga
 */

#include "bam_db.h"

//------------------------------------------------------------------------

bam_query_fields_t *bam_query_fields_new(char *chr, int strand, int start,
					 int end, int flag) {
  bam_query_fields_t *p = calloc(1, sizeof(bam_query_fields_t));

  p->chr = strdup(chr);
  p->strand = strand;
  p->start = start;
  p->end = end;
  p->flag = flag;

  return p;
}

void bam_query_fields_free(bam_query_fields_t *p) {
  if (p) {
    if (p->chr) free(p->chr);
    free(p);
  }
}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int create_bam_query_fields(sqlite3 *db) {

  // create record_stats table for bam files
  int rc;
  char *error_msg;
  char sql[1000];
  sprintf(sql, "CREATE TABLE record_query_fields (chr TEXT, strand CHAR, start INT, end INT, flag CHAR)");

  if (rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg)) {
    LOG_FATAL_F("Stats database failed: %s\n", error_msg);
  }
  return 0;
}

//------------------------------------------------------------------------

int insert_bam_query_fields(void *custom_fields, sqlite3 *db) {
  bam_query_fields_t *fields = (bam_query_fields_t *) custom_fields;

  // insert into record_query_fields
  int rc;
  char sql[1024];
  sprintf(sql, "INSERT INTO record_query_fields VALUES('%s', %i, %i, %i, %i)", 
	  fields->chr, fields->strand, fields->start, fields->end, fields->flag);

  char *error_msg;
  if (rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg)) {
    LOG_DEBUG_F("Stats database failed: %s\n", error_msg);
  }
  return rc;
}


//------------------------------------------------------------------------
//------------------------------------------------------------------------

