#ifndef BAM_DB_H
#define BAM_DB_H

/*
 * bam_db.h
 *
 *  Created on: Apr 23, 2013
 *      Author: jtarraga
 */

//------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#include "commons/log.h"
#include "commons/sqlite/sqlite3.h"


//------------------------------------------------------------------------

typedef struct bam_query_fields {
  char *chr;
  int strand;
  int start;
  int end;
  int flag;
} bam_query_fields_t;

//------------------------------------------------------------------------

bam_query_fields_t *bam_query_fields_new(char *chr, int strand, int start,
					 int end, int flag);

void bam_query_fields_free(bam_query_fields_t *p);

//------------------------------------------------------------------------
// 
//------------------------------------------------------------------------

int create_bam_query_fields(sqlite3 *db);

//------------------------------------------------------------------------

int insert_bam_query_fields(void *custom_fields, sqlite3 *db);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // end of BAM_DB_H
