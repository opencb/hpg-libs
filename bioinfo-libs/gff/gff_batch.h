#ifndef GFF_BATCH_H
#define GFF_BATCH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <list.h>

#include "gff_file_structure.h"

//====================================================================================
//  gff_batch.h
//
//  gff_batch_t structures and prototypes
//====================================================================================

/**
 * Struct which represents a batch of GFF records whose fields have been 
 * already loaded.
 */
typedef list_t gff_batch_t;
// typedef struct gff_batch {
//        size_t num_records;	// Number of records read
//        size_t size;             // Max buffer size (can be greater than num_records)
//        gff_record_t **records;	// Records read
// } gff_batch_t;

gff_batch_t* gff_batch_new(size_t size);

void gff_batch_free(gff_batch_t *gff_batch);

void add_record_to_batch(gff_record_t *record, gff_batch_t *gff_batch);

int batch_is_empty(gff_batch_t *gff_batch);

int batch_is_full(gff_batch_t *gff_batch);

void gff_batch_print(FILE *fd, gff_batch_t *gff_batch);

#endif
