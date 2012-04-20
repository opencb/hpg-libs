#ifndef GFF_READ_H
#define GFF_READ_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <log.h>

#include "gff_file_structure.h"

//====================================================================================
//  gff_read.h
//
//  gff reading functions prototypes
//====================================================================================



/* ************ Header management functions **********************/

gff_header_entry_t* create_gff_header_entry();

void set_gff_header_entry_text(char *text, gff_header_entry_t *entry);

/* ************ Record management functions **********************/

gff_record_t* create_gff_record();

void set_gff_record_sequence(char* sequence, gff_record_t* gff_record);

void set_gff_record_source(char* source, gff_record_t* gff_record);

void set_gff_record_feature(char* feature, gff_record_t* gff_record);

void set_gff_record_start(long start, gff_record_t* gff_record);

void set_gff_record_end(long end, gff_record_t* gff_record);

void set_gff_record_score(int score, gff_record_t* gff_record);

void set_gff_record_strand(char strand, gff_record_t* gff_record);

void set_gff_record_frame(int frame, gff_record_t* gff_record);

void set_gff_record_attribute(char* attribute, gff_record_t* gff_record);

#endif
