#include "gff_read.h"

//====================================================================================
// gff reading functions
//====================================================================================


/* ************ Header management functions **********************/

gff_header_entry_t* create_gff_header_entry() {
    gff_header_entry_t *entry = (gff_header_entry_t*) malloc (sizeof(gff_header_entry_t));
    entry->text = NULL;
    return entry;
}

void set_gff_header_entry_text(char *text, gff_header_entry_t *entry) {
    entry->text = text;
    LOG_DEBUG_F("set text: %s\n", entry->text);
}


/* ************ Record management functions **********************/

gff_record_t* create_gff_record() {
    return malloc (sizeof(gff_record_t));
}

void set_gff_record_sequence(char* sequence, gff_record_t* gff_record) {
    if (starts_with(sequence, "chrom")) {
        gff_record->sequence = (char*) calloc (strlen(sequence)-4, sizeof(char));
        strncat(gff_record->sequence, sequence+5, strlen(sequence)-5);
        free(sequence);
    } else if (starts_with(sequence, "chr")) {
        gff_record->sequence = (char*) calloc (strlen(sequence)-2, sizeof(char));
        strncat(gff_record->sequence, sequence+3, strlen(sequence)-3);
        free(sequence);
    } else {
        gff_record->sequence = sequence;
    }
    LOG_DEBUG_F("set sequence: %s\n", gff_record->sequence);
}

void set_gff_record_source(char* source, gff_record_t* gff_record) {
    gff_record->source = source;
    LOG_DEBUG_F("set source: %s\n", gff_record->source);
}

void set_gff_record_feature(char* feature, gff_record_t* gff_record) {
    gff_record->feature = feature;
    LOG_DEBUG_F("set feature: %s\n", gff_record->feature);
}

void set_gff_record_start(long start, gff_record_t* gff_record) {
    gff_record->start = start;
    LOG_DEBUG_F("set start: %ld\n", gff_record->start);
}

void set_gff_record_end(long end, gff_record_t* gff_record) {
    gff_record->end = end;
    LOG_DEBUG_F("set end: %ld\n", gff_record->end);
}

void set_gff_record_score(int score, gff_record_t* gff_record) {
    gff_record->score = score;
    LOG_DEBUG_F("set score: %d\n", gff_record->score);
}

void set_gff_record_strand(char strand, gff_record_t* gff_record) {
    gff_record->strand = strand;
    LOG_DEBUG_F("set strand: %c\n", gff_record->strand);
}

void set_gff_record_frame(int frame, gff_record_t* gff_record) {
    gff_record->frame = frame;
    LOG_DEBUG_F("set frame: %d\n", gff_record->frame);
}

void set_gff_record_attribute(char* attribute, gff_record_t* gff_record) {
    gff_record->attribute = attribute;
    LOG_DEBUG_F("set attribute: %s\n", gff_record->attribute);
}
