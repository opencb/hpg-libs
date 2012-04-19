#include "gff_read.h"

//====================================================================================
// gff reading functions
//====================================================================================


/* ************ Header management functions **********************/

gff_header_entry_t* create_header_entry() {
    gff_header_entry_t *entry = (gff_header_entry_t*) malloc (sizeof(gff_header_entry_t));
    entry->text = NULL;
    return entry;
}

void set_header_entry_text(char *text, gff_header_entry_t *entry) {
    entry->text = text;
    dprintf("set text: %s\n", entry->text);
}


/* ************ Record management functions **********************/

gff_record_t* create_record() {
    return malloc (sizeof(gff_record_t));
}

void set_record_sequence(char* sequence, gff_record_t* gff_record) {
    gff_record->sequence = sequence;
    dprintf("set sequence: %s\n", gff_record->sequence);
}

void set_record_source(char* source, gff_record_t* gff_record) {
    gff_record->source = source;
    dprintf("set source: %s\n", gff_record->source);
}

void set_record_feature(char* feature, gff_record_t* gff_record) {
    gff_record->feature = feature;
    dprintf("set feature: %s\n", gff_record->feature);
}

void set_record_start(long start, gff_record_t* gff_record) {
    gff_record->start = start;
    dprintf("set start: %ld\n", gff_record->start);
}

void set_record_end(long end, gff_record_t* gff_record) {
    gff_record->end = end;
    dprintf("set end: %ld\n", gff_record->end);
}

void set_record_score(int score, gff_record_t* gff_record) {
    gff_record->score = score;
    dprintf("set score: %d\n", gff_record->score);
}

void set_record_strand(char strand, gff_record_t* gff_record) {
    gff_record->strand = strand;
    dprintf("set strand: %c\n", gff_record->strand);
}

void set_record_frame(int frame, gff_record_t* gff_record) {
    gff_record->frame = frame;
    dprintf("set frame: %d\n", gff_record->frame);
}

void set_record_attribute(char* attribute, gff_record_t* gff_record) {
    gff_record->attribute = attribute;
    dprintf("set attribute: %s\n", gff_record->attribute);
}
