#include "gff_write.h"

int gff_write_to_file(gff_file_t *gff_file, FILE *fd) {
    if(gff_file == NULL || fd == NULL) {
        return -1;
    }

    // Write header entries
    cp_list_iterator *headers_iter = cp_list_create_iterator(gff_file->header_entries, COLLECTION_LOCK_READ);
//     cp_list_iterator_init(headers_iter, gff_file->header_entries, COLLECTION_LOCK_READ);
    gff_header_entry_t *entry = NULL;
    while ((entry = cp_list_iterator_next(headers_iter)) != NULL) {
        write_gff_header_entry(entry, fd);
    }
    
    // Write records
    cp_list_iterator *records_iter = cp_list_create_iterator(gff_file->records, COLLECTION_LOCK_READ);
//     cp_list_iterator_init(records_iter, gff_file->records, COLLECTION_LOCK_READ);
    gff_record_t *record = NULL;
    while ((record = cp_list_iterator_next(records_iter)) != NULL) {
        write_gff_record(record, fd);
    }

    return 0;
}

void write_gff_header_entry(gff_header_entry_t *entry, FILE *fd) {
    if (entry == NULL || fd == NULL) {
        return;
    }

    fprintf(fd, "##%s\n", entry->text);
}

void write_gff_batch(gff_batch_t *gff_batch, FILE *fd) {
    if(gff_batch == NULL || fd == NULL) {
        return;
    }
    for (list_item_t *i = gff_batch->first_p; i != NULL; i = i->next_p) {
        write_gff_record(i->data_p, fd);
    }
}

void write_gff_record(gff_record_t* gff_record, FILE *fd) {
    if(gff_record == NULL || fd == NULL) {
        return;
    }
    
    fprintf(fd, "%s\t%s\t%s\t%ld\t%ld\t", gff_record->sequence, gff_record->source, gff_record->feature, gff_record->start, gff_record->end);
    if (gff_record->score < 0) {
        fprintf(fd, ".\t");
    } else {
        fprintf(fd, "%d\t", gff_record->score);
    }
    fprintf(fd, "%c\t", gff_record->strand);
    if (gff_record->frame < 0) {
        fprintf(fd, ".\t");
    } else {
        fprintf(fd, "%d\t", gff_record->frame);
    }
    fprintf(fd, "%s\t\n", gff_record->attribute);
}
