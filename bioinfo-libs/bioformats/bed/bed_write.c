#include "bed_write.h"

int bed_write_to_file(bed_file_t *bed_file, FILE *fd) {
    assert(bed_file);
    assert(fd);

    // Write header entries
    linked_list_iterator_t *headers_iter = linked_list_iterator_new(bed_file->header_entries);
    bed_header_entry_t *entry = NULL;
    while ((entry = linked_list_iterator_curr(headers_iter))) {
        write_bed_header_entry(entry, fd);
        linked_list_iterator_next(headers_iter);
    }
    linked_list_iterator_free(headers_iter);
    
    // Write records
    linked_list_iterator_t *records_iter = linked_list_iterator_new(bed_file->records);
    bed_record_t *record = NULL;
    while ((entry = linked_list_iterator_curr(records_iter))) {
        write_bed_record(record, fd);
        linked_list_iterator_next(records_iter);
    }
    linked_list_iterator_free(records_iter);

    return 0;
}

void write_bed_header_entry(bed_header_entry_t *entry, FILE *fd) {
    assert(entry);
    assert(fd);

    fprintf(fd, "##%s\n", entry->text);
}

void write_bed_batch(bed_batch_t *bed_batch, FILE *fd) {
    assert(bed_batch);
    assert(fd);
    
    for (int i = 0; i < bed_batch->records->size; i++) {
        write_bed_record(bed_batch->records->items[i], fd);
    }
}

void write_bed_record(bed_record_t* bed_record, FILE *fd) {
    assert(bed_record);
    assert(fd);
    
    fprintf(fd, "%s\t%s\t%s\t%ld\t%ld\t", bed_record->sequence, bed_record->source, bed_record->feature, bed_record->start, bed_record->end);
    if (bed_record->score < 0) {
        fprintf(fd, ".\t");
    } else {
        fprintf(fd, "%d\t", bed_record->score);
    }
    fprintf(fd, "%c\t", bed_record->strand);
    if (bed_record->frame < 0) {
        fprintf(fd, ".\t");
    } else {
        fprintf(fd, "%d\t", bed_record->frame);
    }
    fprintf(fd, "%s\t\n", bed_record->attribute);
}
