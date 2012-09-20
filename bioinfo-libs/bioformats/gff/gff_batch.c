#include "gff_batch.h"
#include "gff_file.h"

gff_batch_t* gff_batch_new(size_t size) 
{
    gff_batch_t *gff_batch_p = (gff_batch_t*) malloc (sizeof(gff_batch_t));
    list_init("batch", 1, size, gff_batch_p);
    return gff_batch_p;
}

void gff_batch_free(gff_batch_t* gff_batch_p) 
{
    list_item_t *item;
    while ((item = list_remove_item_async(gff_batch_p)) != NULL)
    {
        gff_record_free(item->data_p);
        list_item_free(item);
    }
    free(gff_batch_p);
}

void add_record_to_gff_batch(gff_record_t *record, gff_batch_t *gff_batch)
{
    list_item_t *item = list_item_new(gff_batch->length, 1, record);
    list_insert_item(item, gff_batch);
}

inline int gff_batch_is_empty(gff_batch_t *gff_batch)
{
    return gff_batch->length == 0;
}

inline int gff_batch_is_full(gff_batch_t *gff_batch)
{
    return gff_batch->length == gff_batch->max_length;
}

void gff_batch_print(FILE *fd, gff_batch_t *gff_batch)
{
    gff_record_t *first_record = gff_batch->first_p->data_p;
    fprintf(fd, "Batch with %zu/%zu records - %s in %ld\n", 
        gff_batch->length, gff_batch->max_length, 
        first_record->sequence, first_record->start);
}
