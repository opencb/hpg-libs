#include "vcf_batch.h"
#include "vcf_file.h"

vcf_batch_t* vcf_batch_new(size_t size) 
{
    vcf_batch_t *vcf_batch_p = (vcf_batch_t*) malloc (sizeof(vcf_batch_t));
    list_init("batch", 1, size, vcf_batch_p);
    return vcf_batch_p;
}

void vcf_batch_free(vcf_batch_t* vcf_batch_p) 
{
    // Free records
    list_item_t *item;
    while ((item = list_remove_item_async(vcf_batch_p)) != NULL)
    {
        vcf_record_free(item->data_p);
        list_item_free(item);
    }
    free(vcf_batch_p);
}

void add_record_to_vcf_batch(vcf_record_t *record, vcf_batch_t *vcf_batch)
{
    list_item_t *item = list_item_new(vcf_batch->length, 1, record);
    list_insert_item(item, vcf_batch);
}

inline int vcf_batch_is_empty(vcf_batch_t *vcf_batch)
{
    return vcf_batch->length == 0;
}

inline int vcf_batch_is_full(vcf_batch_t *vcf_batch)
{
    return vcf_batch->length == vcf_batch->max_length;
}

void vcf_batch_print(FILE *fd, vcf_batch_t *vcf_batch)
{
    vcf_record_t *first_record = vcf_batch->first_p->data_p;
    LOG_DEBUG_F("Batch with %zu/%zu records - %s in %ld\n", 
                vcf_batch->length, vcf_batch->max_length, 
                first_record->chromosome, first_record->position);
}
