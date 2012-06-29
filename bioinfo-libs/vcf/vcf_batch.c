#include "vcf_batch.h"
#include "vcf_file.h"

vcf_batch_t* vcf_batch_new(size_t size) {
    return (vcf_batch_t*) array_list_new(size, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
}

void vcf_batch_free(vcf_batch_t* vcf_batch) {
    array_list_free(vcf_batch, vcf_record_free);
}

void add_record_to_vcf_batch(vcf_record_t *record, vcf_batch_t *vcf_batch) {
    array_list_insert(record, vcf_batch);
}

inline int vcf_batch_is_empty(vcf_batch_t *vcf_batch) {
    return vcf_batch->size == 0;
}

inline int vcf_batch_is_full(vcf_batch_t *vcf_batch) {
    return vcf_batch->size == vcf_batch->capacity;
}

void vcf_batch_print(FILE *fd, vcf_batch_t *vcf_batch)
{
    vcf_record_t *first_record = array_list_get(0, vcf_batch);
    LOG_DEBUG_F("Batch with %zu/%zu records - %s in %ld\n", 
                vcf_batch->size, vcf_batch->capacity, 
                first_record->chromosome, first_record->position);
}
