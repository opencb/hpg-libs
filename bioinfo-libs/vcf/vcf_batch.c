#include "vcf_batch.h"
#include "vcf_file.h"

vcf_batch_t* vcf_batch_new(size_t size) {
    vcf_batch_t *vcf_batch = calloc (1, sizeof(vcf_batch_t));
    vcf_batch->text = NULL;
    vcf_batch->records = array_list_new(size, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    
    return vcf_batch;
}

void vcf_batch_free(vcf_batch_t* vcf_batch) {
    if (!vcf_batch) {
        return;
    }
    
    if (vcf_batch->text && !mmap_vcf) {
//         printf("text to free = '%.*s'\n", 50, vcf_batch->text);
        free(vcf_batch->text);
    }
    array_list_free(vcf_batch->records, vcf_record_free);
    free(vcf_batch);
}

void add_record_to_vcf_batch(vcf_record_t *record, vcf_batch_t *vcf_batch) {
    array_list_insert(record, vcf_batch->records);
}

inline int vcf_batch_is_empty(vcf_batch_t *vcf_batch) {
    return vcf_batch->records->size == 0;
}

inline int vcf_batch_is_full(vcf_batch_t *vcf_batch) {
//     printf("batch size = %zu\tcapacity = %zu\n", vcf_batch->size, vcf_batch->capacity);
    return vcf_batch->records->size == vcf_batch->records->capacity;
}

void vcf_batch_print(FILE *fd, vcf_batch_t *vcf_batch) {
    vcf_record_t *first_record = array_list_get(0, vcf_batch->records);
    LOG_DEBUG_F("Batch with %zu/%zu records - %s in %ld\n", 
                vcf_batch->records->size, vcf_batch->records->capacity, 
                first_record->chromosome, first_record->position);
}
