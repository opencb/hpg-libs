#include "vcf_batch.h"
#include "vcf_file.h"

vcf_batch_t* vcf_batch_new(size_t size) 
{
	vcf_batch_t *vcf_batch_p = (vcf_batch_t*) malloc (sizeof(vcf_batch_t));
	list_init("batch", 1, size, vcf_batch_p);
	
// 	vcf_batch_p->records = (vcf_record_t**) malloc (size * sizeof(vcf_record_t*));
// 	memset(vcf_batch_p->records, 0, size * sizeof(vcf_record_t*));
// 	
// 	vcf_batch_p->num_records = 0;
// 	vcf_batch_p->size = size;
	
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
	
// 	if (vcf_batch_p->records != NULL)
// 	{
// 		for (int i = 0; i < vcf_batch_p->num_records; i++)
// 		{
// 			dprintf("About to free record %d\n", i);
// 			vcf_record_free(vcf_batch_p->records[i]);
// 		}
// 		free(vcf_batch_p->records);
// 	}
	free(vcf_batch_p);
}

void add_record_to_batch(vcf_record_t *record, vcf_batch_t *vcf_batch)
{
	list_item_t *item = list_item_new(vcf_batch->length, 1, record);
	list_insert_item(item, vcf_batch);
// 	vcf_batch->records[vcf_batch->num_records] = record;
// 	vcf_batch->num_records++;
}

inline int batch_is_empty(vcf_batch_t *vcf_batch)
{
// 	return vcf_batch->num_records == 0;
	return vcf_batch->length == 0;
}

inline int batch_is_full(vcf_batch_t *vcf_batch)
{
// 	return vcf_batch->num_records == vcf_batch->size;
	return vcf_batch->length == vcf_batch->max_length;
}

void vcf_batch_print(FILE *fd, vcf_batch_t *vcf_batch)
{
	vcf_record_t *first_record = vcf_batch->first_p->data_p;
	fprintf(fd, "Batch with %zu/%zu records - %s in %ld\n", 
		vcf_batch->length, vcf_batch->max_length, 
		first_record->chromosome, first_record->position);
}
