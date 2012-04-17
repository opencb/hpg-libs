#include "vcf_read.h"

//====================================================================================
// vcf reading functions
//====================================================================================


/* ************ Header management functions **********************/

void set_file_format(char *fileformat, vcf_file_t *vcf_file)
{
	vcf_file->format = fileformat;
	dprintf("set format = %s\n", vcf_file->format);
}

vcf_header_entry_t* create_header_entry()
{
	vcf_header_entry_t *entry = (vcf_header_entry_t*) malloc (sizeof(vcf_header_entry_t));
	entry->name = NULL;
    entry->num_keys = 0;
	entry->keys = (list_t*) malloc (sizeof(list_t));
	list_init("keys", 1, INT_MAX, entry->keys);
    entry->num_values = 0;
	entry->values = (list_t*) malloc (sizeof(list_t));
	list_init("values", 1, INT_MAX, entry->values);
	return entry;
}

void set_header_entry_name(char *name, vcf_header_entry_t *entry)
{
	entry->name = name;
	dprintf("set name: %s\n", entry->name);
}

void add_header_entry_key(char *key, vcf_header_entry_t *entry)
{
	list_item_t *item = list_item_new(entry->num_keys, 1, key);
	int result = list_insert_item(item, entry->keys);
	if (result) {
		entry->num_keys++;
		dprintf("key %zu = %s\n", entry->num_keys, (char*) item->data_p);
	} else {
		dprintf("key %zu not inserted\n", entry->num_keys);
	}
}

void add_header_entry_value(char *value, vcf_header_entry_t *entry)
{
	list_item_t *item = list_item_new(entry->num_values, 1, value);
	int result = list_insert_item(item, entry->values);
	if (result) {
		entry->num_values++;
		dprintf("value %zu = %s\n", entry->num_values, (char*) item->data_p);
	} else {
		dprintf("value %zu not inserted\n", entry->num_values);
	}
}


/* ************ Record management functions **********************/

vcf_record_t* create_record()
{
	vcf_record_t *record = (vcf_record_t*) malloc (sizeof(vcf_record_t));
	record->samples = (list_t*) malloc (sizeof(list_t));
	list_init("samples", 1, INT_MAX, record->samples);
	return record;
}

void set_record_chromosome(char* chromosome, vcf_record_t* vcf_record)
{
	vcf_record->chromosome = chromosome;
	dprintf("set chromosome: %s\n", vcf_record->chromosome);
}

void set_record_position(long position, vcf_record_t* vcf_record) 
{
	vcf_record->position = position;
	dprintf("set position: %ld\n", vcf_record->position);
}

void set_record_id(char* id, vcf_record_t* vcf_record) 
{
	vcf_record->id = id;
	dprintf("set id: %s\n", vcf_record->id);
}

void set_record_reference(char* reference, vcf_record_t* vcf_record) 
{
	vcf_record->reference = reference;
	dprintf("set reference: %s\n", vcf_record->reference);
}

void set_record_alternate(char* alternate, vcf_record_t* vcf_record)
{
	vcf_record->alternate = alternate;
	dprintf("set alternate: %s\n", vcf_record->alternate);
}

void set_record_quality(float quality, vcf_record_t* vcf_record)
{
	vcf_record->quality = quality;
	dprintf("set quality: %f\n", vcf_record->quality);
}

void set_record_filter(char* filter, vcf_record_t* vcf_record)
{
	vcf_record->filter = filter;
	dprintf("set filter: %s\n", vcf_record->filter);
}

void set_record_info(char* info, vcf_record_t* vcf_record)
{
	vcf_record->info = info;
	dprintf("set info: %s\n", vcf_record->info);
}

void set_record_format(char* format, vcf_record_t* vcf_record)
{
	vcf_record->format = format;
	dprintf("set format: %s\n", vcf_record->format);
}

void add_record_sample(char* sample, vcf_record_t* vcf_record, size_t *sample_idx)
{
	list_item_t *item = list_item_new(*sample_idx, 1, sample);
	int result = list_insert_item(item, vcf_record->samples);
	if (result) {
		(*sample_idx)++;
		dprintf("sample %zu = %s\n", *sample_idx, (char*) item->data_p);
	} else {
		dprintf("sample %zu not inserted\n", *sample_idx);
	}
}
