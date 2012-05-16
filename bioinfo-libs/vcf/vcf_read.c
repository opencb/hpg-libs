#include "vcf_read.h"

//====================================================================================
// vcf reading functions
//====================================================================================


/* ************ Header management functions **********************/

void set_file_format(char *fileformat, vcf_file_t *vcf_file)
{
    vcf_file->format = fileformat;
    LOG_DEBUG_F("set format = %s\n", vcf_file->format);
}

void set_header_entry_name(char *name, vcf_header_entry_t *entry)
{
    entry->name = name;
    LOG_DEBUG_F("set name: %s\n", entry->name);
}

void add_header_entry_key(char *key, vcf_header_entry_t *entry)
{
    list_item_t *item = list_item_new(entry->num_keys, 1, key);
    int result = list_insert_item(item, entry->keys);
    if (result) {
        entry->num_keys++;
        LOG_DEBUG_F("key %zu = %s\n", entry->num_keys, (char*) item->data_p);
    } else {
        LOG_DEBUG_F("key %zu not inserted\n", entry->num_keys);
    }
}

void add_header_entry_value(char *value, vcf_header_entry_t *entry)
{
    list_item_t *item = list_item_new(entry->num_values, 1, value);
    int result = list_insert_item(item, entry->values);
    if (result) {
        entry->num_values++;
        LOG_DEBUG_F("value %zu = %s\n", entry->num_values, (char*) item->data_p);
    } else {
        LOG_DEBUG_F("value %zu not inserted\n", entry->num_values);
    }
}


/* ************ Record management functions **********************/

void set_record_chromosome(char* chromosome, vcf_record_t* vcf_record)
{
    vcf_record->chromosome = chromosome;
    LOG_DEBUG_F("set chromosome: %s\n", vcf_record->chromosome);
}

void set_record_position(long position, vcf_record_t* vcf_record) 
{
    vcf_record->position = position;
    LOG_DEBUG_F("set position: %ld\n", vcf_record->position);
}

void set_record_id(char* id, vcf_record_t* vcf_record) 
{
    vcf_record->id = id;
    LOG_DEBUG_F("set id: %s\n", vcf_record->id);
}

void set_record_reference(char* reference, vcf_record_t* vcf_record) 
{
    vcf_record->reference = reference;
    LOG_DEBUG_F("set reference: %s\n", vcf_record->reference);
}

void set_record_alternate(char* alternate, vcf_record_t* vcf_record)
{
    vcf_record->alternate = alternate;
    LOG_DEBUG_F("set alternate: %s\n", vcf_record->alternate);
}

void set_record_quality(float quality, vcf_record_t* vcf_record)
{
    vcf_record->quality = quality;
    LOG_DEBUG_F("set quality: %f\n", vcf_record->quality);
}

void set_record_filter(char* filter, vcf_record_t* vcf_record)
{
    vcf_record->filter = filter;
    LOG_DEBUG_F("set filter: %s\n", vcf_record->filter);
}

void set_record_info(char* info, vcf_record_t* vcf_record)
{
    vcf_record->info = info;
    LOG_DEBUG_F("set info: %s\n", vcf_record->info);
}

void set_record_format(char* format, vcf_record_t* vcf_record)
{
    vcf_record->format = format;
    LOG_DEBUG_F("set format: %s\n", vcf_record->format);
}

void add_record_sample(char* sample, vcf_record_t* vcf_record, size_t *sample_idx)
{
    list_item_t *item = list_item_new(*sample_idx, 1, sample);
    int result = list_insert_item(item, vcf_record->samples);
    if (result) {
        (*sample_idx)++;
        LOG_DEBUG_F("sample %zu = %s\n", *sample_idx, (char*) item->data_p);
    } else {
        LOG_DEBUG_F("sample %zu not inserted\n", *sample_idx);
    }
}
