#include "vcf_file.h"


//====================================================================================
//  vcf_file.c
//  vcf file management functions
//====================================================================================

//-----------------------------------------------------
// vcf_open
//-----------------------------------------------------


vcf_file_t *vcf_open(char *filename) {
    if (!exists(filename)) {
        return NULL;
    }
    
	vcf_file_t *vcf_file = (vcf_file_t *) malloc(sizeof(vcf_file_t));
    vcf_file->filename = filename;
    
    // Initialize file descriptor or mmap'd buffers
    if (mmap_vcf) {
        size_t len;
        char *data = mmap_file(&len, filename);
//        vcf_file->fd = NULL;
        vcf_file->fd = fopen(filename, "r");
        vcf_file->data = data;
        vcf_file->data_len = len;
    } else {
        vcf_file->fd = fopen(filename, "r");
        vcf_file->data = NULL;
        vcf_file->data_len = 0;
    }

	// Initialize header
    vcf_file->header_entries = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);
	vcf_file->num_header_entries = 0;
	
	// Initialize samples names list
    vcf_file->samples_names = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);
	vcf_file->num_samples = 0;
	
	// Initialize records
    vcf_file->records = array_list_new(100, 1.2, COLLECTION_MODE_SYNCHRONIZED);
	vcf_file->num_records = 0;
	
	return vcf_file;
}


//-----------------------------------------------------
// vcf_close and memory freeing
//-----------------------------------------------------

void vcf_close(vcf_file_t *vcf_file) {
    // Free file format
//     free(vcf_file->format);

    // Free samples names
    array_list_free(vcf_file->samples_names, free);
    // Free header entries
    array_list_free(vcf_file->header_entries, vcf_header_entry_free);
    // Free records list
    array_list_free(vcf_file->records, vcf_record_free);

    if (mmap_file) {
        munmap((void*) vcf_file->data, vcf_file->data_len);
    } 
//    else {
        fclose(vcf_file->fd);
//    }
    free(vcf_file);
}


//-----------------------------------------------------
// Header/records allocation/freeing
//-----------------------------------------------------


vcf_header_entry_t* create_header_entry() {
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

void vcf_header_entry_free(vcf_header_entry_t *vcf_header_entry) {
	// Free entry name
	free(vcf_header_entry->name);
	// Free list of keys
	list_item_t* item = NULL;
	while ( (item = list_remove_item_async(vcf_header_entry->keys)) != NULL )  {
		free(item->data_p);
		list_item_free(item);
	}
	free(vcf_header_entry->keys);
	// Free list of values
	item = NULL;
	while ( (item = list_remove_item_async(vcf_header_entry->values)) != NULL ) {
		free(item->data_p);
		list_item_free(item);
	}
	free(vcf_header_entry->values);
	
	free(vcf_header_entry);
}

vcf_record_t* create_record() {
    vcf_record_t *record = (vcf_record_t*) calloc (1, sizeof(vcf_record_t));
    record->samples = array_list_new(16, 1.5, COLLECTION_MODE_ASYNCHRONIZED);
    return record;
}

void vcf_record_free(vcf_record_t *vcf_record) {
//     free(vcf_record->chromosome);
//     free(vcf_record->id);
//     free(vcf_record->reference);
//     free(vcf_record->alternate);
//     free(vcf_record->filter);
//     free(vcf_record->info);
//     free(vcf_record->format);
    array_list_free(vcf_record->samples, free);
    free(vcf_record);
}


//-----------------------------------------------------
// I/O operations (read and write) in various ways
//-----------------------------------------------------

// int vcf_read(vcf_file_t *vcf_file) {
// 	return vcf_ragel_read(NULL, 1, vcf_file, 0);
// }

int vcf_read_batches(list_t *batches_list, size_t batch_size, vcf_file_t *vcf_file) {
    return vcf_light_read(batches_list, batch_size, vcf_file);
}

int vcf_multiread_batches(list_t **batches_list, size_t batch_size, vcf_file_t **vcf_files, int num_files) {
    return vcf_light_multiread(batches_list, batch_size, vcf_files, num_files);
}

int vcf_parse_batches(list_t *batches_list, size_t batch_size, vcf_file_t *vcf_file, int read_samples) {
	return vcf_read_and_parse(batches_list, batch_size, vcf_file, read_samples);
}

int vcf_write(vcf_file_t *vcf_file, char *filename) {
	FILE *fd = fopen(filename, "w");
	if (fd < 0) {
		fprintf(stderr, "Error opening file: %s\n", filename);
		exit(1);
	}
	
	if (vcf_write_to_file(vcf_file, fd)) {
		fprintf(stderr, "Error writing file: %s\n", filename);
		fclose(fd);
		exit(1);
	}
	
	fclose(fd);
	return 0;
}

//-----------------------------------------------------
// load data into the vcf_file_t
//-----------------------------------------------------

int add_header_entry(vcf_header_entry_t *header_entry, vcf_file_t *vcf_file) {
    int result = array_list_insert(header_entry, vcf_file->header_entries);
    if (result) {
        (vcf_file->num_header_entries)++;
        LOG_DEBUG_F("header entry %zu\n", vcf_file->header_entries->size);
    } else {
        LOG_DEBUG_F("header entry %zu not inserted\n", vcf_file->num_header_entries);
    }
    return result;
}

int add_sample_name(char *name, int length, vcf_file_t *vcf_file) {
    int result = array_list_insert(strndup(name, length), vcf_file->samples_names);
    if (result) {
        (vcf_file->num_samples)++;
        LOG_DEBUG_F("sample %zu is %s\n", vcf_file->samples_names->size, name);
    } else {
        LOG_DEBUG_F("sample %zu not inserted\n", vcf_file->num_samples);
    }
    return result;
}

int add_record(vcf_record_t* record, vcf_file_t *vcf_file) {
    int result = array_list_insert(record, vcf_file->records);
    if (result) {
        (vcf_file->num_records)++;
        LOG_DEBUG_F("record %zu\n", vcf_file->num_records);
    } else {
        LOG_DEBUG_F("record %zu not inserted\n", vcf_file->num_records);
    }
    return result;
}
