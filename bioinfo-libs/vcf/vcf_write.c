#include "vcf_write.h"

int vcf_write_to_file(vcf_file_t *vcf_file, FILE *fd) {
    if(vcf_file == NULL || fd == NULL) {
        return -1;
    }
    
	// Write fileformat
	write_file_format(vcf_file, fd);
	// Write header entries
	for (list_item_t *i = vcf_file->header_entries->first_p; i != NULL; i = i->next_p) {
		vcf_header_entry_t *entry = i->data_p;
		write_header_entry(entry, fd);
	}
	// Write delimiter
	write_delimiter(vcf_file, fd);
	// Write records
	for (list_item_t *i = vcf_file->records->first_p; i != NULL; i = i->next_p) {
		vcf_record_t *record = i->data_p;
		write_record(record, fd);
	}
	
	return 0;
}

void write_file_format(vcf_file_t *vcf_file, FILE *fd) {
    if(vcf_file == NULL || vcf_file->format == NULL || fd == NULL) {
        return;
    }
    
    fprintf(fd, "##fileformat=%s\n", vcf_file->format);
}

void write_header_entry(vcf_header_entry_t *entry, FILE *fd) {
    if (entry == NULL || fd == NULL) {
        return;
    }

	fprintf(fd, "##");
	
	// Entries with the form ##value
	if (entry->name == NULL) {
		char *value = entry->values->first_p->data_p; // Just first (and only) value
		fprintf(fd, "%s\n", value);
		return;
	}
	
	fprintf(fd, "%s", entry->name);
	// Entries with the form ##name=value
	if (entry->num_keys == 0 && entry->num_values > 0) {
		char *value = entry->values->first_p->data_p; // Just first (and only) value
		fprintf(fd, "=%s\n", value);
	}	
	// Entries with the form ##name=<field_id=value,field_id=value,...>
	else if (entry->num_keys > 0 && entry->num_keys == entry->num_values) {
		fprintf(fd, "=<");
		list_item_t *key = entry->keys->first_p;
		list_item_t *value = entry->values->first_p;
		
		if (strcmp("Description", (char*) key->data_p) != 0)
		{
			fprintf(fd, "%s=%s", (char*) key->data_p, (char*) value->data_p);
		} else 
		{
			fprintf(fd, "%s=\"%s\"", (char*) key->data_p, (char*) value->data_p);
		}
		
		// Get next pair key-value
		key = key->next_p;
		value = value->next_p;
		
		while (key != NULL && value != NULL) {
			if (strcmp("Description", (char*) key->data_p) != 0) {
				fprintf(fd, ",%s=%s", (char*) key->data_p, (char*) value->data_p);
			} else {
				fprintf(fd, ",%s=\"%s\"", (char*) key->data_p, (char*) value->data_p);
			}
			
			// Get next pair key-value
			key = key->next_p;
			value = value->next_p;
		} 
		
		fprintf(fd, ">\n");
	}
}

void write_delimiter(vcf_file_t *vcf_file, FILE *fd) {
    if(vcf_file == NULL || fd == NULL) {
        return;
    }
    
	fprintf(fd, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	
	list_item_t *item = vcf_file->samples_names->first_p;
	while (item != NULL) {
		fprintf(fd, "\t%s", (char*) item->data_p);
		item = item->next_p;
	}
	fprintf(fd, "\n");
}

void write_batch(vcf_batch_t *vcf_batch, FILE *fd) {
    if(vcf_batch == NULL || fd == NULL) {
        return;
    }
    
    for (list_item_t *i = vcf_batch->first_p; i != NULL; i = i->next_p) {
        write_record(i->data_p, fd);
    }
}

void write_record(vcf_record_t* vcf_record, FILE *fd) {
    if(vcf_record == NULL || fd == NULL) {
        return;
    }
    
	fprintf(fd, "%s\t%ld\t%s\t%s\t%s\t", vcf_record->chromosome, vcf_record->position, vcf_record->id, vcf_record->reference, vcf_record->alternate);
	if (vcf_record->quality < 0) {
		fprintf(fd, ".\t");
	} else {
		fprintf(fd, "%.2f\t", vcf_record->quality);
	}
	fprintf(fd, "%s\t%s\t%s", vcf_record->filter, vcf_record->info, vcf_record->format);
	
	list_item_t *item = vcf_record->samples->first_p;
	while (item != NULL) {
		fprintf(fd, "\t%s", (char*) item->data_p);
		item = item->next_p;
	}
	
	fprintf(fd, "\n");
}
