#include "vcf_write.h"

int write_vcf_file(vcf_file_t *file, FILE *fd) {
    assert(file);
    assert(fd);
    
    // Write header: file format, header entries and delimiter
    if (write_vcf_header(file, fd) > 0) {
        return 1;
    }
    
    // Write records: grouped in batches
    vcf_batch_t *batch;
    while ((batch = fetch_vcf_batch(file)) != NULL) {
        if (write_vcf_batch(batch, fd) > 0) {
            return 1;
        }
    }
    
    return 0;
}

int write_vcf_header(vcf_file_t* file, FILE* fd) {
    assert(file);
    assert(fd);
    
    // Write fileformat
    if (write_vcf_fileformat(file, fd) > 0) {
        return 1;
    }
    
    // Write header entries
    for (int i = 0; i < file->header_entries->size; i++) {
        if (write_vcf_header_entry((vcf_header_entry_t*) array_list_get(i, file->header_entries), fd) > 0) {
            return 1;
        }
    }
    
    // Write delimiter
    if (write_vcf_delimiter(file, fd) > 0) {
        return 1;
    }
    
    return 0;
}


int write_vcf_fileformat(vcf_file_t *file, FILE *fd) {
    assert(file);
    assert(fd);
    
    if (fprintf(fd, "##fileformat=%s\n", file->format) < 0) {
        return 1;
    }
    
    return 0;
}

int write_vcf_header_entry(vcf_header_entry_t *entry, FILE *fd) {
    assert(entry);
    assert(fd);
    
    if (fprintf(fd, "##") < 0) {
        return 1;
    }

    // Entries with the form ##value
    if (entry->name == NULL) {
        char *value = entry->values->first_p->data_p; // Just first (and only) value
        if (fprintf(fd, "%s\n", value) < 0) {
            return 1;
        }
        return 0;
    }
    
    // Name for entries with the form ##name=...
    if (fprintf(fd, "%s", entry->name) < 0) {
        return 1;
    }
    
    // Entries with the form ##name=value
    if (entry->num_keys == 0 && entry->num_values > 0) {
        char *value = entry->values->first_p->data_p; // Just first (and only) value
        if (fprintf(fd, "=%s\n", value) < 0) {
            return 1;
        }
    }
    
    // Entries with the form ##name=<field_id=value,field_id=value,...>
    else if (entry->num_keys > 0 && entry->num_keys == entry->num_values) {
        if (fprintf(fd, "=<") < 0) {
            return 1;
        }
        
        list_item_t *key = entry->keys->first_p;
        list_item_t *value = entry->values->first_p;
        
        if (strcmp("Description", (char*) key->data_p) != 0) {
            if (fprintf(fd, "%s=%s", (char*) key->data_p, (char*) value->data_p) < 0) {
                return 1;
            }
        } else {
            if (fprintf(fd, "%s=\"%s\"", (char*) key->data_p, (char*) value->data_p) < 0) {
                return 1;
            }
        }
        
        // Get next pair key-value
        key = key->next_p;
        value = value->next_p;
        
        while (key != NULL && value != NULL) {
            if (strcmp("Description", (char*) key->data_p) != 0) {
                if (fprintf(fd, ",%s=%s", (char*) key->data_p, (char*) value->data_p) < 0) {
                    return 1;
                }
            } else {
                if (fprintf(fd, ",%s=\"%s\"", (char*) key->data_p, (char*) value->data_p) < 0) {
                    return 1;
                }
            }
            
            // Get next pair key-value
            key = key->next_p;
            value = value->next_p;
        } 
        
        if (fprintf(fd, ">\n") < 0) {
            return 1;
        }
    }
    
    return 0;
}

int write_vcf_delimiter(vcf_file_t *file, FILE *fd) {
    assert(file);
    assert(fd);
    
    if (fprintf(fd, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT") < 0) {
        return 1;
    }
    
    for (int i = 0; i < file->samples_names->size; i++) {
        if (fprintf(fd, "\t%s", (char*) file->samples_names->items[i]) < 0) {
            return 1;
        }
    }
       
    if (fprintf(fd, "\n") < 0) {
        return 1;
    }
    
    return 0;
}

int write_vcf_batch(vcf_batch_t *batch, FILE *fd) {
    assert(batch);
    assert(fd);
    
    vcf_record_t *record;
    for (int i = 0; i < batch->records->size; i++) {
//         record = array_list_get(i, batch->records);
        record = batch->records->items[i];
        if (write_vcf_record(record, fd) > 0) {
            return 1;
        }
    }
    
    return 0;
}

int write_vcf_record(vcf_record_t* record, FILE *fd) {
    assert(record);
    assert(fd);
    
    if (fprintf(fd, "%s\t%ld\t%s\t%s\t%s\t", record->chromosome, record->position, record->id, record->reference, record->alternate) < 0) {
        return 1;
    }
    
    if (record->quality < 0) {
        if (fprintf(fd, ".\t") < 0) {
            return 1;
        }
    } else {
        if (fprintf(fd, "%.2f\t", record->quality) < 0) {
            return 1;
        }
    }
    if (fprintf(fd, "%s\t%s\t%s", record->filter, record->info, record->format) < 0) {
        return 1;
    }

    for (int i = 0; i < record->samples->size; i++) {
        if (fprintf(fd, "\t%s", (char*) array_list_get(i, record->samples)) < 0) {
            return 1;
        }
    }

    if (fprintf(fd, "\n") < 0) {
        return 1;
    }
    
    return 0;
}
