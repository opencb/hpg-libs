#include "ped_file.h"


//====================================================================================
//  ped_file.c
//  ped file management functions
//====================================================================================

//-----------------------------------------------------
// ped_open
//-----------------------------------------------------


ped_file_t *ped_open(char *filename) {
    size_t len;
    char *data = mmap_file(&len, filename);

    ped_file_t *ped_file = (ped_file_t *) malloc(sizeof(ped_file_t));
    ped_file->filename = filename;
    ped_file->data = data;
    ped_file->data_len = len;
    
    ped_file->families = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                                       10,
                                                       cp_hash_istring,         // Hash function
                                                       (cp_compare_fn) strcasecmp,     // Key comparison function
                                                       NULL,                    // Key copy function
                                                       NULL,                    // Key destructor function
                                                       NULL,                    // Value copy function
                                                       (cp_destructor_fn) family_free // Value destructor function
                                                      );
    return ped_file;
}


//-----------------------------------------------------
// ped_close and memory freeing
//-----------------------------------------------------

void ped_close(ped_file_t *ped_file, int free_families) {
    // Free string members
    free(ped_file->filename);
   
    // Free records list if asked to
    if (free_families) {
        cp_hashtable_destroy(ped_file->families);
    }
    
    munmap((void*) ped_file->data, ped_file->data_len);
    free(ped_file);
}

void ped_record_free(ped_record_t* ped_record) {
    free(ped_record->family_id);
    free(ped_record->individual_id);
    free(ped_record->father_id);
    free(ped_record->mother_id);
    free(ped_record);
}

//-----------------------------------------------------
// I/O operations (read and write) in various ways
//-----------------------------------------------------

int ped_read(ped_file_t *ped_file) {
    int ret_code = 0;
    list_t *ped_batches = (list_t*) malloc (sizeof(list_t));
    list_init("batches", 1, 100, ped_batches);
    
#pragma omp parallel sections
{
#pragma omp section
    {
        ret_code = ped_read_batches(ped_batches, 1000, ped_file);
        if (ret_code) {
            LOG_FATAL_F("Error %d while reading the file %s\n", ret_code, ped_file->filename);
        }
        
        list_decr_writers(ped_batches);
    }
    
#pragma omp section
    {
        list_item_t *item = NULL;
        list_item_t *record_item = NULL;
        while ((item = list_remove_item(ped_batches)) != NULL) {
            ped_batch_t *batch = (ped_batch_t*) item->data_p;
            while ((record_item = list_remove_item_async(batch)) != NULL) {
                ret_code &= add_ped_record(record_item->data_p, ped_file);
                ped_record_free(record_item->data_p);
                list_item_free(record_item);
            }
            ped_batch_free(batch);
        }
    }
}
    return ret_code;
}

int ped_read_batches(list_t *batches_list, size_t batch_size, ped_file_t *ped_file) {
    return ped_ragel_read(batches_list, batch_size, ped_file);
}

int ped_write(ped_file_t *ped_file, char *filename) {
    FILE *fd = fopen(filename, "w");
    if (fd < 0) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return 1;
    }
    
    if (ped_write_to_file(ped_file, fd)) {
        fprintf(stderr, "Error writing file: %s\n", filename);
        fclose(fd);
        return 2;
    }
    
    fclose(fd);
    return 0;
}

//-----------------------------------------------------
// load data into the ped_file_t
//-----------------------------------------------------

int add_family(family_t* family, ped_file_t* ped_file) {
    return cp_hashtable_put(ped_file->families, family->id, family) == NULL;
}

int get_num_families(ped_file_t* ped_file) {
    if (ped_file == NULL) {
        return -1;
    }
    return cp_hashtable_count(ped_file->families);
}

int add_ped_record(ped_record_t* record, ped_file_t *ped_file) {
    if (record == NULL) {
        return -1;
    }
    if (ped_file == NULL) {
        return -2;
    }
    
    int result = 0;
    char *aux_buffer;
    
    // Get family or, should it not exist yet, create it
    family_t *family = cp_hashtable_get(ped_file->families, record->family_id);
    if (family == NULL) {
        aux_buffer = (char*) calloc (strlen(record->family_id)+1, sizeof(char));
        strncat(aux_buffer, record->family_id, strlen(record->family_id));
        family = family_new(aux_buffer);
        if (add_family(family, ped_file)) {
            return ALREADY_EXISTING_FAMILY;
        }
    }
    
    LOG_DEBUG_F("father id = %s\tmother id = %s\n", record->father_id, record->mother_id);
    
    // Get parents from family or, should they not exist yet, create them
    individual_t *father = NULL, *mother = NULL;
    if (record->father_id != NULL) {
        if (family->father != NULL) {
            if (strcasecmp(family->father->id, record->father_id) == 0) {
                father = family->father;
            } else {
                return FATHER_APPEARS_MORE_THAN_ONCE;
            }
        } else {
            aux_buffer = (char*) calloc (strlen(record->father_id)+1, sizeof(char));
            strncat(aux_buffer, record->father_id, strlen(record->father_id));
            father = individual_new(aux_buffer, -9, MALE, NULL, NULL, family);
            family_set_parent(father, family);
        }
    }
    
    if (record->mother_id != NULL) {
        if (family->mother != NULL) {
            if (strcasecmp(family->mother->id, record->mother_id) == 0) {
                mother = family->mother;
            } else {
                return MOTHER_APPEARS_MORE_THAN_ONCE;
            }
        } else {
            aux_buffer = (char*) calloc (strlen(record->mother_id)+1, sizeof(char));
            strncat(aux_buffer, record->mother_id, strlen(record->mother_id));
            mother = individual_new(aux_buffer, -9, FEMALE, NULL, NULL, family);
            family_set_parent(mother, family);
        }
    }
    
    // Create individual with the information extracted from the PED record
    aux_buffer = (char*) calloc (strlen(record->individual_id)+1, sizeof(char));
    strncat(aux_buffer, record->individual_id, strlen(record->individual_id));
    individual_t *individual = individual_new(aux_buffer, record->phenotype, record->sex, father, mother, family);
    if (father != NULL || mother != NULL) {
        LOG_DEBUG_F("** add family %s child (id %s)\n", family->id, individual->id);
        family_add_child(individual, family);
    } else {
        LOG_DEBUG_F("** set family %s parent of sex %d (id %s)\n", family->id, individual->sex, individual->id);
        result = family_set_parent(individual, family);
        if (result == 1) {
            result = FATHER_APPEARS_MORE_THAN_ONCE;
        } else if (result == 2) {
            result = MOTHER_APPEARS_MORE_THAN_ONCE;
        }
    }

    return result;
}

