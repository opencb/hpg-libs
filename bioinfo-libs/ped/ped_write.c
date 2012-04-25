#include "ped_write.h"

int ped_write_to_file(ped_file_t *ped_file, FILE *fd) {
    if(ped_file == NULL || fd == NULL) {
        return -1;
    }
    
    family_t **families = cp_hashtable_get_values(ped_file);
    int num_families = get_num_families(ped_file);
    cp_list_iterator *children_iterator;
    
    // Write members of each family
    family_t *family;
    for (int i = 0; i < num_families; i++) {
        family = families[i];
        // Write mother and father
        write_ped_individual(family->father);
        write_ped_individual(family->mother);
        // Write children
        cp_list_iterator *iterator = cp_list_create_iterator(family->children, COLLECTION_LOCK_READ);
        individual *child = NULL;
        while ((child = cp_list_iterator_next(iterator)) != NULL) {
            write_ped_individual(child);
        }
    }
    
    return 0;
}


void write_ped_individual(individual_t* individual, FILE* fd) {
    if (individual == NULL || fd == NULL) {
        return;
    }
    
    fprintf(fd, "%s\t%s\t", individual->family->id, individual->id);
    fprintf(fd, "%s\t%s\t%d\t", (individual->father == NULL) ? 0 : individual->father->id,
                                (individual->mother == NULL) ? 0 : individual->mother->id,
                                individual->sex);
    if (individual->phenotype - ((int) individual->phenotype) > 1e3) {
        fprintf(fd, "%f\t\n", individual->phenotype);
    } else {
        fprintf(fd, "%d\t\n", (int) individual->phenotype);
    }
}



void write_ped_batch(ped_batch_t *ped_batch, FILE *fd) {
    if(ped_batch == NULL || fd == NULL) {
        return;
    }
    for (list_item_t *i = ped_batch->first_p; i != NULL; i = i->next_p) {
        write_ped_record(i->data_p, fd);
    }
}

void write_ped_record(ped_record_t* ped_record, FILE *fd) {
    if(ped_record == NULL || fd == NULL) {
        return;
    }
    
    fprintf(fd, "%s\t%s\t%s\t%s\t%d\t", ped_record->family_id, ped_record->individual_id, ped_record->father_id, ped_record->mother_id, ped_record->sex);
    if (ped_record->phenotype - ((int) ped_record->phenotype) > 1e3) {
        fprintf(fd, "%f\t\n", ped_record->phenotype);
    } else {
        fprintf(fd, "%d\t\n", (int) ped_record->phenotype);
    }
}
