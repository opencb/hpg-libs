#ifndef PED_FILE_STRUCTURE_H
#define PED_FILE_STRUCTURE_H

#include <sys/types.h>


#include <commons/file_utils.h>

#include <containers/cprops/hashtable.h>
#include <containers/khash.h>

#include <bioformats/family/family.h>

KHASH_MAP_INIT_STR(str, int);

/**
 * Entry in the PED document body.
 */
typedef struct ped_record {
    char *family_id;
    char *individual_id;
    char *father_id;
    char *mother_id;
    enum Sex sex;
    int phenotype;
} ped_record_t;

/**
 * PED file structure. The physical file is defined by its file descriptor, its 
 * filename and the mode it has been open.
 *
 * It contains a collection of families with a unique ID, which are specified via 
 * a group of records which contain their members.
 */
typedef struct ped_file {
    char* filename;
    char* mode;
    char *data;
    size_t data_len;
    
    cp_hashtable *families;
    
    //TODO: Will have extra value for the number of the field to compare
    int unaffected_id;
    int affected_id;
    khash_t(str) *phenotypes;
} ped_file_t;

#endif
