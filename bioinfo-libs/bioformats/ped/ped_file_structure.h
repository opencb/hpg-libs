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
    char* phenotype;
    char* custom_field;     /**< Value in the variable column. */
    
    int var_index;        /**< Variable index in the variables' khash */
    
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

    char* unaffected;
    char* affected;
    khash_t(str) *variables;   /**<  Differents values in the variable field */
    int num_variables;         /**<  Number of differents phenotypes */
    int accept_new_values;      /**<  Boolean. New values will be accepted in the khash while reading the ped file  */

    char* variable_field;         /**<  Name of the variable field */
    int num_field;              /**<  Number of the column from the variable field */

} ped_file_t;

#endif
