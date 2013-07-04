#ifndef PED_UTIL_H
#define PED_UTIL_H

#include <bioformats/family/family.h>
/*TODO Like in ped_file.h:101, it'll take the condition from a selected field.
 * 		enum Condition get_condition(individual_t *individual, ped_file_t *ped_file);
 * Also, will change from a getter to setter like:
 * 		int set_condition(individual_t *individual, ped_file_t *ped_file);
 * 
 * Notice that it's not an "util". Will move to ped_file.h/c. Meanwhile, #include "ped_file.h"
 * */
#include "ped_file.h"
enum Condition get_condition_from_phenotype(int phenotype, ped_file_t *ped_file);

#endif
