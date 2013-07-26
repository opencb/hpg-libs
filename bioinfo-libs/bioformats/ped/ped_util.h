#ifndef PED_UTIL_H
#define PED_UTIL_H

#include <bioformats/family/family.h>
#include "ped_file.h"
enum Condition get_condition_from_phenotype(char* phenotype, ped_file_t *ped_file);

#endif
