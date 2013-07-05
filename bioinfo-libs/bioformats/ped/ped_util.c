#include "ped_util.h"

enum Condition get_condition_from_phenotype(int phenotype, ped_file_t *ped_file) {
	if(ped_file->unaffected_id == -1)
        ped_file->unaffected_id = kh_get(str, ped_file->phenotypes, "1");
    if(ped_file->affected_id == -1)
        ped_file->affected_id = kh_get(str, ped_file->phenotypes, "2");
    
    
    
    if (phenotype == ped_file->unaffected_id) {//printf("ped_util.c : Phenotype %d es unaffected \n", phenotype);
        return UNAFFECTED;
    } else if (phenotype == ped_file->affected_id) {//printf("ped_util.c : Phenotype %d es affected\n", phenotype);
        return AFFECTED;
    } else {//printf("ped_util.c : Phenotype %d es desconocido\n", phenotype);
        return UNKNOWN_CONDITION;
    }
}
