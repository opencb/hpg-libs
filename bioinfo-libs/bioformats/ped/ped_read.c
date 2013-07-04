#include "ped_read.h"

//====================================================================================
// ped reading functions
//====================================================================================

ped_record_t* create_ped_record() {
    ped_record_t *record = (ped_record_t*) malloc (sizeof(ped_record_t));
    record->family_id = NULL;
    record->individual_id = NULL;
    record->father_id = NULL;
    record->mother_id = NULL;
    record->sex = UNKNOWN_SEX;
    record->phenotype = -9;
    return record;
}

void set_ped_record_family_id(char* family_id, ped_record_t* ped_record) {
    ped_record->family_id = family_id;
    LOG_DEBUG_F("set family_id: %s\n", ped_record->family_id);
}

void set_ped_record_individual_id(char* individual_id, ped_record_t* ped_record) {
    ped_record->individual_id = individual_id;
    LOG_DEBUG_F("set individual_id: %s\n", ped_record->individual_id);
}

void set_ped_record_father_id(char* father_id, ped_record_t* ped_record) {
    ped_record->father_id = father_id;
    LOG_DEBUG_F("set father_id: %s\n", ped_record->father_id);
}

void set_ped_record_mother_id(char* mother_id, ped_record_t* ped_record) {
    ped_record->mother_id = mother_id;
    LOG_DEBUG_F("set mother_id: %s\n", ped_record->mother_id);
}

void set_ped_record_sex(enum Sex sex, ped_record_t* ped_record) {
    ped_record->sex = sex;
    LOG_DEBUG_F("set sex: %ld\n", ped_record->sex);
}

void set_ped_record_phenotype(char* phenotype, ped_record_t* ped_record, ped_file_t *ped_file) {
    
    int ret = 888;
    int k = kh_get(str, ped_file->phenotypes, phenotype);
    if(k == kh_end(ped_file->phenotypes))	//If k == kh_end, is missing
    {
		//printf("Added phenotype. k = %d\n", k);
		k = kh_put(str, ped_file->phenotypes, strdup(phenotype), &ret);
		kh_value(ped_file->phenotypes, k) = kh_size(ped_file->phenotypes)-1;//ped_file->num_phenotypes;
    }
    ped_record->phenotype = kh_value(ped_file->phenotypes, k);
    //printf("ped_read.c:48: El string %s tiene el khiter %d. Value %d.    RET = %d\n", phenotype, k, ped_record->phenotype, ret);
	
    LOG_DEBUG_F("set phenotype: %f\n", ped_record->phenotype);
}
