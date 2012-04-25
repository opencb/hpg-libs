#include "family.h"


/*
 * Individual management functions
 */

individual_t *individual_new(char *id, float phenotype, enum Sex sex, individual_t *father, individual_t *mother, family_t *family) {
    individual_t *individual = (individual_t*) malloc (sizeof(invididual_t));
    individual_init(id, phenotype, sex, father, mother, family, individual);
    return individual;
}

void individual_init(char *id, float phenotype, enum Sex sex, individual_t *father, individual_t *mother, family_t *family, individual_t *individual) {
    if (individual == NULL) {
        return;
    }
    
    individual->id = id;
    individual->phenotype = phenotype;
    individual->sex = sex;
    individual->father = father;
    individual->mother = mother;
    individual->family = family;
}

void individual_free(individual_t *individual) {
    if (individual == NULL) {
        return;
    }
    
    free(individual->id);
    free(individual);
}

int individual_compare(individual *a, individual *b) {
    int result = strcasecmp(a->id, b->id);
    if (a->family == NULL && b->family != NULL) {
        result = -1;
    } else if (a->family != NULL && b->family == NULL) {
        result = 1;
    }
    result &= strcasecmp(a->family->id, b->family->id);
    
    return result;
}


/*
 * Family management functions
 */

family_t *family_new(char *id) {
    family_t *family = (family_t*) malloc (sizeof(family);
    family->id = id;
    family->children = cp_list_create_list(COLLECTION_MODE_DEEP,
                                           (cp_compare_fn) individual_compare,
                                           NULL,
                                           individual_free
                                          );
    return family;
}

int family_set_parent(invididual_t *parent, family_t *family) {
    if (parent == NULL) {
        return 1;
    }
    if (family == NULL) {
        return 2;
    }
    if (parent->sex == UNKNOWN) {
        return 3;
    }
    if (family_contains_individual(parent, family)) {
        return 4;
    }
    
    if (parent->sex == MALE) {
        family->father = parent;
    } else if (parent->sex == FEMALE) {
        family->mother = parent;
    }
    return 0;
}

int family_add_child(individual_t *child, family_t *family) {
    if (child == NULL) {
        return 1;
    }
    if (family == NULL) {
        return 2;
    }
    if (family_contains_individual(child, family)) {
        return 3;
    }
    
    return cp_list_append(family->children, child) == NULL;
}

void family_free(family_t *family) {
    if (family == NULL) {
        return;
    }
    
    free(family->id);
    individual_free(family->father);
    individual_free(family->mother);
    cp_list_destroy(family->children);
}

individual *family_contains_individual(invididual_t *individual, family_t *family) {
    if (individual_compare(individual, family->father) == 0) {
        return family->father;
    } 
    if (individual_compare(individual, family->mother) == 0) {
        return family->mother;
    }
    
    cp_list_iterator *iterator = cp_list_create_iterator(family->children, COLLECTION_LOCK_READ);
    individual *child = NULL;
    while ((child = cp_list_iterator_next(iterator)) != NULL) {
        if (individual_compare(individual, child) == 0) {
            return child;
        }
    }
    
    return NULL;
}
