#ifndef FAMILY_H
#define FAMILY_H

#include <stdlib.h>

#include <cprops/linked_list.h>

#include <string_utils.h>


enum Sex { UNKNOWN, MALE, FEMALE };


/**
 * Entry in the PED document body, representing an individual and member of a family.
 */
typedef struct individual {
    char *id;   /**< Unique ID of the individual **/
    float phenotype;    /**< Numerical descriptor for the affection of the individual */
    enum Sex sex;   /**< Sex of the individual */
    struct individual *father;  /**< Father of the individual (NULL in case he's parent in a family) */
    struct individual *mother;  /**< Mother of the individual (NULL in case he's parent in a family) */
    struct family *family;  /**< Family of the individual **/
} individual_t;

/**
 * Family described in a PED document
 */
typedef struct family {
    char *id;   /**< Unique ID of the family **/
    individual_t *father;  /**< Man in the root of the genealogical tree */
    individual_t *mother;  /**< Woman in the root of the genealogical tree */
    cp_list *children;  /**< Children of the main roots in the genealogical tree */
} family_t;


/*
 * Individual management functions
 */

/**
 * Creates a new individual with the characteristics provided as arguments.
 * 
 * @param id unique ID of the individual
 * @param phenotype numerical descriptor for the affection of the individual
 * @param sex sex of the individual
 * @param father father of the individual, if apply
 * @param mother mother of the individual, if apply
 * @param family family of the individual
 * @return The newly created individual
 */
individual_t *individual_new(char *id, float phenotype, enum Sex sex, individual_t *father, individual_t *mother, family_t *family);

/**
 * Fills member of an already existing individual with the characteristics provided as arguments.
 * 
 * @param id unique ID of the individual
 * @param phenotype numerical descriptor for the affection of the individual
 * @param sex sex of the individual
 * @param father father of the individual, if apply
 * @param mother mother of the individual, if apply
 * @param family family of the individual
 * @param individual individual whose members are being filled
 */
void individual_init(char *id, float phenotype, enum Sex sex, individual_t *father, individual_t *mother, family_t *family, individual_t *individual);

/**
 * Free memory associated to an individual.
 * 
 * @param individual individual to be freed
 */
void individual_free(individual_t *individual);


/*
 * Family management functions
 */

/**
 * Creates a new family with the characteristics provided as arguments.
 * 
 * @param id unique ID of the family
 * @return The newly created family
 * 
 */
family_t *family_new(char *id);

/**
 * Sets an individual as parent of a family (its position as father or mother is automatically 
 * assigned based on its sex). An individual can only be set as parent if his sex is known and 
 * he is not part of the family yet.
 * 
 * @param parent the individual to set as parent
 * @param family the family where the individual could belong to
 * @return 0 if the individual was succcessfully set, 1-2 if one of the arguments is NULL, 
 * 3 if the parent sex in unknown, or 4 if the family already contains the individual
 */
int family_set_parent(individual_t *parent, family_t *family);

/**
 * Sets an individual as child of a family. An individual can only be set as child if he is 
 * not part of the family yet.
 * 
 * @param child the individual to set as child
 * @param family the family where the individual could belong to
 * @return 0 if the individual was succcessfully set, 1-2 if one of the arguments is NULL, 
 * 3 if the family already contains the individual
 */
int family_add_child(individual_t *child, family_t *family);

/**
 * Free memory associated to a family and its individuals.
 * 
 * @param family the family to be freed
 */
void family_free(family_t *family);

/**
 * Checks whether an individual is part of a family.
 * 
 * @param individual the individual to check
 * @param family the family that could contain the individual
 */
individual_t *family_contains_individual(individual_t *individual, family_t *family);

#endif
