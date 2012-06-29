#ifndef VCF_STATS_H
#define VCF_STATS_H

#include <stdlib.h>
#include <string.h>

#include <containers/array_list.h>
#include <containers/list.h>
#include <commons/log.h>

#include "vcf_file_structure.h"
#include "vcf_util.h"


typedef struct {
    int variants_count;
    int samples_count;
    
    int snps_count;
    int indels_count;
    
    int transitions_count;
    int transversions_count;
    
    int biallelics_count;
    int multiallelics_count;
    
    int pass_count;
    float accum_quality;
    float mean_quality;
} file_stats_t;


typedef struct {
    char *chromosome;
    unsigned long position;
    
    char *ref_allele;
    char **alternates;
    
    int num_alleles;
    int *alleles_count;
    int *genotypes_count;
    float *alleles_freq;
    float *genotypes_freq;
    
    int missing_alleles;
    int missing_genotypes;
} variant_stats_t;


/* ********************************
 * Initialization and destruction *
 * ********************************/

/**
 * Initialize a global_stats_t structure mandatory fields.
 */
file_stats_t *new_file_stats();

/**
 * Free memory associated to a global_stats_t structure.
 */
void free_file_stats(file_stats_t *file_stats);


/**
 * Initialize a variant_stats_t structure mandatory fields.
 */
variant_stats_t *new_variant_stats(char *chromosome, unsigned long position, char *ref_allele);

/**
 * Free memory associated to a variant_stats_t structure.
 */
void free_variant_stats(variant_stats_t *variant_stats);



/* ******************************
 *           Execution          *
 * ******************************/

int get_variants_stats(vcf_record_t **variants, int num_variants, list_t *output_list, file_stats_t *file_stats);

/**
 * Given a file_stats_t with some values, perform their sum with the arguments.
 * 
 */
void update_file_stats(int variants_count, int samples_count, int snps_count, int transitions_count, int transversions_count, 
                       int indels_count, int biallelics_count, int multiallelics_count, int pass_count, float accum_quality, 
                       file_stats_t *stats);


#endif
