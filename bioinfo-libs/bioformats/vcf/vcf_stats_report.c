#include "vcf_stats_report.h"


/* ***********************************************
 *                   Global report               *
 * ***********************************************/

char *get_vcf_file_stats_output_filename(char *vcf_filename, char *out_filename, char *outdir) {
    char *summary_filename;
    
    if (out_filename == NULL || strlen(out_filename) == 0) {
        char suffix_filename[strlen(vcf_filename)];
        get_filename_from_path(vcf_filename, suffix_filename);
        
        summary_filename = malloc ((strlen(outdir) + strlen(suffix_filename) + strlen(".summary-stats") + 2) * sizeof(char));
        sprintf(summary_filename, "%s/%s.stats-summary", outdir, suffix_filename);
    } else {
        summary_filename = (char*) calloc ((strlen(outdir) + strlen(out_filename) + strlen(".summary-stats") + 2), sizeof(char));
        sprintf(summary_filename, "%s/%s.stats-summary", outdir, out_filename);
    }
    
    return summary_filename;
}

static void report_summary(FILE *summary_fd, file_stats_t *file_stats) {
    // Write whole file stats (data only got when launching variant stats)
    fprintf(summary_fd, 
            "Number of variants = %d\nNumber of samples = %d\nNumber of biallelic variants = %d\nNumber of multiallelic variants = %d\n\n",
            file_stats->variants_count, file_stats->samples_count, file_stats->biallelics_count, file_stats->multiallelics_count);
    
    fprintf(summary_fd, 
            "Number of SNP = %d\nNumber of indels = %d\n\n",
            file_stats->snps_count, file_stats->indels_count);

    fprintf(summary_fd, 
            "Number of transitions = %d\nNumber of transversions = %d\nTi/TV ratio = %.3f\n\n",
            file_stats->transitions_count, file_stats->transversions_count,
            (float) file_stats->transitions_count / file_stats->transversions_count);
    
    fprintf(summary_fd,
            "Percentage of PASS = %.3f%%\nAverage quality = %.3f\n",
            ((float) file_stats->pass_count / file_stats->variants_count) * 100.0,
            file_stats->accum_quality / file_stats->variants_count);

}


static void report_summary_sqlite3(sqlite3 *db, file_stats_t *output) {
    sqlite3_stmt* stmt;
    char* error_message;
    char aux[64];

    prepare_statement_global_stats(db, &stmt);

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &error_message);

    sprintf(aux, "%d", output->variants_count);
    insert_statement_global_stats("NUM_VARIANTS", "Number of variants", aux, stmt, db);

    sprintf(aux, "%d", output->samples_count);
    insert_statement_global_stats("NUM_SAMPLES", "Number of samples", aux, stmt, db);

    sprintf(aux, "%d", output->biallelics_count);
    insert_statement_global_stats("NUM_BIALLELIC", "Number of biallelic variants", aux, stmt, db);

    sprintf(aux, "%d", output->multiallelics_count);
    insert_statement_global_stats("NUM_MULTIALLELIC", "Number of multiallelic variants", aux, stmt, db);

    sprintf(aux, "%d", output->snps_count);
    insert_statement_global_stats("NUM_SNPS", "Number of SNPs", aux, stmt, db);

    sprintf(aux, "%d", output->indels_count);
    insert_statement_global_stats("NUM_INDELS", "Number of indels", aux, stmt, db);

    sprintf(aux, "%d", output->transitions_count);
    insert_statement_global_stats("NUM_TRANSITIONS", "Number of transitions", aux, stmt, db);

    sprintf(aux, "%d", output->transversions_count);
    insert_statement_global_stats("NUM_TRANSVERSIONS", "Number of transversions", aux, stmt, db);

    sprintf(aux, "%.3f", (float) output->transitions_count / output->transversions_count);
    insert_statement_global_stats("TITV_RATIO", "Ti/TV ratio", aux, stmt, db);

    sprintf(aux, "%.3f", ((float) output->pass_count / output->variants_count) * 100.0);
    insert_statement_global_stats("PERCENT_PASS", "Percentage of PASS", aux, stmt, db);

    sprintf(aux, "%.3f", output->accum_quality / output->variants_count);
    insert_statement_global_stats("AVG_QUALITY", "Average quality", aux, stmt, db);

    sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &error_message);

    finalize_statement_global_stats(stmt);
}


void report_vcf_summary_stats(FILE *stats_fd, void *db, file_stats_t *stats) {
    // Write to plain text file
    report_summary(stats_fd, stats);

    // Write to database (optional)
    if (db) {
        report_summary_sqlite3((sqlite3 *) db, stats);
    }
}


/* ***********************************************
 *                 Variants report               *
 * ***********************************************/

char *get_variant_stats_output_filename(char *vcf_filename, char *out_filename, char *outdir) {
    char *stats_filename;
    
    if (out_filename == NULL || strlen(out_filename) == 0) {
        char suffix_filename[strlen(vcf_filename) + 1];
        get_filename_from_path(vcf_filename, suffix_filename);
        
        stats_filename = (char*) calloc ((strlen(outdir) + strlen(suffix_filename) + strlen(".stats-variants") + 2), sizeof(char));
        sprintf(stats_filename, "%s/%s.stats-variants", outdir, suffix_filename);
    } else {
        stats_filename = (char*) calloc ((strlen(outdir) + strlen(out_filename) + strlen(".stats-variants") + 2), sizeof(char));
        sprintf(stats_filename, "%s/%s.stats-variants", outdir, out_filename);
    }
    
    return stats_filename;
}

static int report_variant_alleles_stats(variant_stats_t *var_stats, FILE *stats_fd) {
    int written = 0;
    
    // Is indel?
    written += (var_stats->is_indel) ? fprintf(stats_fd, "Y\t") : fprintf(stats_fd, "N\t");
    
    // Reference allele
    written += fprintf(stats_fd, "%s\t%d\t%.4f\t",
                       var_stats->ref_allele,
                       var_stats->alleles_count[0],
                       var_stats->alleles_freq[0]);

    // Alternate alleles
    for (int i = 1; i < var_stats->num_alleles; i++) {
        written += fprintf(stats_fd, "%s\t%d\t%.4f\t",
                           var_stats->alternates[i-1],
                           var_stats->alleles_count[i],
                           var_stats->alleles_freq[i]);
    }
    
    return written;
} 
    
static int report_variant_genotypes_stats(variant_stats_t *var_stats, FILE *stats_fd) {
    int written = 0;
    int gt_count = 0;
    float gt_freq = 0;
    
    for (int i = 0; i < var_stats->num_alleles; i++) {
        for (int j = i; j < var_stats->num_alleles; j++) {
            int idx1 = i * var_stats->num_alleles + j;
            if (i == j) {
                gt_count = var_stats->genotypes_count[idx1];
                gt_freq = var_stats->genotypes_freq[idx1];
            } else {
                int idx2 = j * var_stats->num_alleles + i;
                gt_count = var_stats->genotypes_count[idx1] + var_stats->genotypes_count[idx2];
                gt_freq = var_stats->genotypes_freq[idx1] + var_stats->genotypes_freq[idx2];
            }

            written += fprintf(stats_fd, "%s|%s\t%d\t%.4f\t",
                               i == 0 ? var_stats->ref_allele : var_stats->alternates[i-1],
                               j == 0 ? var_stats->ref_allele : var_stats->alternates[j-1],
                               gt_count, gt_freq);
        }
    }
    
    return written;
}

static inline int report_variant_missing_data(variant_stats_t *var_stats, FILE *stats_fd) {
    return fprintf(stats_fd, "%d\t%d\t",
                   var_stats->missing_alleles,
                   var_stats->missing_genotypes);
}

static inline int report_variant_inheritance_data(variant_stats_t *var_stats, FILE *stats_fd) {
    return fprintf(stats_fd, "%d\t%.2f | %.2f\t%.2f | %.2f\n",
                   var_stats->missing_alleles,
                   var_stats->cases_percent_dominant,
                   var_stats->controls_percent_dominant,
                   var_stats->cases_percent_recessive,
                   var_stats->controls_percent_recessive);
}

void report_vcf_variant_stats(FILE *stats_fd, void *db, variant_stats_t *stats) {
    // Write to plain text file
    fprintf(stats_fd, "%s\t%ld\t", stats->chromosome, stats->position);

    report_variant_alleles_stats(stats, stats_fd);
    report_variant_genotypes_stats(stats, stats_fd);
    report_variant_missing_data(stats, stats_fd);
    report_variant_inheritance_data(stats, stats_fd);

    // Write to database (optional)
    if (db) {
        //report_summary_sqlite3((sqlite3 *) db, stats);
    }
}

inline void report_vcf_variant_stats_header(FILE *stats_fd) {
    fprintf(stats_fd, 
            "#CHROM\tPOS\tINDEL?\tList of [ALLELE  COUNT  FREQ]\t\t\tList of [GT  COUNT  FREQ]\t\t\t\t\t\tMISS_AL\tMISS_GT\tMEND_ER\t%% AFF | UNAFF dominant\t%% AFF | UNAFF recessive\n");
}


/* ***********************************************
 *                 Samples report                *
 * ***********************************************/

char *get_sample_stats_output_filename(char *vcf_filename, char *out_filename, char *outdir) {
    char *stats_filename;
    
    if (out_filename == NULL || strlen(out_filename) == 0) {
        char suffix_filename[strlen(vcf_filename) + 1];
        get_filename_from_path(vcf_filename, suffix_filename);
        
        stats_filename = (char*) calloc ((strlen(outdir) + strlen(suffix_filename) + strlen(".stats-samples") + 2), sizeof(char));
        sprintf(stats_filename, "%s/%s.stats-samples", outdir, suffix_filename);
    } else {
        stats_filename = (char*) calloc ((strlen(outdir) + strlen(out_filename) + strlen(".stats-variants") + 2), sizeof(char));
        sprintf(stats_filename, "%s/%s.stats-samples", outdir, out_filename);
    }
    
    return stats_filename;
}

void report_sample_stats(FILE *stats_fd, void *db, size_t num_samples, sample_stats_t **stats) {
    sample_stats_t *sam_stats;
    for (int i = 0; i < num_samples; i++) {
        sam_stats = stats[i];
        fprintf(stats_fd, "%s\t\t%zu\t\t%zu\n", sam_stats->name, sam_stats->missing_genotypes, sam_stats->mendelian_errors);
        sample_stats_free(sam_stats);
    }
}

inline void report_sample_variant_stats_header(FILE *stats_fd) {
    fprintf(stats_fd, "#SAMPLE\t\tMISS GT\t\tMENDEL ERR\n");
}
