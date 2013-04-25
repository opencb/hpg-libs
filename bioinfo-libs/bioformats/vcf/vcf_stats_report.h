#ifndef VCF_STATS_REPORT_H
#define VCF_STATS_REPORT_H

#include <stdlib.h>
#include <string.h>

#include <bioformats/db/db_utils.h>
#include <commons/file_utils.h>
#include <commons/log.h>
#include <commons/sqlite/sqlite3.h>

#include "vcf_stats.h"


char *get_vcf_file_stats_output_filename(char *vcf_filename, char *out_filename, char *outdir);

void report_vcf_summary_stats(FILE *stats_fd, void *db, file_stats_t *stats);


char *get_variant_stats_output_filename(char *vcf_filename, char *out_filename, char *outdir);

void report_vcf_variant_stats(FILE *stats_fd, void *db, variant_stats_t *stats);

void report_vcf_variant_stats_header(FILE *stats_fd);


void report_sample_stats_output(char *stats_fd, void *db, sample_stats_t *stats);

#endif

