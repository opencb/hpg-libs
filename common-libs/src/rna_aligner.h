#ifndef RNA_ALIGNER_H
#define RNA_ALIGNER_H
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include "commons/log.h"
#include "error.h"
#include "timing.h"
#include "buffers.h"
#include "bioformats/fastq/fastq_batch_reader.h"
#include "bwt_server.h"
#include "batch_writer.h"
#include "region_seeker.h"
#include "cal_seeker.h"
#include "sw_server.h"
#include "pair_server.h"
#include "batch_aligner.h"
#include "options.h"
#include "statistics.h"
// rna server
#include "rna_splice.h"
#include "rna_server.h"


void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		     pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     options_t *options);

#endif
