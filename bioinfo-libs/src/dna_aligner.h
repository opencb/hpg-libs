#ifndef DNA_ALIGNER_H
#define DNA_ALIGNER_H
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include "commons/log.h"
#include "commons/file_utils.h"
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

void run_dna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     pair_mng_t *pair_mng, options_t *options);

#endif
