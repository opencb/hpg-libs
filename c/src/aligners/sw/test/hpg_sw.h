#ifndef HPG_SW_H
#define HPG_SW_H

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "macros.h"
#include "sse.h"
#include "smith_waterman.h"

#ifdef SW_WORKFLOW
#include "sw_workflow.h"
#endif

#ifdef SW_OMP
#include "sw_omp.h"
#endif

void run_sw(char *q_filename, char *r_filename,
	    float match, float mismatch,
	    float gap_open, float gap_extend, char *matrix_filename,
	    int batch_size, int num_threads, char *out_filename);

void display_usage(char *msg);
int file_exists(const char *filename);

#endif // HPG_SW_H
