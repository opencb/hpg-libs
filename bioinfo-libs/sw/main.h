#ifndef MAIN_H
#define MAIN_H

#include "macros.h"
#include "emboss.h"
#include "sse.h"
#include "smith_waterman.h"
//#include "avx.h"

void run_emboss(char *q_filename, char *r_filename, 
		float gap_open, float gap_extend,
		float match, float mismatch,
		int batch_size, int num_threads, char *out_filename);

void run_sse(char *q_filename, char *r_filename, 
	     float gap_open, float gap_extend, char *matrix_filename,
	     int batch_size, int num_threads, char *out_filename);

#endif // MAIN_H
