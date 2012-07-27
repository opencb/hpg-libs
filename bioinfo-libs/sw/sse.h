#ifndef SSE_H
#define SSE_H

#include <immintrin.h>

#include "macros.h"

//------------------------------------------------------------------------
// SSE functions
//------------------------------------------------------------------------
/*
void sw_sse(char **q, char **r, float gap_open, float gap_extend, 
	    float *score, char **m, char **n, int *start1, int *start2);

void sse_matrix(int num_seqs, 
		char **q, int *q_len, int max_q_len,
		char **r, int *r_len, int max_r_len,
		float match, float mismatch,
		float gap_open, float gap_extend,
		float *H, int *C, float *max_score);

void sse_traceback(int num_seqs, 
		   char **q, int *q_len, int max_q_len,
		   char **r, int *r_len, int max_r_len,
		   float gap_open, float gap_extend,
		   float *H, int *C, float *max_score,
		   char **q_alig, int *q_start,
		   char **r_alig, int *r_start, 
		   int *len_alig);
*/
//------------------------------------------------------------------------

void sw_sse(char **q, char **r, float gap_open, float gap_extend, 
	    float *score, char **m, char **n, int *start1, int *start2);

void sse_matrix(int num_seqs, 
		char **q, int *q_len, int max_q_len,
		char **r, int *r_len, int max_r_len,
		float profile[128][128], float gap_open, float gap_extend,
		float *H, float *D, float *E, float *F, float *max_score);

void sw_sse1(char **q, char **r, float gap_open, float gap_extend, 
	     float *score, char **m, char **n, int *start1, int *start2);

void sse1_matrix(int num_seqs, 
		 char **q, int *q_len, int max_q_len,
		 char **r, int *r_len, int max_r_len,
		 float profile[128][128], float gap_open, float gap_extend,
		 float *H, float *F, int *C, float *max_score);

/*
void sse_traceback2(int num_seqs, 
		    char **q, int *q_len, int max_q_len,
		    char **r, int *r_len, int max_r_len,
		    float gap_open, float gap_extend,
		    float *H, float *D, float *E, float *F, 
		    float *max_score,
		    char **q_alig, int *q_start,
		    char **r_alig, int *r_start, 
		    int *len_alig);

//------------------------------------------------------------------------

void sse_find_position(int index, char *q, int q_len, char *r, int r_len,
		       float *H, int cols, int row, float score, 
		       int *q_pos, int *r_pos);
*/
//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // SSE_H
