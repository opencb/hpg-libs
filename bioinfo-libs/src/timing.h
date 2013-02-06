#ifndef TIMING_H
#define TIMING_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>


typedef struct timing {
  int num_sections;

  char** section_labels_p;
  int* num_threads_p;

  double* section_times_p;
  double** thread_times_p;
  struct timeval** start_times_p;
  struct timeval** stop_times_p;
} timing_t;


timing_t* timing_new(char** section_labels, int* num_threads_p, int num_sections);
void timing_free(timing_t* timing_p);

void timing_start(int section_id, int thread_id, timing_t* timing_p);
void timing_stop(int section_id, int thread_id, timing_t* timing_p);
void timing_display(timing_t* timing_p);


extern char time_on;
extern timing_t *timing_p;

// for debugging
extern double bwt_time[100], seeding_time[100], cal_time[100], sw_time[100];
extern int thr_bwt_items[100], thr_seeding_items[100], thr_cal_items[100], thr_sw_items[100];
extern int thr_batches[100];

#endif // end of if TIMING






