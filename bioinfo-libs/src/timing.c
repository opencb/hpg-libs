#include "timing.h"

//-----------------------------------------------------

double bwt_time[100], seeding_time[100], cal_time[100], sw_time[100];
int thr_bwt_items[100], thr_seeding_items[100], thr_cal_items[100], thr_sw_items[100];
int thr_batches[100];

//-----------------------------------------------------

timing_t* timing_new(char** section_labels_p, int* num_threads_p, int num_sections) {

  // for debugging
  for (int i = 0; i < 100; i++) {
    bwt_time[i] = 0;
    seeding_time[i] = 0;
    cal_time[i] = 0;
    sw_time[i] = 0;
    thr_batches[i] = 0;
    thr_bwt_items[i] = 0;
    thr_seeding_items[i] = 0;
    thr_cal_items[i] = 0;
    thr_sw_items[i] = 0;
  }

  timing_t* t_p = (timing_t*) calloc(1, sizeof(timing_t));

  t_p->num_sections = num_sections;

  t_p->section_labels_p = (char**) calloc(num_sections, sizeof(char*));
  t_p->num_threads_p = (int*) calloc(num_sections, sizeof(int));
  t_p->section_times_p = (double*) calloc(num_sections, sizeof(double));

  t_p->thread_times_p = (double**) calloc(num_sections, sizeof(double*));
  t_p->start_times_p = (struct timeval**) calloc(num_sections, sizeof(struct timeval*));
  t_p->stop_times_p = (struct timeval**) calloc(num_sections, sizeof(struct timeval*));

  int i, j;
  for(i=0 ; i<num_sections ; i++) {
    t_p->section_labels_p[i] = strdup(section_labels_p[i]);
    t_p->num_threads_p[i] = num_threads_p[i];

    t_p->thread_times_p[i] = (double*) calloc(num_threads_p[i], sizeof(double));
    t_p->start_times_p[i] = (struct timeval*) calloc(num_threads_p[i], sizeof(struct timeval));
    t_p->stop_times_p[i] = (struct timeval*) calloc(num_threads_p[i], sizeof(struct timeval));
  }
  
  return t_p;
} 

//-----------------------------------------------------

void timing_free(timing_t* t_p) {
  if (t_p == NULL) return;

  int i, j;
  for(i=0 ; i<t_p->num_sections ; i++) {
    if (t_p->section_labels_p[i] != NULL) { free(t_p->section_labels_p[i]); } 

    if (t_p->thread_times_p[i] != NULL) { free(t_p->thread_times_p[i]); }
    if (t_p->start_times_p[i] != NULL) {free(t_p->start_times_p[i]); } 
    if (t_p->stop_times_p[i] != NULL) { free(t_p->stop_times_p[i]); }
  }
  
  if (t_p->section_labels_p != NULL) { free(t_p->section_labels_p); }

  if (t_p->thread_times_p != NULL) { free(t_p->thread_times_p); }
  if (t_p->start_times_p != NULL) { free(t_p->start_times_p); }
  if (t_p->stop_times_p != NULL) { free(t_p->stop_times_p); }

  if (t_p->num_threads_p != NULL) { free(t_p->num_threads_p); }
  if (t_p->section_times_p != NULL) { free(t_p->section_times_p); }
  
  free(t_p);
} 

//-----------------------------------------------------

void timing_start(int section_id, int thread_id, timing_t* t_p) {
  if (t_p == NULL) return;
  gettimeofday(&t_p->start_times_p[section_id][thread_id], NULL);
} 

//-----------------------------------------------------

void timing_stop(int section_id, int thread_id, timing_t* t_p) {
  if (t_p == NULL) return;
  
  gettimeofday(&t_p->stop_times_p[section_id][thread_id], NULL);
  t_p->thread_times_p[section_id][thread_id] += ((t_p->stop_times_p[section_id][thread_id].tv_sec - t_p->start_times_p[section_id][thread_id].tv_sec) * 1e6 + (t_p->stop_times_p[section_id][thread_id].tv_usec - t_p->start_times_p[section_id][thread_id].tv_usec));
} 

//-----------------------------------------------------

void timing_display(timing_t* t_p) {
  if (t_p == NULL) return;

  printf("\n");
  printf("===========================================================\n");
  printf("=                    T i m i n g                          =\n");
  printf("===========================================================\n");
  int i, j;
  for(i=0 ; i < t_p->num_sections ; i++) {
    if(i == t_p->num_sections - 1) {
      printf("-----------------------------------------------------------\n");
    }
    t_p->section_times_p[i] = 0;
    for(j=0 ; j < t_p->num_threads_p[i] ; j++) {
      if (t_p->section_times_p[i] < t_p->thread_times_p[i][j]) {
	t_p->section_times_p[i] = t_p->thread_times_p[i][j];
      }
    }
    if (t_p->num_threads_p[i] > 1) {
      printf("%s\t: %4.04f sec (maximum time)\n", t_p->section_labels_p[i], t_p->section_times_p[i] / 1000000);
      for(j=0 ; j < t_p->num_threads_p[i] ; j++) {
	printf("\t%i : %4.04f sec\n", j, t_p->thread_times_p[i][j] / 1000000);
      }
    } else {
      printf("%s\t: %4.04f sec\n", t_p->section_labels_p[i], t_p->section_times_p[i] / 1000000);
    }
  }
  printf("===========================================================\n");
  printf("===========================================================\n");
  printf("\n");
}

//-----------------------------------------------------
//-----------------------------------------------------
