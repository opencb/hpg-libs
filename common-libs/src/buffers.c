
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "buffers.h"

#define MAXLINE 2048

//====================================================================================

void region_batch_init(array_list_t **allocate_mapping_p, fastq_batch_t *unmapped_batch_p, region_batch_t *region_batch_p){
  region_batch_p->allocate_mapping_p = allocate_mapping_p;
  region_batch_p->unmapped_batch_p = unmapped_batch_p;
}

void region_batch_free(region_batch_t *region_batch_p){
  for(int i = 0; i < region_batch_p->unmapped_batch_p->num_reads; i++){
    array_list_free(region_batch_p->allocate_mapping_p[i], region_free);
  }
  free(region_batch_p->allocate_mapping_p);
  fastq_batch_free(region_batch_p->unmapped_batch_p);
  free(region_batch_p);
  
}

//====================================================================================

void sw_batch_init(unsigned int num_reads, array_list_t **allocate_cals_p, fastq_read_t **allocate_reads_p, sw_batch_t *sw_batch_p){
  sw_batch_p->num_reads = num_reads;
  sw_batch_p->allocate_reads_p = allocate_reads_p;
  sw_batch_p->allocate_cals_p = allocate_cals_p;
}

void sw_batch_free(sw_batch_t *sw_batch_p){
  
  for(int i = 0; i < sw_batch_p->num_reads; i++){
    array_list_free(sw_batch_p->allocate_cals_p[i], cal_free);
    fastq_read_free(sw_batch_p->allocate_reads_p[i]);
    //fastq_read_free(sw_batch_p->allocate_reads_p[i]);
  }

  free(sw_batch_p->allocate_cals_p);
  free(sw_batch_p->allocate_reads_p);
  free(sw_batch_p);
  
}

//====================================================================================
//  write_batch functions
//====================================================================================

write_batch_t* write_batch_new(unsigned int allocate_size, unsigned char flag) {
  write_batch_t* write_batch_p = (write_batch_t*) calloc(1, sizeof(write_batch_t));
  
  write_batch_p->flag = flag;
  write_batch_p->size = 0;
  
  if(flag != MATCH_FLAG){
    write_batch_p->allocated_size = allocate_size;
    write_batch_p->buffer_p = (void *) calloc(allocated_size, sizeof(char));
  }else{
    write_batch_p->allocated_size = allocate_size/sizeof(alignment_t *);
    write_batch_p->buffer_p = (void *) calloc(write_batch_p->allocated_size, sizeof(alignment_t *));
  }
  return write_batch_p;
}

//------------------------------------------------------------------------------------

void write_batch_free(write_batch_t* write_batch_p) {
 if (write_batch_p == NULL) return;
 
 if (write_batch_p->buffer_p != NULL) free(write_batch_p->buffer_p);
 
 free(write_batch_p);
}

//====================================================================================

pair_mng_t *pair_mng_new(int pair_mode, size_t min_distance, size_t max_distance) {
  pair_mng_t *p = (pair_mng_t*) calloc(1, sizeof(pair_mng_t));

  p->pair_mode = pair_mode;
  p->min_distance = min_distance;
  p->max_distance = max_distance;

  return p;
}

//------------------------------------------------------------------------------------

void pair_mng_free(pair_mng_t *p) {
  if (p != NULL)
    free(p);
}

//====================================================================================

cal_batch_t* cal_batch_new(array_list_t **allocate_mapping_p, fastq_batch_t *unmapped_batch_p){
  cal_batch_t* cal_batch_p = (cal_batch_t *)malloc(sizeof(cal_batch_t));
  
  cal_batch_p->allocate_mapping_p = allocate_mapping_p;
  cal_batch_p->unmapped_batch_p = unmapped_batch_p;
  
  return cal_batch_p;
}

void cal_batch_free(cal_batch_t *cal_batch_p){
  for(int i = 0; i < cal_batch_p->unmapped_batch_p->num_reads; i++){
    array_list_free(cal_batch_p->allocate_mapping_p[i], region_free);
  }
  free(cal_batch_p->allocate_mapping_p);
  fastq_batch_free(cal_batch_p->unmapped_batch_p);
  free(cal_batch_p);
  
}

//====================================================================================

sw_batch_t* sw_batch_new(unsigned int num_reads, array_list_t **allocate_cals_p, fastq_read_t **allocate_reads_p){
  sw_batch_t* sw_batch_p = (sw_batch_t *)malloc(sizeof(sw_batch_t));

  sw_batch_p->num_reads = num_reads;
  sw_batch_p->allocate_reads_p = allocate_reads_p;
  sw_batch_p->allocate_cals_p = allocate_cals_p;
  
  return sw_batch_p;
}

void sw_batch_free(sw_batch_t *sw_batch_p){
  
  for(int i = 0; i < sw_batch_p->num_reads; i++){
    array_list_free(sw_batch_p->allocate_cals_p[i], cal_free);
    fastq_read_free(sw_batch_p->allocate_reads_p[i]);
  }

  free(sw_batch_p->allocate_cals_p);
  free(sw_batch_p->allocate_reads_p);
  free(sw_batch_p);
}

//-----------------------------------------------------------------------------------------

unsigned int pack_junction(unsigned int chromosome, unsigned int strand, unsigned int start, unsigned int end, unsigned int junction_id, unsigned int num_reads, char* buffer_p){
  int len;
  char str[1024];
  char *chr_p, *p = buffer_p;
  char strand_char[2] = {'+', '-'};

  if (chromosome == 23) { sprintf(str, "%c\0", 'X'); }
  else if (chromosome == 24) { sprintf(str, "%c\0", 'Y'); }
  else if (chromosome == 25) { sprintf(str, "%s\0", "MT"); }
  else { sprintf(str, "%i\0", chromosome); }
 
  len = strlen(str);
  memcpy(p, str, len);
  p += len;
  *p = '\t';
  p++;
  
  sprintf(str, "%i", start);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  sprintf(str, "%i", end);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  sprintf(str, "JUNCTION_%i", junction_id);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  sprintf(str, "%i", num_reads);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  *p = strand_char[strand]; 
  p++;
  *p = '\n'; 
  p++;
  

  return (p - buffer_p);
}

//=====================================================================================
//=====================================================================================

aligner_batch_t *aligner_batch_new(fastq_batch_t *fq_batch) {
  aligner_batch_t *p = (aligner_batch_t *) calloc(1, sizeof(aligner_batch_t));

  size_t num_reads = fq_batch->num_reads;

  p->action = BWT_ACTION;
  p->all_targets = 1;
  p->num_targets = 0;
  p->num_allocated_targets = num_reads;
  p->num_mapping_lists = num_reads;

  p->num_done = 0;
  p->num_to_do = 0;

  p->fq_batch = fq_batch;
  p->targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(500, 
					 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED); 
  }
    
  return p;
}

//------------------------------------------------------------------------------------

void aligner_batch_free(aligner_batch_t *p) {
  if (p == NULL) return;
  
  if (p->fq_batch != NULL) fastq_batch_free(p->fq_batch);
  if (p->targets != NULL) free(p->targets);
  if (p->mapping_lists != NULL) { free(p->mapping_lists); }
  
  free(p);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
