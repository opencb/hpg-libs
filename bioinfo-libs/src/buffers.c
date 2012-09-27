
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "buffers.h"

#define MAXLINE 2048

//====================================================================================
//  write_batch functions
//====================================================================================

write_batch_t* write_batch_new(unsigned int allocate_size, unsigned char flag) {
  write_batch_t* write_batch_p = (write_batch_t*) calloc(1, sizeof(write_batch_t));
  
  write_batch_p->flag = flag;
  write_batch_p->size = 0;
  
  if(flag != MATCH_FLAG){
    write_batch_p->allocated_size = allocate_size;
    write_batch_p->buffer_p = (void *) calloc(write_batch_p->allocated_size, sizeof(char));
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
