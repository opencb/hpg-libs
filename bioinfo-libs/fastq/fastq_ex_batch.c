#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fastq_ex_batch.h"
#include "string_utils.h"

/*-----------------------------------------------------
 * constructor: allocate memory for the batch
 *----------------------------------------------------*/

fastq_ex_batch_t* fastq_ex_batch_new(unsigned long size) {

  // minimum expected sequence length is 35
  //
  static int min_read_size = 35;
  static int max_header_length = 128;

  fastq_ex_batch_t* batch_p = (fastq_ex_batch_t*) calloc(1, sizeof(fastq_ex_batch_t));

  batch_p->num_reads = 0;
  batch_p->data_size = size / 2; // seq + quality = size
  batch_p->data_indices_size = sizeof(int) * batch_p->data_size / min_read_size;

  batch_p->header_indices = (int*) calloc(1, batch_p->data_indices_size);
  batch_p->header = (char*) calloc(1, max_header_length * batch_p->data_size / min_read_size);
 
  batch_p->data_indices = (int*) calloc(1, batch_p->data_indices_size);
  batch_p->seq = (char*) calloc(1, batch_p->data_size);
  batch_p->quality = (char*) calloc(1, batch_p->data_size);

  return batch_p;
}

/*-----------------------------------------------------
 * destructor: free memory
 *----------------------------------------------------*/

void fastq_ex_batch_free(fastq_ex_batch_t* batch_p) {
  if (batch_p == NULL) return;

  if (batch_p->header_indices != NULL) free(batch_p->header_indices);
  if (batch_p->header != NULL)         free(batch_p->header);

  if (batch_p->data_indices != NULL)   free(batch_p->data_indices);
  if (batch_p->seq != NULL)            free(batch_p->seq);
  if (batch_p->quality != NULL)        free(batch_p->quality);

  free(batch_p);
}

/*-----------------------------------------------------
 * read batch from file
 *----------------------------------------------------*/

int fastq_ex_batch_read(fastq_ex_batch_t* batch_p, int encode, fastq_file_t* file_p) {
  static int MAX_HEADER_LENGTH = 1024;
  static int MAX_SEQUENCE_LENGTH = 4096;

  unsigned long accumulated_size = 0;

  int* p;

  char header1[MAX_HEADER_LENGTH];
  char sequence[MAX_SEQUENCE_LENGTH];
  char header2[MAX_HEADER_LENGTH];
  char quality[MAX_SEQUENCE_LENGTH];
  int header_length, sequence_length, quality_length;

  register unsigned long max_size = batch_p->data_size;
  register int count = 0;

  batch_p->header_indices[count] = 0;
  batch_p->data_indices[count] = 0;
  
  while (accumulated_size  <= (max_size-1024)                  && 
	 count*sizeof(int) <= (batch_p->data_indices_size-100) &&
	 fgets(header1, MAX_HEADER_LENGTH, file_p->fd) != NULL  ) {
    
    fgets(sequence, MAX_SEQUENCE_LENGTH, file_p->fd);
    fgets(header2, MAX_HEADER_LENGTH, file_p->fd);
    fgets(quality, MAX_SEQUENCE_LENGTH, file_p->fd);
    
    header_length = strlen(header1);
    sequence_length = strlen(sequence);
    quality_length = strlen(quality);
    
    // remove '\n' character, now length includes '\0' character
    //
    chomp_at(header1, header_length - 1);
    chomp_at(sequence, sequence_length - 1);
    chomp_at(quality, quality_length - 1);

    if (sequence_length != quality_length) {
      printf("ERROR: sequence length (%i) and quality length (%i) do not match, read: %s\n", sequence_length, quality_length, header1);
      exit(-1);
    }
    
    count++;
    
    //	  printf("******** a = (%i)\n", table['a']);
    //	  printf("******** c = (%i)\n", table['c']);
    
    //printf("header1 = %s\n", header1);
    //encodeBases(qualities, sequence, sequence_length);
    if (encode) {
      encodeBases(&(batch_p->seq[batch_p->data_indices[count-1]]), sequence, sequence_length);
    } else {
      memcpy(&(batch_p->seq[batch_p->data_indices[count-1]]), sequence, sequence_length);
    }
    //int i=0;
    //for(i=0 ; i<sequence_length ; i++) {
    //	    printf("d(%i)", (int)batch_p->data[batch_p->data_indices[count-1]+i]);
    //	    printf("q(%c)",qualities[i]);
    //	  }
    //	  printf("\n");
    
    //memcpy(&(batch_p->data[batch_p->data_indices[count-1]]), sequence, sequence_length);
    memcpy(&(batch_p->quality[batch_p->data_indices[count-1]]), quality, quality_length);
    memcpy(&(batch_p->header[batch_p->header_indices[count-1]]), header1, header_length);
    //printf("sequence = %s\n", sequence);
    
    //printf("reader, header = %s\n", &(batch_p->header[batch_p->header_indices[count-1]]));
    //printf("header = %s\n", &(batch_p->header[batch_p->header_indices[count-1]]));
    
    if (count*sizeof(int) >= batch_p->data_indices_size) {			  
      
      printf("reallocate %i >= %i, accumulated size = %i, count = %i\n", (count*sizeof(int)), batch_p->data_indices_size, accumulated_size, count);
      
      // maybe we should use the realloc function ??
      //
      int size = (count + 100) * sizeof(int);
      
      // that's for data indices...
      //
      printf("000000, reallocating size %i\n", size);
      p = (int*) calloc(size, sizeof(int));
      printf("111111\n");
      memcpy((void*) p, (void*)batch_p->data_indices, count * sizeof(int));
      printf("222222\n");
      
      free(batch_p->data_indices);
      batch_p->data_indices = p;
      
      batch_p->data_indices_size = size;
    }

    batch_p->data_indices[count] = batch_p->data_indices[count-1] + sequence_length;		
    batch_p->header_indices[count] = batch_p->header_indices[count-1] + header_length;		
    accumulated_size += sequence_length;
  }
  
  batch_p->num_reads = count;
  
  return batch_p->num_reads;
}


/*-----------------------------------------------------
 * write batch to file
 *----------------------------------------------------*/

int fastq_ex_batch_write(fastq_ex_batch_t* batch_p, int decode, fastq_file_t* file_p) {
  static int MAX_SEQUENCE_LENGTH = 4096;
  register int length, i = 0;

  char seq[MAX_SEQUENCE_LENGTH];
  
  printf("fastq_write_batch_ex\n");
  
  for(i=0 ; i<batch_p->num_reads ; i++) {
    
    printf("write, header = %s\n", &(batch_p->header[batch_p->header_indices[i]]));
    
    // writing header...
    //
    length = batch_p->header_indices[i+1] - batch_p->header_indices[i];
    fwrite(&(batch_p->header[batch_p->header_indices[i]]), 1, length, file_p->fd);
    fprintf(file_p->fd, "\n");
    
    // writing sequence...
    //
    length = batch_p->data_indices[i+1] - batch_p->data_indices[i];
    if (decode) {
      decodeBases(seq, &(batch_p->seq[batch_p->data_indices[i]]), length);
      fwrite(seq, 1, length, file_p->fd);
    } else {
      fwrite(&(batch_p->seq[batch_p->data_indices[i]]), 1, length, file_p->fd);
    }
    fprintf(file_p->fd, "\n+\n");
    
    // writing quality...
    //
    fwrite(&(batch_p->quality[batch_p->data_indices[i]]), 1, length, file_p->fd);
    fprintf(file_p->fd, "\n");;

    /*
    seq[length] = '\0';
    printf("write, seq = %s, length = %i\n", seq, length);
    fprintf(file_p->fd, "%s\n", seq);
    printf("%s\n", seq);
    */
  }
  return batch_p->num_reads;
}











/*
void fastq_ex_batch_print(fastq_ex_batch_t* fastq_ex_batch_p, FILE* fd) {
    int i;
  
    for(i=0; i<fastq_ex_batch_p->num_reads; i++) {
	fprintf(fd, "%s\n", &(fastq_ex_batch_p->header[fastq_ex_batch_p->header_indices[i]]));
	fprintf(fd, "%s\n", &(fastq_ex_batch_p->data[fastq_ex_batch_p->data_indices[i]]));
	fprintf(fd, "+\n");
	fprintf(fd, "%s\n", &(fastq_ex_batch_p->data[(fastq_ex_batch_p->data_indices[i]+fastq_ex_batch_p->data_indices[i+1])/2]));
    }
}

void fprintf_read(FILE* fd, fastq_ex_batch_t* batch_p, int index) {
    fprintf(fd, "%s\n", &(batch_p->header[batch_p->header_indices[index]]));
    fprintf(fd, "%s\n", &(batch_p->data[batch_p->data_indices[index]]));
    fprintf(fd, "+\n");
    fprintf(fd, "%s\n", &(batch_p->data[(batch_p->data_indices[index]+batch_p->data_indices[index+1])/2]));
}

void fprintf_rtrim_read(FILE* fd, fastq_ex_batch_t* batch_p, int index, int rtrim_length) {
    fprintf(fd, "%s\n", &(batch_p->header[batch_p->header_indices[index]]));
    fprintf(fd, "%s\n", rtrim(&(batch_p->data[batch_p->data_indices[index]]), rtrim_length));
    fprintf(fd, "+\n");
    fprintf(fd, "%s\n", rtrim(&(batch_p->data[(batch_p->data_indices[index]+batch_p->data_indices[index+1])/2]), rtrim_length));
}

void fprintf_ltrim_read(FILE* fd, fastq_ex_batch_t* batch_p, int index, int ltrim_length) {
    fprintf(fd, "%s\n", &(batch_p->header[batch_p->header_indices[index]]));
    fprintf(fd, "%s\n", ltrim(&(batch_p->data[batch_p->data_indices[index]]), ltrim_length));
    fprintf(fd, "+\n");
    fprintf(fd, "%s\n", ltrim(&(batch_p->data[(batch_p->data_indices[index]+batch_p->data_indices[index+1])/2]), ltrim_length));
}

void fprintf_trim_read(FILE* fd, fastq_ex_batch_t* batch_p, int index, int rtrim_length, int ltrim_length) {
    char* str;
    fprintf(fd, "%s\n", &(batch_p->header[batch_p->header_indices[index]]));
    //fprintf(fd, "%s\n", cut_last_nts(&(batch_p->data[batch_p->data_indices[index]]), trim_length));
    str = &(batch_p->data[batch_p->data_indices[index]]);
    if (rtrim_length > 0) {  str = rtrim(str, rtrim_length); }
    if (ltrim_length > 0) {  str = ltrim(str, ltrim_length); }
    fprintf(fd, "%s\n", str);
    fprintf(fd, "+\n");
    //fprintf(fd, "%s\n", cut_last_nts(&(batch_p->data[(batch_p->data_indices[index]+batch_p->data_indices[index+1])/2]), trim_length));
    str = &(batch_p->data[(batch_p->data_indices[index]+batch_p->data_indices[index+1])/2]);
    if (rtrim_length > 0) {  str = rtrim(str, rtrim_length); }
    if (ltrim_length > 0) {  str = ltrim(str, ltrim_length); }
    fprintf(fd, "%s\n", str);
}
*/


/*-----------------------------------------------------
 *----------------------------------------------------*/
