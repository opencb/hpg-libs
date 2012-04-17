
#ifndef BAM_FILE_H_
#define BAM_FILE_H_

#include <stdio.h>

#ifdef THRUST-GPU
  #include <thrust/host_vector.h>
  #include <thrust/device_vector.h>
  #include <thrust/sort.h>
  #include <thrust/copy.h>
#endif

#include "alignment.h"
#include "bam.h"
#include "bam_commons.h"
#include "commons.h"
#include "file_utils.h"
#include "string_utils.h"
#include "system_utils.h"


//====================================================================================
//  bam_file.h
//
//  bam_file structure and prototypes
//====================================================================================

#define MAX_NUM_PRODUCERS		10
#define BAM_BATCH_READ_SIZE			50000000
#define BAM_BATCH_WRITE_SIZE			50000000

#define SINGLE_CHROM_BATCH	0
#define MULTIPLE_CHROM_BATCH	1

// bam_file_t and bam_batch_t structures
//
typedef struct bam_file {
  unsigned int num_alignments;
  char* filename;
  char* mode;
  bamFile bam_fd;
  bam_header_t* bam_header_p;
} bam_file_t;


typedef struct bam_batch {
  int type; // SINGLE or MULTIPLE CHROMOSOMES (alignments)
  int num_alignments;
  
  int allocated_alignments;
  bam1_t** alignments_p;
} bam_batch_t;

//----------------------------------------------------------------------------------
//  bam file functions
//----------------------------------------------------------------------------------

// open and close functions
//
bam_file_t* bam_fopen(char* filename);
bam_file_t* bam_fopen(char* filename, bam_header_t* bam_header_p, char* mode);
void bam_fclose(bam_file_t *bam_file);

// read functions
//
int bam_fread_max_size(bam_batch_t* batch_p, size_t batch_size, int base_quality, bam_file_t* bam_file_p);
int bam_fread_max_size_no_duplicates(bam_batch_t* batch_p, size_t batch_size, int base_quality, bam_file_t* bam_file_p, uint8_t* prev_seq, int* prev_seq_length, int* prev_seq_start_coordinate);
int bam_fread_max_size_by_chromosome(bam_batch_t* batch_p, size_t batch_size, int chromosome, bam_file_t* bam_file_p);

// write functions
//
void bam_fwrite_header(bam_header_t* bam_header_p, bam_file_t* bam_file_p);

int bam_fwrite(bam1_t* alignment_p, bam_file_t* bam_file_p);
int bam_fwrite_array(bam1_t** alignment_p, int length, bam_file_t* bam_file_p);
int bam_fwrite_batch(bam_batch_t* batch_p, bam_file_t* bam_file_p);

#ifdef THRUST-GPU
  int bam_fwrite_sorted_array(bam1_t** alignment_p, thrust::host_vector<int> indices, int length, bam_file_t* bam_file_p);
#else
  int bam_fwrite_sorted_array(bam1_t** alignment_p, int* indices, int length, bam_file_t* bam_file_p);
#endif

//high level alignment functions
int alignment_fwrite(alignment_t* alignment_p, bam_file_t* bam_file_p);
int alignment_fwrite_array(alignment_t** alignment_p, int length, bam_file_t* bam_file_p);
int alignment_fwrite_batch(bam_batch_t* batch_p, bam_file_t* bam_file_p);  

//validation of the header
int bam_validate_header(bam_file_t* bam_file_p);
  
unsigned int bam_fcount(bam_file_t* bam_file);


//----------------------------------------------------------------------------------
//  bam batch functions
//----------------------------------------------------------------------------------

bam_batch_t* bam_batch_new(size_t batch_size, int type);
void bam_batch_free(bam_batch_t* batch_p, int free_alignments);
void bam_batch_print(bam_batch_t* batch_p, FILE* fd);
int bam_batch_compare_seq(uint8_t* data1, int length_seq1, int start_seq1, uint8_t* data2, int length_seq2, int start_seq2);

//----------------------------------------------------------------------------------
//  bam1_t functions
//----------------------------------------------------------------------------------

void free_bam1(bam1_t** alignment_p, int num_alignments);
void print_bam1(bam1_t* alignment_p, FILE* fd);

#endif    // __BAM_FILE_H_