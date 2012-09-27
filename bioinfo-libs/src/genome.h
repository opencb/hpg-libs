
#ifndef GENOME_H
#define GENOME_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct genome {
  unsigned long int nX;
  unsigned int num_chromosomes;
  unsigned int chr_name_length[200];
  unsigned int chr_size[200];
  unsigned int chr_offset[200];
  unsigned int chr_name[200][100];
  char* X;
} genome_t;

genome_t* genome_new(char* sequence_filename, char* chromosome_filename);
void genome_free(genome_t* genome_p);

void genome_read_sequence(char* sequence, unsigned int strand, char* chromosome, unsigned long int* start_p, unsigned long int* end_p, genome_t* genome_p);
void genome_read_sequence_by_chr_index(char* sequence, unsigned int strand, unsigned int chr, unsigned long int* start_p, unsigned long int* end_p, genome_t* genome_p);
char* genome_get_chr_name(unsigned int chr, unsigned int* len, genome_t* genome_p);


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif  // GENOME_H
