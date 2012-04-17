
#ifndef MAP_ALIGNMENTS_H_
#define MAP_ALIGNMENTS_H_

#include <stdio.h>

#include "bam.h"
#include "bam_commons.h"
#include "commons.h"


//====================================================================================
//  map_alignments.h
//
//  bam_file structure and prototypes
//====================================================================================

#define MIN_ALLOCATED_SIZE_FOR_CIGAR_STRING		5

#define HUMAN			1	//Specie: Human

#define NCBI37			1	//Assembly: NCBI 3.7


//bam_record structure

typedef struct alignment {
  short int optional_fields_length;	//Length of optional fields
  short int chromosome;		//Chromosome
  int position;			//Mapping position 
  int mate_position;		//Mapping position of the mate
  short int mate_chromosome;	//Chromosome of the mate
  short int template_length;	//Template length
  int map_quality;		//Map quality
  int num_cigar_operations;	//number of CIGAR operations
  
  //flags
  uint8_t is_paired_end;	//0: single end, 1: paired end
  uint8_t is_paired_end_mapped;	//0: pair not mapped, 1: pair mapped
  uint8_t is_seq_mapped;	//0: seq not mapped, 1: seq mapped
  uint8_t is_mate_mapped;	//0: mate unmapped, 1, mate mapped
  uint8_t seq_strand;		//0 (forward) or 1 (reverse)
  uint8_t mate_strand;		//Strand of the mate read, same values
  uint8_t pair_num;		//1 or 2
  uint8_t primary_alignment;	//0: not primary, 1 primary
  uint8_t fails_quality_check;	//0: meet quality checks, 1: quality checks not meeted
  uint8_t pc_optical_duplicate;	//0: not duplicate, 1: ducplicate

  char* query_name;		//Query template name
  char* sequence;		//Sequence of nts
  char* quality;		//Quality of nts
  char* cigar;			//CIGAR string
  uint8_t* optional_fields;	//Optional fields
} alignment_t;

typedef struct alignment_batch {
  int type; // SINGLE or MULTIPLE CHROMOSOMES (alignments)
  int num_alignments;
  
  int allocated_alignments;
  alignment_t** alignments_p;
} alignment_batch_t;


//----------------------------------------------------------------------------------
//  map alignments functions
//----------------------------------------------------------------------------------

alignment_t* alignment_new();
void alignment_init_single_end(char* query_name, char* sequence, char* quality, short int strand, short int chromosome, int position, char* cigar, short int num_cigar_operations, short int map_quality, alignment_t* alignment_p);
void alignment_init_paired_end(char* query_name, char* sequence1, char* sequence2, char* quality1, char* quality2, short int strand1, short int strand2, short int chromosome1, int position1, int position2, short int chromosome2, char* cigar1, char* cigar2, short int num_cigar_operations1, short int num_cigar_operations2, short int map_quality1, short int map_quality2, alignment_t* alignment1_p, alignment_t* alignment2_p);
alignment_t* alignment_new_by_bam(bam1_t* bam_p, int base_quality);
void alignment_free(alignment_t* alignment_p);

bam1_t* convert_to_bam(alignment_t* alignment_p, int base_quality);

void alignment_print(alignment_t* aligment_p);
void bam_print(bam1_t* bam_p, int base_quality);

bam_header_t* bam_header_new(int specie, int assembly);

#endif /* MAP_ALIGNMENTS_ */

