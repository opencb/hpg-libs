#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <omp.h>

#include "genome.h"

//===============================================================
//  main
//===============================================================

int main(int argc, char* argv[]) {

  // chromosome
  // start
  // end
  // genome-filename     = "/scratch2/jtarraga/appl/ngs-gpu/mapping-pipe/BWT";
  // chromosome-filename = "/scratch2/jtarraga/appl/ngs-gpu/mapping-pipe/sort.exome.index";

  printf("start main: %s chromosome start end genome-filename chromosome-filename\n", argv[0]);
  if (argc != 6) {
    printf("Usage: %s chromosome start end genome-filename chromosome-filename\n", argv[0]);
    printf("start main: done !!\n");
    exit(-1);
  }

  int chromosome = 1;
  long unsigned int start = 3000000L, end = 3000050L;

  sscanf(argv[1], "%i", &chromosome);
  sscanf(argv[2], "%lu", &start);
  sscanf(argv[3], "%lu", &end);
  char* genome_filename     = argv[4];
  char* chromosome_filename = argv[5];

  char* sequence = (char *) calloc(1, end - start + 1);
  genome_t* genome_p = genome_new(genome_filename, chromosome_filename);
  genome_read_sequence_by_chr_index(sequence, chromosome, &start, &end, genome_p);

  printf("%i:%li-%li, %s\n", chromosome, start, end, sequence);
}

//-----------------------------------------------------
//-----------------------------------------------------
