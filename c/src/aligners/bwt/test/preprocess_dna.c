#include <stdio.h>

#include "genome.h"

#define MIN_OPTIONS 3

int main(int argc, char** argv) {
  if (argc < MIN_OPTIONS) {
    printf("Usage: %s Dna_reference_Path Output_Filename\n", argv[0]);
    exit(-1);
  }
  char *dna_binary_filename = argv[2];
  char *dna_filename = argv[1];

  generate_codes(dna_binary_filename, dna_filename);
}
