#include <stdio.h>
#include <stdlib.h>

#include "bwt.h"

int main(int argc, char *argv[]) {

  if (argc != 6) {
    printf("Error.\n");
    printf("Usage: %s index-dirname seq num-errors start-seed end-seed\n", argv[0]);
    exit(-1);
  }

  char *index_dirname = argv[1];
  char *seq = argv[2];
  int num_errors = atoi(argv[3]);
  int start_seed = atoi(argv[4]);
  int end_seed = atoi(argv[5]);

  // initializations
  initReplaceTable();

  size_t len = strlen(seq);

  bwt_optarg_t *bwt_optarg = bwt_optarg_new(num_errors, 1, 100000, 1, 0, 0);
  bwt_index_t *bwt_index = bwt_index_new(index_dirname);

  // seed
  {
    array_list_t *mapping_list = array_list_new(100000, 1.25f, 
						COLLECTION_MODE_SYNCHRONIZED);
    
    size_t num_mappings;
    
    char *code_seq = (char *) calloc(len + 10, sizeof(char));
    replaceBases(seq, code_seq, len);
    
    num_mappings = bwt_map_inexact_seed(code_seq, len, start_seed, end_seed,
					bwt_optarg, bwt_index, mapping_list);
    
    region_t *region;
    for (size_t i = 0; i < num_mappings; i++) {
      region = array_list_get(i, mapping_list);
      printf("Region: chr = %lu, strand = %d, start = %lu, end = %lu\n", 
	     region->chromosome_id, region->strand, region->start, region->end);
    }
  }



  // seq
  {
    char aux[len + 1];
    memset(aux, 0, len + 1);
    memcpy(aux, seq + start_seed, end_seed - start_seed + 1);

    alignment_t *alig;
    array_list_t *mapping_list = array_list_new(100000, 1.25f, 
						COLLECTION_MODE_SYNCHRONIZED);
    
    size_t num_mappings;
    
    num_mappings = bwt_map_seq(aux, bwt_optarg, 
			       bwt_index, mapping_list);
    printf("aux seq: %s\n", aux);
    printf("num_mappings = %lu\n", num_mappings);
    for (size_t i = 0; i < num_mappings; i++) {
      alig = array_list_get(i, mapping_list);
      printf("%lu\t---------------------\n", i);
      printf("\tstrand = %i, chromosome = %i, position = %i\n", 
	     alig->seq_strand, alig->chromosome, alig->position);
    }
  }


  printf("Done.\n");

}
