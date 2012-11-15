#include <stdio.h>
#include <stdlib.h>

#include "bwt.h"

int main(int argc, char *argv[]) {

  if (argc != 4) {
    printf("Error.\n");
    printf("Usage: %s index-dirname seq num-errors\n", argv[0]);
    exit(-1);
  }

  char *index_dirname = argv[1];
  char *seq = argv[2];
  int num_errors = atoi(argv[3]);

  // initializations
  initReplaceTable();

  bwt_optarg_t *bwt_optarg = bwt_optarg_new(num_errors, 1, 10000, 1, 0, 0);
  bwt_index_t *bwt_index = bwt_index_new(index_dirname);

  // seq
  {
    array_list_t *mapping_list = array_list_new(100000, 1.25f, 
						COLLECTION_MODE_SYNCHRONIZED);
    
    size_t num_mappings;
    
    num_mappings = bwt_map_seq(seq, bwt_optarg, 
			       bwt_index, mapping_list);
    printf("seq: %s\n", seq);
    printf("num_mappings = %lu\n", num_mappings);
    for (size_t i = 0; i < num_mappings; i++) {
      printf("%lu\t---------------------\n", i);
      alignment_print(array_list_get(i, mapping_list));
    }
  }


  // seed
  {
    array_list_t *mapping_list = array_list_new(100000, 1.25f, 
						COLLECTION_MODE_SYNCHRONIZED);
    
    size_t num_mappings;
    
    size_t len = strlen(seq);
    char *code_seq = (char *) calloc(len + 10, sizeof(char));
    replaceBases(seq, code_seq, len);
    
    num_mappings = bwt_map_exact_seeds_seq(code_seq, 18, 16,
					   bwt_optarg, bwt_index, mapping_list);
    
    region_t *region;
    for (size_t i = 0; i < num_mappings; i++) {
      region = array_list_get(i, mapping_list);
      printf("Region: chr = %lu, strand = %d, start = %lu, end = %lu\n", 
	     region->chromosome_id, region->strand, region->start, region->end);
    }
  }

  printf("Done.\n");

}
