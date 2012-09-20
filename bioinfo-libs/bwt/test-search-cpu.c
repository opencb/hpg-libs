#include "bwt.h"

void main(int argc, char *argv[]) {

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

  bwt_optarg_t *bwt_optarg = bwt_optarg_new(num_errors, 1, 10000);
  bwt_index_t *bwt_index = bwt_index_new(index_dirname);

  array_list_t *mapping_list = array_list_new(100000, 1.25f, 
					      COLLECTION_MODE_SYNCHRONIZED);

  size_t num_mappings = bwt_map_seq(seq, bwt_optarg, 
				    bwt_index, mapping_list);

  printf("seq: %s\n", seq);
  printf("num_mappings = %d\n", num_mappings);
  for (size_t i = 0; i < num_mappings; i++) {
    printf("%i\t---------------------\n", i);
    alignment_print(array_list_get(i, mapping_list));
  }
  
  printf("Done.\n");
}
