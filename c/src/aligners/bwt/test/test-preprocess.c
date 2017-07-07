#include "bwt.h"

int main(int argc, char *argv[]) {

  if (argc != 4) {
    printf("Error.\n");
    printf("Usage: %s ref-filename output-dir s-ratio\n", argv[0]);
    exit(-1);
  }

  char *ref_filename = argv[1];
  char *out_dirname = argv[2];
  unsigned int s_ratio = atoi(argv[3]);

  // initializations
  initReplaceTable();

  bwt_generate_index_files(ref_filename, out_dirname, s_ratio);

  printf("Done.\n");
}
