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

  /*
  printf("input_filename = %s, index_dirname = %s\n", input_filename, index_dirname);

  bwt_index_t *bwt_index = bwt_index_new(index_dirname);

  printf("**** index_dirname = %s\n", bwt_index->dirname);

  unsigned int n_reads = 8;

  array_list_t *cal_list_p = array_list_new(100000, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  //fastq_read_t *fq_read;
  unsigned int num_mappings;
  
  char *seq = (char *)malloc(sizeof(char )*66);
  seq = "TATTTATTCGTAAATGCATAACTTGAAAGGCTATTTTGTAGATATTAAATGCAACTAATTTGCC\0"; //chromosome 2L 11177498_11177561/1
  //     TATTTATTCGCAAATGCATAACTTGAAAGGTTATTTTGTAGATATTAATTGCAACTAATTTGCC
  //               X                   X                 X   
  
  //seq = "TATTTATTCGCAAATCGCATAACTTGAAAGGCTATTTTGTAGATATTAAATGCAACTAATTTGCC\0";
  //     TATTTATTCGCAAATCGCATAACTTGAAAGGCTATTTTGTAGATATTAAATGCAACTAATTTGCC  //chromosoma 2L 11177498_11177561/1
  
  //seq = "AAAATGGAGCAAATAAATAGATCGATCAGTGCAGAAATAAACTTTAAACCTAAATTTCTAGT\0";
  //     AAAATGTAGCAAATAAATAGATAGATCAGTGCAGAAATAAACTTAAAACCTAAATTTCTAGT    //2L:21677681-21677742/1
  //           X                X                     X                 
  bwt_optarg_t * bwt_optarg_p = bwt_optarg_new(1, 1, 50);
  cal_optarg_t * bwt_caloptarg_p = cal_optarg_new(50, 50, 18, 10);
 

  
  
  num_mappings = bwt_find_cals_from_seq(seq, bwt_optarg_p, bwt_index, bwt_caloptarg_p, cal_list_p);
  cal_t *cal_p;
  unsigned int cals_found = array_list_size(cal_list_p);
  printf("Found %d CALs:\n", cals_found);
  
  for (int i = 0; i < cals_found; i++) {
    cal_p = (cal_t *) array_list_get(i, cal_list_p);
    printf("CAL:%i - chromosome %i - strand %i - len %d - [%i-%i] \n", i, cal_p->chromosome_id, cal_p->strand, cal_p->end - cal_p->start, cal_p->start, cal_p->end);
  }
  
  
  
  /*
  fq_read = fastq_read_new("@read_4 2L_253017_253081/1",
			   "CTTGTCGAAGTACTGCCGCTTGGTTCGTTGGGATCCTCAGCCTCGCGCACCTTTTCGGCAATCT",
			   "################################################################");

  num_mappings = bwt_single_read(fq_read, NULL, bwt_index, NULL);
  printf("num_mappings %i found for %s\n", num_mappings, fq_read->id);
  fastq_read_free(fq_read);
  */
  /*


  fq_read = fastq_read_new("@read_5 2L_11177498_11177561/1",
			   "TATTTATTCGCAAATCGCATAACTTGAAAGGCTATTTTGTAGATATTAAATGCAACTAATTTGCC",
			   "#################################################################");

  */
  


  /*
  int num_reads;
  unsigned long batch_size = 64 * 1024;
  fastq_batch_t *fq_batch;

  // open FastQ files
  fastq_file_t *fq_file = fastq_fopen(input_filename);

  fastq_read_t read;


  while (1) {

    fastq_fread(
    //fq_batch = fastq_batch_new();
    //fastq_batch_init(fq_batch, batch_size);

    //    num_reads = fastq_fread_batch_max_size(fq_batch, batch_size, fq_file);

    printf("num reads = %i\n", num_reads);

    // if no reads, free memory and go out...
    if (num_reads == 0) {
      fastq_batch_free(fq_batch);
      break;
    }


    


    fastq_batch_free(fq_batch);
  }

  // close FastQ file and free memory
  fastq_fclose(fq_file);
  */

  /*
  bwt_index_free(bwt_index);
}
  */
