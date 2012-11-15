#include <stdio.h>
#include <stdlib.h>

#include "bwt.h"
//#include "bioformats/fastq/fastq_batch_reader.h"
#include "bioformats/fastq/fastq_batch.h"
#include "bioformats/fastq/fastq_file.h"

#include "common-libs/commons/commons.h"

void main(int argc, char *argv[]) {

  if (argc != 5) {
    printf("Error.\n");
    printf("Usage: %s index-dirname file-name num-errors batch-size\n", argv[0]);
    exit(-1);
  }

  char *index_dirname = argv[1];
  char *file_name = argv[2];
  int num_errors = atoi(argv[3]);
  size_t batch_size = atoi(argv[4]);

  // initializations
  initReplaceTable();

  bwt_optarg_t *bwt_optarg = bwt_optarg_new(num_errors, 1, 10000, 1, 0, 0);
  //printf("Loading genome...\n");
  bwt_index_t *bwt_index = bwt_index_new(index_dirname);
  //printf("Loading done!\n");
    

  size_t num_mappings;
  fastq_batch_t *batch, *unmapped_batch;
  int num_reads = -1;
  size_t total_reads = 0;
  size_t num_batches = 0;
  array_list_t *mapping_list;
  size_t reads_unmapped = 0;
  fastq_file_t *file = fastq_fopen(file_name);
  
  //printf("Process File: %s\n", file_name);
  // seq
  /*
  {
    array_list_t *mapping_list = array_list_new(100000, 1.25f,
                                                COLLECTION_MODE_SYNCHRONIZED);

    size_t num_mappings;
    
    num_mappings = bwt_map_seq(seq, bwt_optarg,
                               bwt_index, mapping_list);
    printf("seq: %s\n", seq);
    printf("num_mappings = %d\n", num_mappings);
    for (size_t i = 0; i < num_mappings; i++) {
      printf("%i\t---------------------\n", i);
      alignment_print(array_list_get(i, mapping_list));
    }
    exit(-1);
  }
 */
  while(1) {
    
    batch = fastq_batch_new(batch_size);
    unmapped_batch = fastq_batch_new(batch_size);
    mapping_list = array_list_new(100000, 1.25f, 
				  COLLECTION_MODE_SYNCHRONIZED);
    
    // read reads from file #1
    if (num_reads != 0) {
      //printf("Loading batch from file...\n");
      num_reads = fastq_fread_batch_max_size(batch, batch_size, file);
      //printf("Load done!\n");
    }
    
    if (num_reads > 0) {
      num_batches++;
      total_reads += num_reads;
      //Process bwt batch
      //printf("Process BWT batch %i ...\n", num_batches);
      //bwt_map_seq(seq, bwt_optarg, index, mapping_list);   
      num_mappings = bwt_map_exact_batch(batch,
					 bwt_optarg, 
					 bwt_index, 
					 unmapped_batch,
					 mapping_list);
      //printf("Process BWT done!\n");
      // allocationg memory for the current fastq batch
      batch = fastq_batch_new(batch_size);
    }
    
    reads_unmapped += unmapped_batch->num_reads;
    
    fastq_batch_free(batch);
    fastq_batch_free(unmapped_batch);
    //show results
    //printf("num_mappings = %d\n", num_mappings);
        
    //	  printf("seq: %s\n", seq);
    /*
    for (size_t i = 0; i < num_mappings; i++) {
      //   printf("%i\t---------------------\n", i);
      alignment_t *alignment = array_list_get(i, mapping_list);
      /*printf("@DNA chromosome:%i | strand:%i | position:%i\n", 
	     alignment->chromosome,
	     alignment->seq_strand,
	     alignment->position);
      printf("%s\n", alignment->sequence);
      printf("+\n");
      printf("%s\n", alignment->quality);
    }
  */
    
    
    array_list_free(mapping_list, alignment_free);
    if (num_reads == 0) {
	break;
    }
  }
  
  /*printf("Reads Process %i\n", total_reads);
  printf("Reads Unmapped %i\n", reads_unmapped);
  printf("Reads Mapped %i\n", total_reads - reads_unmapped);

  printf("Done.\n");
  */
}
