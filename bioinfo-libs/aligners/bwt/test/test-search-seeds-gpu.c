#include <stdio.h>
#include <stdlib.h>

#include "bwt.h"
#include "bwt_gpu.h"
#include "gpu.h"

#include "bioformats/fastq/fastq_file.h"
#include "bioformats/bam-sam/alignment.h"

#include "bioformats/fastq/fastq_batch.h"

#include "common-libs/commons/commons.h"

void main(int argc, char *argv[]) {

  if (argc != 6) {
    printf("Error.\n");
    printf("Usage: %s index-dirname file-name num-errors batch-size gpu-threads\n", argv[0]);
    exit(-1);
  }

  char *index_dirname = argv[1];
  char *file_name = argv[2];
  int num_errors = atoi(argv[3]);
  size_t batch_size = atoi(argv[4]);
  int num_gpu_threads = atoi(argv[5]);

  // initializations
  initReplaceTable();

  bwt_optarg_t *bwt_optarg = bwt_optarg_new(num_errors, 1, 10000, 1, 0, 0);
  printf("Loading genome...\n");
  bwt_index_t *bwt_index = bwt_index_new(index_dirname);
  printf("Loading done!\n");

  gpu_context_t *context = gpu_context_new(0, num_gpu_threads, bwt_index);    
  cal_optarg_t *cal_optarg = cal_optarg_new(60, 60, 18, 18, 0);

  size_t num_mappings;
  fastq_batch_t *batch, *unmapped_batch;
  int num_reads = -1;
  size_t total_reads = 0;
  size_t num_batches = 0;
  array_list_t **mapping_list;
  size_t reads_unmapped = 0;
  fastq_file_t *file = fastq_fopen(file_name);
  
  while(1) {
    
    batch = fastq_batch_new(batch_size);
    unmapped_batch = fastq_batch_new(batch_size);
    
    // read reads from file #1
    if (num_reads != 0) {
      //printf("Loading batch from file...\n");
      num_reads = fastq_fread_batch_max_size(batch, batch_size, file);
      //printf("Load done!\n");
    }else {
      break;
    }
    
    if (num_reads > 0) {
      num_batches++;
      total_reads += num_reads;
      //Process bwt batch
      //printf("Process BWT batch %i ...\n", num_batches);
      //bwt_map_seq(seq, bwt_optarg, index, mapping_list);   
      mapping_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
      for (size_t i = 0; i < num_reads; i++) {
	mapping_list[i] = array_list_new(1000, 
					 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED);
      }

      num_mappings = bwt_map_exact_seed_batch_gpu(batch,
						  bwt_optarg, 
						  cal_optarg,
						  bwt_index,
						  context,
						  mapping_list);
      
      //printf("Process BWT done!\n");
      // allocationg memory for the current fastq batch
      batch = fastq_batch_new(batch_size);
      
    
      reads_unmapped += unmapped_batch->num_reads;
      
      fastq_batch_free(batch);
      fastq_batch_free(unmapped_batch);
      //show results
      printf("===================================\n");
      printf("num_mappings = %d\n", num_mappings);
      
      //printf("seq: %s\n", seq);
 
      for (size_t i = 0; i < num_reads; i++) {
	num_mappings = array_list_size(mapping_list[i]);
	printf("Read %i (%i):\n", i, num_mappings);
	for (size_t j = 0; j < num_mappings; j++) {
	  printf("\tmapping: %i\n", j);
	  region_t *region = array_list_get(j, mapping_list[i]);
	  printf("\t@DNA chromosome:%i | strand:%i | position:%i\n", 
		 region->chromosome_id,
		 region->strand,
		 region->start);
	}
      }

      printf("===================================\n\n");
      /*printf("%s\n", alignment->sequence);
	printf("+\n");
	printf("%s\n", alignment->quality);
	}
      */
      if (mapping_list) {
	for (size_t i = 0; i < num_reads; i++) {
	  array_list_free(mapping_list[i], region_free);
	}
	free(mapping_list);
      }else {
	break;
      }

    }
  
  /*printf("Reads Process %i\n", total_reads);
  printf("Reads Unmapped %i\n", reads_unmapped);
  printf("Reads Mapped %i\n", total_reads - reads_unmapped);

  printf("Done.\n");
  */
  }
}
