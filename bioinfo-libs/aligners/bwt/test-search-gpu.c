
#include "bwt.h"
#include "bwt_gpu.h"
#include "gpu.h"

#include "bioformats/fastq/fastq_file.h"
#include "bioformats/bam-sam/alignment.h"

extern double kl_time, s_time;

void main(int argc, char *argv[]) {

  printf("Starting %s..\n", argv[0]);

  if (argc <= 4) {
    printf("Error.\n");
    printf("Usage: %s index-dirname num-gpu-threads mode sequence/filename batch_size\n", argv[0]);
    printf("The 'mode' parameter must be 'seq' or 'file'\n");
    exit(-1);
  }
  char *index_dirname = argv[1];
  int num_gpu_threads = atoi(argv[2]);
  char *mode = argv[3];
  char *filename = NULL, *seq = NULL;
  unsigned int batch_size;

  if (strcmp(mode, "file") == 0) {
    filename = argv[4];
    batch_size = atoi(argv[5]);
  } else if (strcmp(mode, "seq") == 0) {
    seq = argv[4];
    batch_size = 1000;
  } else {
    printf("Error.\n");
    printf("Usage: %s index-dirname num-gpu-threads mode sequence/filename batch_size\n", argv[0]);
    printf("   - mode: it must be seq or file\n");
    exit(-1);
  }
  
  int num_errors = 0;
  
  // initializations
  initReplaceTable();
  bwt_optarg_t *bwt_optarg = bwt_optarg_new(num_errors, 1, 10000);

  printf("reading BWT index from %s...\n", index_dirname);
  bwt_index_t *bwt_index = bwt_index_new(index_dirname);
  printf("...done !!\n");

  gpu_context_t *context = gpu_context_new(0, num_gpu_threads, bwt_index);

  fastq_batch_t *batch = fastq_batch_new(batch_size);
  fastq_batch_t *unmapped_batch = fastq_batch_new(batch_size);

  size_t num_mappings = 0;
  array_list_t *mapping_list = array_list_new(100000, 1.25f, 
					      COLLECTION_MODE_SYNCHRONIZED);

  int num_reads, num_batches = 0;
  size_t total_reads = 0, total_unmapped_reads = 0;

  if (seq != NULL) {
    // searching for input sequence

    total_reads = 1;

    batch->num_reads = 1;
    batch->data_indices_size = 2 * sizeof(int);
    batch->data_size = strlen(seq) + 1;
    batch->data_indices = (int*) calloc(2, sizeof(int));
    batch->data_indices[0] = 0;
    batch->data_indices[1] = strlen(seq) + 1;
    batch->seq = seq;

  
    num_mappings = bwt_map_exact_batch_gpu(batch, bwt_optarg,
					   bwt_index, context,
					   unmapped_batch, mapping_list);

    total_unmapped_reads += unmapped_batch->num_reads;

    array_list_free(mapping_list, alignment_free);
    fastq_batch_free(batch);

  } else {
    // searching from file
    fastq_file_t *file = fastq_fopen(filename);
    while ((num_reads = fastq_fread_batch_max_size(batch, batch_size, file)) > 0) {

     
      num_mappings += bwt_map_exact_batch_gpu(batch, bwt_optarg,
					      bwt_index, context,
					      unmapped_batch, mapping_list);
      num_batches++;
      total_reads += num_reads;    
      total_unmapped_reads += unmapped_batch->num_reads;


      array_list_free(mapping_list, alignment_free);
      mapping_list = array_list_new(100000, 1.25f, 
				    COLLECTION_MODE_SYNCHRONIZED);

      fastq_batch_free(batch);
      batch = fastq_batch_new(batch_size);
    }
    
    fastq_fclose(file);
  }

  printf("Total batches : %lu (%0.2f reads/batch)\n", num_batches, 1.0 * total_reads / num_batches);
  printf("\n");
  printf("Total reads   : %lu\n", total_reads);
  printf("Mapped reads  : %lu (%0.2f %)\n", 
	 total_reads - total_unmapped_reads, (total_reads - total_unmapped_reads) * 100.0 / total_reads);
  printf("Unmapped reads: %lu (%0.2f %)\n", 
	 total_unmapped_reads, total_unmapped_reads * 100.0 / total_reads);

  printf("\n");
  printf("Total mappings: %lu\n", num_mappings); 


  printf("\n");
  printf("KL time: %0.2f s\n", kl_time / 1000000.0);
  printf("S time: %0.2f s\n", s_time / 1000000.0);


  /*
  for (size_t i = 0; i < num_mappings; i++) {
    printf("alignment #%i\n", i);
    alignment_print(array_list_get(i, mapping_list));
  }
  */

  // free memory
  gpu_context_free(context);
  bwt_index_free(bwt_index);
  bwt_optarg_free(bwt_optarg);

  printf("Done.\n");
}
