#include "test2.h"


void main(int argc, char *argv[]) {

  char *input_filename, *index_dirname, *fastq_filename;

  if (argc != 6) {
    printf("Error.\n");
    printf("Usage: %s input-filename index-dirname fastq-filename threads-number batch-size\n", argv[0]);
    exit(-1);
  }
  
  struct timeval start_time, end_time;
  double timer;
  unsigned int threads = atoi(argv[4]);
  input_filename = argv[1];
  index_dirname = argv[2];
  fastq_filename = argv[3];
  
  // initializations
  initReplaceTable();

  printf("input_filename = %s, index_dirname = %s\n", input_filename, index_dirname);

  bwt_index_t *bwt_index = bwt_index_new(index_dirname);

  printf("**** index_dirname = %s\n", bwt_index->dirname);

  //unsigned int n_reads = 8;
  array_list_t *mappings;
  unsigned int num_mappings, total_mappings = 0;
  
  bwt_optarg_t * bwt_optarg_p = bwt_optarg_new(1, threads, 50);
  /*char *seq = (char *)malloc(sizeof(char)*26);
  seq = "TTTATTCGCAAATGCATAACTTGAAA";
  fastq_read_t *fq_read;*/
  
  /*fq_read = fastq_read_new("@read_0 2L_11177498_11177561/1", 
			   "TTTATTCGCAAATGCATAACTTGAAA", 
			   "##########################");*/
  /*printf("BWT EXACT +\n");
  num_mappings = bwt_map_exact_seq_cpu(seq, bwt_optarg_p, bwt_index, mappings);
  printf("\n\nBWT INEXACT +\n");
  num_mappings = bwt_map_exact_read_cpu(fq_read, bwt_optarg_p, bwt_index, mappings);*/
  
  unsigned int batch_size = atoi(argv[5]); //30MB
  fastq_batch_t *fastq_batch_p;
  fastq_batch_t *unmapped_batch_p;
  
  unsigned int num_reads = 10;
  unsigned int count = 0, length_header, length_seq;
  char header[1024];
  char seq_err[1024]; //Correct read: TTTATTCGCAAATGCATAACTTGAAA
  char seq_ok[1024];
  char quality[1024];
  char seq[1024];

  list_t read_list;
  list_init("read", 1, 10, &read_list);
  list_item_t *item_p = NULL;
  
  printf("Batch Init...\n");
  fastq_batch_reader_input_t input;
  fastq_batch_reader_input_init(fastq_filename, batch_size, &read_list, &input);
  
  omp_set_nested(1);
  omp_set_num_threads(2);
  printf("Parallel Time\tSequential Time\tTotal Time\n");
  #pragma omp parallel sections
  {
    #pragma omp section
    {
      fastq_batch_reader(&input);
    }
    #pragma omp section
    {
      	while ( (item_p = list_remove_item(&read_list)) != NULL ) {
      		fastq_batch_p = (fastq_batch_t *)item_p->data_p;
      		//printf("Extract one item ...\n");
      		unmapped_batch_p = fastq_batch_new(batch_size);
      
      		//printf("Batch Init done. Processing batch ...\n");
      		mappings = array_list_new(100000, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
      		//start_timer(start_time);
      		num_mappings = bwt_map_inexact_batch_cpu(fastq_batch_p, bwt_optarg_p, bwt_index, unmapped_batch_p, mappings);
      		total_mappings += num_mappings;
      		//sleep(5);
      		//stop_timer(start_time, end_time, timer);
      
     		/*printf("Process batch end!\n\n");
     		printf("=========================Time============================\n");
     		printf("Total Reads Mapped and Process %d/%d: %4.06f secs\n", fastq_batch_p->num_reads - unmapped_batch_p->num_reads, fastq_batch_p->num_reads, timer / 1000000);
      		printf("=========================================================\n");
		*/
		
     		 //printf("Process batch done! (Total Reads=%d)/(Reads Mapped=%d) - Num Mappings=%d\n", , fastq_batch_p->num_reads - unmapped_batch_p->num_reads, num_mappings);
     		array_list_free(mappings, alignment_free);
     		fastq_batch_free(fastq_batch_p);
     		fastq_batch_free(unmapped_batch_p);
     		list_item_free(item_p);
      /*alignment_t *mapping;
      printf("==========================Output=============================\n");
      printf("Found %d solution:\n", num_mappings);
      for (int i = 0; i < num_mappings; i++) {
	mapping = (alignment_t *) array_list_get(i, mappings);
	printf("Header: %s | Sequence:%s | Strand:%i | chromosome:%i | Cigar:%s | Read(%d/%d)\n", mapping->query_name, mapping->sequence, mapping->seq_strand, mapping->chromosome, mapping->cigar, i, num_mappings);
	//printf("Header: %s\n", mapping->query_name);
      }
      printf("--------------------------------------------------------------\n");
      printf("Reads in unmmaped batch %d:\n", unmapped_batch_p->num_reads);
      for (int i = 0; i < unmapped_batch_p->num_reads; i++) {
	    printf("%s\n", &(unmapped_batch_p->header[unmapped_batch_p->header_indices[i]]));
	    printf("%s\n", &(unmapped_batch_p->seq[unmapped_batch_p->data_indices[i]]));
	    printf("+\n");
	    printf("%s\n\n", &(unmapped_batch_p->quality[unmapped_batch_p->data_indices[i]]));
      }
      printf("================================================================\n");*/
  	}
    }
  }
  printf("-------------------------------------------------\n");
  printf("%4.06f\t%4.06f\t%4.06f\n", global_parallel / 1000000, global_sequential / 1000000, (global_parallel + global_sequential)/ 1000000);
  printf("-------------------------------------------------\n");
  printf("Total Mappings = %d \n", total_mappings);
  bwt_index_free(bwt_index);
}
