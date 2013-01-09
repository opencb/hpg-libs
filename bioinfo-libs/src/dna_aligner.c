#include "dna_aligner.h"


//--------------------------------------------------------------------
void run_dna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     pair_mng_t *pair_mng, options_t *options) {

  // lists to communicate/synchronize the different threads
  list_t read_list, write_list;
  list_init("read", 1, 24, &read_list);
  list_init("write", options->num_cpu_threads, 24, &write_list);
  
  printf("writers to read_list = %d\n", list_get_writers(&read_list));
  printf("writers to write_list = %d\n", list_get_writers(&write_list));

  omp_set_nested(1);  
  bwt_server_input_t bwt_input;
  bwt_server_input_init(NULL, 0, bwt_optarg, bwt_index, 
			NULL, 0, NULL, &bwt_input);
  
  region_seeker_input_t region_input;
  region_seeker_input_init(NULL, cal_optarg, bwt_optarg, 
			   bwt_index, NULL, 0, options->gpu_process, 
			   &region_input);

  cal_seeker_input_t cal_input;
  cal_seeker_input_init(NULL, cal_optarg, NULL, 0, 
			NULL, NULL, &cal_input);
  
  pair_server_input_t pair_input;
  pair_server_input_init(pair_mng, bwt_optarg->report_best, bwt_optarg->report_n_hits, 
			 bwt_optarg->report_all, NULL, NULL, NULL, &pair_input);
  
  sw_server_input_t sw_input;
  sw_server_input_init(NULL, NULL, 0, options->match, options->mismatch, 
		       options->gap_open, options->gap_extend, options->min_score, 
		       options->flank_length, genome, 0, 0, 0,  bwt_optarg, &sw_input);
   

  double t_total;
  struct timeval t1, t2;
  
  #pragma omp parallel sections num_threads(3) //options->num_cpu_threads)
  {
    printf("Principal Sections %d threads\n", omp_get_num_threads());
    
    // fastq batch reader
    #pragma omp section
    {
      fastq_batch_reader_input_t input;
      fastq_batch_reader_input_init(options->in_filename, options->in_filename2, 
				    options->pair_mode, options->batch_size, 
				    &read_list, &input);
      fastq_batch_reader_aligner(&input);
    }

    // batch aligner
    #pragma omp section
    {
      gettimeofday(&t1, NULL);
      #pragma omp parallel num_threads(options->num_cpu_threads)
      {
	batch_aligner_input_t input;
	batch_aligner_input_init(&read_list, &write_list,
				 &bwt_input, &region_input, &cal_input,
				 &pair_input, &sw_input, &input);
	
	batch_aligner(&input);
      }
      gettimeofday(&t2, NULL);
      printf("\n\naligner time (no IO): %0.5f sec\n\n", (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6);
    }

    // batch writer
    #pragma omp section
    {
      batch_writer_input_t input;
      batch_writer_input_init(options->output_filename, NULL, NULL, &write_list, genome, &input);
      batch_writer2(&input);
    }
  }

  if (statistics_on) {
    size_t total_item = 0;
    double max_time = 0, total_throughput = 0;
    printf("\nBWT time:\n");
    for (int i = 0; i < options->num_cpu_threads; i++) {
      printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f BWT/s (reads)\n", 
	     i, bwt_time[i] / 1e6, thr_batches[i], thr_bwt_items[i], 1e6 * thr_bwt_items[i] / bwt_time[i]);
      total_item += thr_bwt_items[i];
      total_throughput += (1e6 * thr_bwt_items[i] / bwt_time[i]);
      if (max_time < bwt_time[i]) max_time = bwt_time[i];
    }
    printf("\n\tTotal BWTs: %lu, Max time = %0.4f, Throughput = %0.2f BWT/s\n", total_item, max_time / 1e6, total_throughput);
    
    total_item = 0; max_time = 0; total_throughput = 0;
    printf("\nSeeding time:\n");
    for (int i = 0; i < options->num_cpu_threads; i++) {
      printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f BWT/s (seeds)\n", 
	     i, seeding_time[i] / 1e6, thr_batches[i], thr_seeding_items[i], 1e6 * thr_seeding_items[i] / seeding_time[i]);
      total_item += thr_seeding_items[i];
      total_throughput += (1e6 * thr_seeding_items[i] / seeding_time[i]);
      if (max_time < seeding_time[i]) max_time = seeding_time[i];
    }
    printf("\n\tTotal BWTs: %lu, Max time = %0.4f, Throughput = %0.2f BWT/s\n", total_item, max_time / 1e6, total_throughput);
    
    total_item = 0; max_time = 0; total_throughput = 0;
    printf("\nCAL time:\n");
    for (int i = 0; i < options->num_cpu_threads; i++) {
      printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f CAL/s)\n", 
	     i, cal_time[i] / 1e6, thr_batches[i], thr_cal_items[i], 1e6 * thr_cal_items[i] / cal_time[i]);
      total_item += thr_cal_items[i];
      total_throughput += (1e6 * thr_cal_items[i] / cal_time[i]);
      if (max_time < cal_time[i]) max_time = cal_time[i];
    }
    printf("\n\tTotal CALs: %lu, Max time = %0.4f, Throughput = %0.2f CAL/s\n", total_item, max_time / 1e6, total_throughput);
    
    total_item = 0; max_time = 0; total_throughput = 0;
    printf("\nSW time:\n");
    for (int i = 0; i < options->num_cpu_threads; i++) {
      printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f SW/s)\n", 
	     i, sw_time[i] / 1e6, thr_batches[i], thr_sw_items[i], 1e6 * thr_sw_items[i] / sw_time[i]);
      total_item += thr_sw_items[i];
      total_throughput += (1e6 * thr_sw_items[i] / sw_time[i]);
      if (max_time < sw_time[i]) max_time = sw_time[i];
    }
    printf("\n\tTotal SWs: %lu, Max time = %0.4f, Throughput = %0.2f SW/s\n", total_item, max_time / 1e6, total_throughput);
  }
}
