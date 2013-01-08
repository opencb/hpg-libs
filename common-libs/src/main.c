#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <omp.h>

#include "commons/log.h"

#include "bioformats/fastq/fastq_batch_reader.h"

#include "error.h"

#include "timing.h"
#include "buffers.h"
#include "bwt_server.h"
#include "batch_writer.h"
#include "region_seeker.h"
#include "cal_seeker.h"
#include "sw_server.h"
#include "pair_server.h"
#include "batch_aligner.h"

#ifdef HPG_GPU
#include "aligners/bwt/gpu.h"
#include "aligners/bwt/bwt_gpu.h"
#endif

// rna server
#include "rna_splice.h"
#include "rna_server.h"


//#include "sw.h"
//#include "smith_waterman.h"

#include "options.h"
#include "statistics.h"



double emboss_matrix_t = 0.0f, emboss_tracking_t = 0.0f;
double sse_matrix_t = 0.0f, sse_tracking_t = 0.0f;
double sse1_matrix_t = 0.0f, sse1_tracking_t = 0.0f;
double avx_matrix_t = 0.0f, avx_tracking_t = 0.0f;
double avx1_matrix_t = 0.0f, avx1_tracking_t = 0.0f;

//char *qq3, *qq2, *qq1, *qq0;
//char *rr3, *rr2, *rr1, *rr0;

//--------------------------------------------------------------------
// constants
//--------------------------------------------------------------------

#define OPTIONS 			30
#define MIN_ARGC  			5
#define NUM_SECTIONS_TIME 		7
#define NUM_SECTIONS_STATISTICS 	5
#define NUM_SECTIONS_STATISTICS_SB	21
//#define REQUIRED 1
//#define NO_REQUIRED 0

//--------------------------------------------------------------------
// global variables for log functions
//--------------------------------------------------------------------

//int log_level = DEFAULT_LOG_LEVEL;
//bool verbose = true;

//--------------------------------------------------------------------
// global variables for timing and capacity meausures
//--------------------------------------------------------------------

char time_on = 0;
char statistics_on = 0;
timing_t* timing_p = NULL;
int num_of_chromosomes;
statistics_t *statistics_p;
basic_statistics_t basic_st;

double kl_time;

#ifdef HPG_GPU
void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     gpu_context_t *context, options_t *options);

void run_dna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     pair_mng_t *pair_mng, gpu_context_t *context, options_t *options);

#else
   void run_dna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		        bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		        pair_mng_t *pair_mng, options_t *options);


void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg,
		     options_t *options);

#endif

//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------

int main(int argc, char* argv[]) {

  //  qq3 = (char *) malloc(18192); qq2 = (char *) malloc(18192); qq1 = (char *) malloc(18192); qq0 = (char *) malloc(18192);
  //  rr3 = (char *) malloc(18192); rr2 = (char *) malloc(18192); rr1 = (char *) malloc(18192); rr0 = (char *) malloc(18192);
	log_level = LOG_DEBUG_LEVEL;
	log_verbose = 1;

	char *command = argv[1];

		// We need to consume command: {qc | filter | prepro}
	argc -= 1;
	argv += 1;
	LOG_DEBUG_F("Command Mode: %s\n", command);
	if(strcmp()...) {

	}
  // parsing options
  options_t *options = parse_options(argc, argv);

  LOG_DEBUG_F("Command Mode: %s\n", command);
  //************** Set Threads to sections **************//
  size_t cpu_threads = options->num_cpu_threads; 
//  if (options->rna_seq) {
  if (!strcmp(command, "rna")) {
	LOG_DEBUG_F("Command Mode: %s\n", command);
    if (!options->bwt_set && 
	!options->reg_set && 
	!options->cal_set &&
	!options->sw_set &&
	cpu_threads > 4) {
      printf("Auto Thread configuration ...\n");
      if (cpu_threads == 5) { options->region_threads++; }
      else if (cpu_threads == 6) { 
	options->region_threads++; 
	options->num_sw_servers++; 
      }
      else {
	options->region_threads = options->num_cpu_threads / 2;
	cpu_threads -= options->region_threads;
	cpu_threads -= options->bwt_threads;
	
	options->num_sw_servers = (cpu_threads / 2) + 1;
	cpu_threads -= options->num_sw_servers;
	
	options->num_cal_seekers = cpu_threads;
      }
      printf("Set %d Threads successful\n", options->num_cpu_threads);
    }
  }
  //****************************************************//

  printf("done !\n");
  printf("displaying options...\n");
  // display selected options
  options_display(options);

  time_on =  (unsigned int) options->timming;
  statistics_on =  (unsigned int) options->statistics;


  // timing
  if (time_on) { 
    char* labels_time[NUM_SECTIONS_TIME] = {"Index Initialization       ", 
					    "BWT server                 ", 
					    "Region Seeker              ", 
					    "CAL seeker                 ", 
					    "Rna server                 ", 
					    "Free main memory           ", 
					    "Total time                 "};
    
    int num_threads[NUM_SECTIONS_TIME] = {1, 1, 1, options->num_cal_seekers, 
					  options->num_sw_servers, 1, 1};
    timing_p = timing_new((char**) labels_time, (int*) num_threads, NUM_SECTIONS_TIME);
    
    //timing_start(INIT_INDEX, 0, timing_p);    
  }

  //if (statistics_on) {
    /*
     * FastqQ reader      : Num Batches %d 
     * BWT server         : Num Reads process %d, Num reads mapped %d, Num reads unmapped %d, Num mappings report %d 
     * Region Seeker      : Total num regions %d
     * CAL seeker         : Total cals %d, Num Reads unmapped %d
     * Rna server         : Num reads mapped %d, Num reads unmapped %d 
     * Batch writer       : Total num mappings report %d
     * Total 		  : Total reads %d, Total reads mapped %d, Total reads unmapped %d, Total splice junctions %d 	
     */				
    /*    char* labels_statistics[NUM_SECTIONS_STATISTICS] = { "BWT server                 ", 
							 "Region seeker              ", 
							 "CAL seeker                 ", 
							 "Rna server                 ",
							 "Total 	Statistics	    "};
    
    char *sub_labels_statistics[NUM_SECTIONS_STATISTICS_SB] = {	"Num batches process     ",
								"Reads process           ",
								"   Reads mapped         ", 
								"   Reads unmapped       ",
								"Mappings report         ",
								
								"Num batches process     ",
								"Num reads process       ",
								
								"Num batches process     ",
								"Total reads process     ",
								"Num reads unmapped      ",
								
								"Num batches process     ",
								"Num reads               ",
								"Num reads mapped        ",
								"Num reads unmapped      ",
								"Num sw process          ",
								"   Num sw valids        ",
								"   Num sw no valids     ",
								
								"Total reads             ",
								"   Total reads mapped   ",
								"   Total reads unmapped ",
								"Total splice junctions  "
								};
				      
    unsigned int num_values[NUM_SECTIONS_STATISTICS] = {5, 2, 3, 7, 4};
    statistics_p = statistics_new((char **)labels_statistics, (char **)sub_labels_statistics, 
				  (unsigned int *)num_values, NUM_SECTIONS_STATISTICS, NUM_SECTIONS_STATISTICS_SB);
    */
  //}

  // initialize some structures: Burrow-Wheeler objects, genome, nucletotide table...
  // (all these initializations could be performed in parallel)

  // genome parameters
  printf("Reading genome...\n");
  genome_t* genome = genome_new(options->genome_filename, options->bwt_dirname);
  printf("Done !!\n");

  /*
  {
    int num_targets = genome->num_chromosomes;

    char *filename = "./header.bam";
    bam_header_t *bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));
    bam_header->n_targets = num_targets;
    bam_header->target_name = (char **) calloc(num_targets, sizeof(char *));
    bam_header->target_len = (uint32_t*) calloc(num_targets, sizeof(uint32_t));
    for (int i = 0; i < num_targets; i++) {
      bam_header->target_name[i] = strdup(genome->chr_name[i]);
      bam_header->target_len[i] = genome->chr_size[i];
      printf("target %i: %s, len = %i\n", i, bam_header->target_name[i], bam_header->target_len[i]);
    }
    bam_header->text = strdup("@PG\tID:hpg-aligner\tVN:1.0\n");
    bam_header->l_text = strlen(bam_header->text);

    bam_file_t* bamf = bam_fopen_mode(filename, bam_header, "w");
    bam_fwrite_header(bam_header, bamf);
    bam_fclose(bamf);
    exit(-1);
  }
  */

  // BWT parameters
  printf("Reading bwt index...\n");
  if (time_on) { timing_start(INIT_BWT_INDEX, 0, timing_p); }
  bwt_index_t *bwt_index = bwt_index_new(options->bwt_dirname);
  printf("Reading bwt index done !!\n");

  // bwt_optarg_new(errors, threads, max aligns) 
  bwt_optarg_t *bwt_optarg = bwt_optarg_new(1, options->bwt_threads, 500,
					    options->report_best,
					    options->report_n_hits, 
					    options->report_all);
  //bwt_optarg_t *bwt_optarg = bwt_optarg_new(1, options->bwt_threads, 200);
  // CAL parameters
  //GOOD LUCK(20, 60, 18, 16, 0)
  cal_optarg_t *cal_optarg = cal_optarg_new(options->min_cal_size, options->seeds_max_distance, 
					    options->min_num_seeds, options->max_num_seeds,
					    options->seed_size, options->min_seed_size, 
					    options->cal_seeker_errors);
  #ifdef HPG_GPU
  gpu_context_t *context;
  if (options->gpu_process) {
    context = gpu_context_new(0, options->num_gpu_threads, bwt_index);      
  }
  #endif

  // paired mode parameters
  pair_mng_t *pair_mng = pair_mng_new(options->pair_mode, options->pair_min_distance, 
				      options->pair_max_distance);
    
  /*
  {
    //@6_13362459_13363092_0_1_0_0_3:0:0_2:0:0_13195/1
    //@3_17488330_17489008_0_1_0_0_3:0:0_1:0:0_7101/1
    char seq[1000];
    int chr = 3;
    int strand = 0;
    size_t start = 17488330;
    size_t end = 17488429;
    genome_read_sequence_by_chr_index(seq, strand, chr - 1, &start, &end, genome_p);
    printf("%s\n", seq);
    exit(-1);
  }
  */

  printf("init table...\n");
  initTable();
  printf("init table done !!\n");


 /* if (cuda) {
    int num_gpus = 2;
    context_p = (void *) gpu_context_new(num_gpus, num_gpu_threads);
  }*/
  
  if (time_on) { 
    timing_start(MAIN_INDEX, 0, timing_p);
  }

  // launch threads in parallel
  omp_set_nested(1);

//  if (options->rna_seq) {
  if (!strcmp(command, "rna")) {
    // RNA version
    #ifdef HPG_GPU
    run_rna_aligner(genome, bwt_index, pair_mng, bwt_optarg, cal_optarg, context, options);
    #else
    run_rna_aligner(genome, bwt_index, pair_mng, bwt_optarg, cal_optarg, options);
    #endif
  } else {
    // DNA version
    #ifdef HPG_GPU
      run_dna_aligner(genome, bwt_index, bwt_optarg, cal_optarg, pair_mng, context, options);
    #else
      run_dna_aligner(genome, bwt_index, bwt_optarg, cal_optarg, pair_mng, options);
    #endif
  }

  printf("\nmain done !!\n");
  if (time_on) { 
    timing_stop(MAIN_INDEX, 0, timing_p);
  }

  // free memory
  if (time_on) { 
    timing_start(FREE_MAIN, 0, timing_p); 
  }
  
  bwt_index_free(bwt_index);
  genome_free(genome);
  bwt_optarg_free(bwt_optarg);
  cal_optarg_free(cal_optarg);
  pair_mng_free(pair_mng);
  
  if (time_on) { timing_stop(FREE_MAIN, 0, timing_p); }
  
  /*if (cuda) {
    gpu_context_free((gpu_context_t*) context_p);
  }
  if (time_on) { timing_stop(FREE_INDEX, 0, timing_p); }*/

  if (time_on) { 
    timing_display(timing_p);
  }
  
  //if (statistics_on) { 
    //statistics_display(statistics_p); 
//  basic_statistics_display(basic_st, options->rna_seq);
  basic_statistics_display(basic_st, !strcmp(command, "rna"));
    //}
  
  //  if (statistics_on && time_on) { timing_and_statistics_display(statistics_p, timing_p); }

  if (time_on){ timing_free(timing_p); }

  //if (statistics_on) { statistics_free(statistics_p); }

  options_free(options);
  
  //  printf("CPU :: Time BWT Search %.3fs\n", time_bwt_seed/1000000);
  //  printf("CPU :: Time S Search %.3fs\n", time_search_seed/1000000);
  //  printf("GPU :: Time BWT Search %.3fs\n", kl_time/1000000);

  return 0;
}

//--------------------------------------------------------------------
//extern int mapped_by_bwt[100];

//extern int unmapped_by_max_cals_counter[100];
//extern int unmapped_by_zero_cals_counter[100];
//extern int unmapped_by_score_counter[100];

#ifdef HPG_GPU
   void run_dna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		        bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		        pair_mng_t *pair_mng, gpu_context_t *context, options_t *options) {
#else
   void run_dna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		        bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		        pair_mng_t *pair_mng, options_t *options) {
#endif

     /*
  // for debugging
  for (int i = 0; i < options->num_cpu_threads; i++) {
    mapped_by_bwt[i] = 0;
    unmapped_by_max_cals_counter[i] = 0;
    unmapped_by_zero_cals_counter[i] = 0;
    unmapped_by_score_counter[i] = 0;
  }
     */
  // lists to communicate/synchronize the different threads
  list_t read_list, write_list;
  list_init("read", 1, 24, &read_list);
  list_init("write", options->num_cpu_threads, 24, &write_list);
  
  printf("writers to read_list = %d\n", list_get_writers(&read_list));
  printf("writers to write_list = %d\n", list_get_writers(&write_list));

  
  bwt_server_input_t bwt_input;
  bwt_server_input_init(NULL, 0, bwt_optarg, bwt_index, 
			NULL, 0, NULL, &bwt_input);
  
  region_seeker_input_t region_input;
  #ifdef HPG_GPU
     region_seeker_input_init(NULL, cal_optarg, bwt_optarg, 
			      bwt_index, NULL, 0, options->gpu_process, 
			      context, &region_input);
  #else 
     region_seeker_input_init(NULL, cal_optarg, bwt_optarg, 
			      bwt_index, NULL, 0, options->gpu_process, 
			      &region_input);
  #endif

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
      //fastq_batch_reader(&input);
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
    
    
  //if (cuda) {
  //   gpu_context_free((gpu_context_t*) context_p);
  //   }
  //   if (time_on) { timing_stop(FREE_INDEX, 0, timing_p); }

    /*    
    int reads_by_bwt = 0, by_max_cals = 0, by_zero_cals, by_score = 0;
    for (int i = 0; i < options->num_cpu_threads; i++) {
      by_max_cals += unmapped_by_max_cals_counter[i];
      by_zero_cals += unmapped_by_zero_cals_counter[i];
      by_score += unmapped_by_score_counter[i];
      reads_by_bwt += mapped_by_bwt[i];
    }
    printf("mapped by BWT = %d\n", reads_by_bwt);
    
    printf("unmapped by MAX_CALS = %d\n", by_max_cals);
    printf("unmapped by ZERO_CALS = %d\n", by_zero_cals);
    //  printf("unmapped by SW score = %d\n", by_score);
    */
  }
}

//--------------------------------------------------------------------
#ifdef HPG_GPU
   void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng, 
			bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
			gpu_context_t *context, options_t *options) {
#else
     void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng,
			  bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
			  options_t *options) {
#endif

  list_t read_list;
  list_init("read", 1, 10, &read_list);
  
  list_t unmapped_reads_list;
  list_init("unmmaped reads", 1, 10, &unmapped_reads_list);
  
  list_t regions_list;
  list_init("regions", 1, 10, &regions_list);
  
  list_t sw_list;
  list_init("sw", options->num_cal_seekers, 10, &sw_list);

  list_t alignment_list;
  list_init("alignments", options->num_sw_servers, 10, &alignment_list);

  list_t write_list;
  list_init("write", 1, 10, &write_list);

  allocate_splice_elements_t chromosome_avls[CHROMOSOME_NUMBER];
  init_allocate_splice_elements(chromosome_avls);

  #pragma omp parallel sections num_threads(7)
  {
    printf("Principal Sections %d threads\n", omp_get_num_threads());
    #pragma omp section
    {
      fastq_batch_reader_input_t input;
      fastq_batch_reader_input_init(options->in_filename, options->in_filename2, 
				    options->pair_mode, options->batch_size, 
				    &read_list, &input);
      //fastq_batch_reader(&input);
      fastq_batch_reader_aligner(&input);
      //fastq_batch_reader_array_list_single(&input);
    }
    #pragma omp section
    {
	bwt_server_input_t input;
        bwt_server_input_init(&read_list,  options->batch_size,  bwt_optarg, 
			      bwt_index, &write_list,  options->write_size, 
			      &unmapped_reads_list, &input);
	bwt_server_cpu(&input, pair_mng);
    }
    #pragma omp section
    {
      
      region_seeker_input_t input;
      #ifdef HPG_GPU
         region_seeker_input_init(&unmapped_reads_list, cal_optarg, 
			          bwt_optarg, bwt_index, &regions_list, 
			          options->region_threads, options->gpu_process, context, &input);
      #else 
         region_seeker_input_init(&unmapped_reads_list, cal_optarg, 
			          bwt_optarg, bwt_index, &regions_list, 
			          options->region_threads, options->gpu_process, &input);
      #endif

      region_seeker_server(&input);
    }
    #pragma omp section
    {
      cal_seeker_input_t input;
      cal_seeker_input_init(&regions_list, cal_optarg, &write_list, options->write_size, &sw_list, NULL, &input);
      #pragma omp parallel num_threads(options->num_cal_seekers)
      {
	cal_seeker_server(&input);
      }
      
    }
    #pragma omp section
    {
      sw_server_input_t input;
      sw_server_input_init(&sw_list, &alignment_list, options->write_size,  options->match,  
			   options->mismatch,  options->gap_open, options->gap_extend,  
			   options->min_score,  options->flank_length, genome,  
			   options->max_intron_length, options->min_intron_length,  
			   options->seeds_max_distance,  bwt_optarg, &input);
      //list_incr_writers(&write_list);
      #pragma omp parallel num_threads( options->num_sw_servers)
      {
	rna_server_omp_smith_waterman(&input, chromosome_avls);
      }

      write_chromosome_avls(chromosome_avls, &write_list,  
			    options->splice_extend_filename, 
			    options->splice_exact_filename,
			    options->write_size);
    }
     #pragma omp section
    {
	 pair_server_input_t input;
	 pair_server_input_init(pair_mng, 0, 0, 0, &alignment_list, NULL, &write_list, &input);
	 prepare_pair_server(&input);
     }
     
     #pragma omp section
    {
      batch_writer_input_t input;
      batch_writer_input_init( options->output_filename,  
			       options->splice_exact_filename,  
			       options->splice_extend_filename,  
			       &write_list, genome, &input);
      batch_writer2(&input);
      }
  }  
  
  
     }

//--------------------------------------------------------------------
//--------------------------------------------------------------------
