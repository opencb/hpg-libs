#include "dna_aligner.h"
#include "rna_aligner.h"

double emboss_matrix_t = 0.0f, emboss_tracking_t = 0.0f;
double sse_matrix_t = 0.0f, sse_tracking_t = 0.0f;
double sse1_matrix_t = 0.0f, sse1_tracking_t = 0.0f;
double avx_matrix_t = 0.0f, avx_tracking_t = 0.0f;
double avx1_matrix_t = 0.0f, avx1_tracking_t = 0.0f;

//--------------------------------------------------------------------
// constants
//--------------------------------------------------------------------

#define OPTIONS 			30
#define MIN_ARGC  			5
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

void run_dna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     pair_mng_t *pair_mng, options_t *options);


void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg,
		     options_t *options);

//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------
int main(int argc, char* argv[]) {
  const char HEADER_FILE[1024] = "Human_NCBI37.hbam\0";

  //log_level = LOG_DEBUG_LEVEL;
  log_verbose = 1;

  if (argc <= 1) {
    usage_cli();
  }

  char *command = argv[1];  
  // We need to consume command: {dna | rna | bs | build-index}
  argc -= 1;
  argv += 1;

  if(strcmp(command, "dna") != 0 && 
     strcmp(command, "rna") != 0 &&
     strcmp(command, "bs" ) != 0 && 
     strcmp(command, "build-index") != 0) {
    LOG_FATAL("Command Mode Unknown.");
  }

  // parsing options
  options_t *options = parse_options(argc, argv);
  log_level = options->log_level;
  validate_options(options, command);
  LOG_DEBUG_F("Command Mode: %s\n", command);

  if (!strcmp(command, "build-index")) {
       run_index_builder(options->genome_filename, options->bwt_dirname, options->index_ratio);
       LOG_DEBUG("Done !!\n");
       exit(0);
  }


  time_on =  (unsigned int) options->timming;
  statistics_on =  (unsigned int) options->statistics;

  // genome parameters
  LOG_DEBUG("Reading genome...");
  genome_t* genome = genome_new("dna_compression.bin", options->bwt_dirname);
  LOG_DEBUG("Done !!");
  
  // BWT parameters
  LOG_DEBUG("Reading bwt index...");
  //if (time_on) { timing_start(INIT_BWT_INDEX, 0, timing_p); }
  bwt_index_t *bwt_index = bwt_index_new(options->bwt_dirname);
  LOG_DEBUG("Reading bwt index done !!");
  
  //BWT parameters
  bwt_optarg_t *bwt_optarg = bwt_optarg_new(1, options->bwt_threads, 100,
					    options->report_best,
					    options->report_n_hits, 
					    options->report_all);
  
  // CAL parameters
  cal_optarg_t *cal_optarg = cal_optarg_new(options->min_cal_size, options->seeds_max_distance, 
					    options->min_num_seeds, options->max_num_seeds,
					    options->seed_size, options->min_seed_size, 
					    options->cal_seeker_errors);
  
  // paired mode parameters
  pair_mng_t *pair_mng = pair_mng_new(options->pair_mode, options->pair_min_distance, 
				      options->pair_max_distance);
  
  LOG_DEBUG("init table...");
  initTable();
  LOG_DEBUG("init table done !!");
  
  if (!strcmp(command, "rna")) {
    run_rna_aligner(genome, bwt_index, pair_mng, bwt_optarg, cal_optarg, options);
  } else {
    // DNA version
    run_dna_aligner(genome, bwt_index, bwt_optarg, cal_optarg, pair_mng, options);
  }

  LOG_DEBUG("main done !!");

  // Free memory
  if (time_on) { 
    timing_start(FREE_MAIN, 0, timing_p); 
  }
  
  bwt_index_free(bwt_index);
  genome_free(genome);
  bwt_optarg_free(bwt_optarg);
  cal_optarg_free(cal_optarg);
  pair_mng_free(pair_mng);
  
  if (time_on) { timing_stop(FREE_MAIN, 0, timing_p); }

  if (time_on) { 
    timing_display(timing_p);
  }

  basic_statistics_display(basic_st, !strcmp(command, "rna"));

  if (time_on){ timing_free(timing_p); }

  options_free(options);

  return 0;
}

//--------------------------------------------------------------------
