#include "options.h"

const char DEFAULT_OUTPUT_NAME[30] = "hpg-aligner_output";
//const char SPLICE_EXACT_FILENAME[30]   = "exact_junctions.bed";
//const char SPLICE_EXTEND_FILENAME[30]  = "extend_junctions.bed";
//const char INDEX_NAME[30]  = "index";

options_t *options_new(void) {
  options_t *options = (options_t*) calloc (1, sizeof(options_t));
  size_t num_cores = 0;
  
  //======================= COMMON OPTIONS ====================
  options->in_filename = NULL;
  options->in_filename2 = NULL;
  options->report_all =  0;
  options->log_level = 0;
  options->output_name = strdup(DEFAULT_OUTPUT_NAME);
  options->num_gpu_threads = DEFAULT_GPU_THREADS;
  //GET Number System Cores
  //----------------------------------------------
  if (num_cores = get_optimal_cpu_num_threads()) {
    options->num_cpu_threads = num_cores;
  }else {
    options->num_cpu_threads = DEFAULT_CPU_THREADS;
  }
  //----------------------------------------------
  options->max_intron_length = DEFAULT_MAX_INTRON_LENGTH;
  options->min_num_seeds = DEFAULT_MIN_NUM_SEEDS;
  options->max_num_seeds = DEFAULT_MAX_NUM_SEEDS;
  options->cal_seeker_errors = DEFAULT_CAL_SEEKER_ERRORS;
  options->write_size = DEFAULT_WRITE_BATCH_SIZE;
  options->bwt_threads = DEFAULT_BWT_THREADS;
  options->region_threads = DEFAULT_REGION_THREADS;
  options->num_cal_seekers = DEFAULT_NUM_CAL_SEEKERS;
  options->num_sw_servers = DEFAULT_NUM_SW_THREADS;
  options->min_score = DEFAULT_SW_MIN_SCORE;
  options->match = DEFAULT_SW_MATCH;
  options->mismatch = DEFAULT_SW_MISMATCH;
  options->gap_open = DEFAULT_SW_GAP_OPEN;
  options->gap_extend = DEFAULT_SW_GAP_EXTEND;
  options->min_intron_length = DEFAULT_MIN_INTRON_LENGTH;
  options->pair_mode = DEFAULT_PAIR_MODE;
  options->pair_min_distance = DEFAULT_PAIR_MIN_DISTANCE;
  options->pair_max_distance = DEFAULT_PAIR_MAX_DISTANCE;
  options->timming = 0;
  options->statistics = 0;
  options->report_best = 0;
  options->report_n_hits = 0;
  options->gpu_process = 0;
  options->bwt_set = 0;
  options->reg_set = 0;
  options->cal_set = 0;
  options->sw_set = 0;
  //=========================================================
  options->min_cal_size = 0; 
  options->seeds_max_distance = 0;
  options->batch_size = 0;
  options->min_seed_size = 0;
  options->seed_size = 0;
  options->flank_length = 0;
  
  return options;
}

void validate_options(options_t *options, char *mode) {
  int value_dir = exists(options->output_name);
  int DEFAULT_READ_BATCH_SIZE;
  int DEFAULT_SEED_SIZE;
  int DEFAULT_FLANK_LENGTH;
  int DEFAULT_MIN_SEED_SIZE;
  int DEFAULT_MIN_CAL_SIZE;
  int DEFAULT_SEEDS_MAX_DISTANCE;

  if (strcmp("rna", mode) == 0) {
    DEFAULT_READ_BATCH_SIZE = 20000;
    DEFAULT_SEED_SIZE	= 20;
    DEFAULT_FLANK_LENGTH = 5;
    DEFAULT_MIN_SEED_SIZE = 16;
    DEFAULT_MIN_CAL_SIZE = 30;
    DEFAULT_SEEDS_MAX_DISTANCE = 100;
  }else if (strcmp("dna", mode) == 0) {
    DEFAULT_READ_BATCH_SIZE = 200000;
    DEFAULT_SEED_SIZE = 15;
    DEFAULT_FLANK_LENGTH = 30;
    DEFAULT_MIN_SEED_SIZE = 15;
    DEFAULT_MIN_CAL_SIZE = 20;
    DEFAULT_SEEDS_MAX_DISTANCE = 60;
  }

  if (!value_dir) {
    create_directory(options->output_name);
    exit(-1);
  }

  if (!options->min_cal_size) {
    options->min_cal_size = DEFAULT_MIN_CAL_SIZE;
  }
 
  if (!options->seeds_max_distance) {
    options->seeds_max_distance = DEFAULT_SEEDS_MAX_DISTANCE;
  }
  
  if (!options->batch_size) {
    options->batch_size = DEFAULT_READ_BATCH_SIZE;
  }
   
  if (!options->min_seed_size) {
    options->min_seed_size = DEFAULT_MIN_SEED_SIZE;
  }
  
  if (!options->seed_size) {
    options->seed_size = DEFAULT_SEED_SIZE;
  }

  if (!options->flank_length) {
    options->flank_length = DEFAULT_FLANK_LENGTH;
  }
}


void options_free(options_t *options) {
	if(options == NULL) { return; }
	//if (options->splice_exact_filename != NULL)	{ free(options->splice_exact_filename); }
	//if (options->splice_extend_filename  != NULL)	{ free(options->splice_extend_filename); }
	if (options->in_filename  != NULL)		{ free(options->in_filename); }
	if (options->in_filename2  != NULL)		{ free(options->in_filename2); }
	if (options->bwt_dirname  != NULL)		{ free(options->bwt_dirname); }
	if (options->output_name  != NULL)		{ free(options->output_name); }

	free(options);
}


void options_display(options_t *options) {

  char* in_filename = strdup(options->in_filename);
  char* in_filename2 = NULL;
  if (options->in_filename2 != NULL) {
    in_filename2 = strdup(options->in_filename2);
  }
  char* bwt_dirname =  strdup(options->bwt_dirname);
  //char* genome_filename =  strdup(options->genome_filename);
  unsigned int  report_all = (unsigned int)options->report_all;
  unsigned int  report_best = (unsigned int)options->report_best;
  unsigned int  report_n_hits = (unsigned int)options->report_n_hits;

  if ((report_best == 0) && (report_n_hits == 0)) {
    report_all = 1;
  }

  char* output_name =  strdup(options->output_name);
  unsigned int num_gpu_threads =  (unsigned int)options->num_gpu_threads;
  unsigned int num_cpu_threads =  (unsigned int)options->num_cpu_threads;
  unsigned int rna_seq =  (unsigned int)options->rna_seq; 
  unsigned int cal_seeker_errors =  (unsigned int)options->cal_seeker_errors; 
  unsigned int min_cal_size =  (unsigned int)options->min_cal_size; 
  unsigned int seeds_max_distance =  (unsigned int)options->seeds_max_distance; 
  unsigned int bwt_threads =  (unsigned int)options->bwt_threads; 
  unsigned int batch_size =  (unsigned int)options->batch_size; 
  unsigned int write_size =  (unsigned int)options->write_size;  
  unsigned int num_cal_seekers =  (unsigned int)options->num_cal_seekers;
  unsigned int region_threads =  (unsigned int)options->region_threads;
  unsigned int num_sw_servers =  (unsigned int)options->num_sw_servers;
  unsigned int min_seed_size =  (unsigned int)options->min_seed_size;
  unsigned int seed_size =  (unsigned int)options->seed_size;
  unsigned int min_num_seeds =  (unsigned int)options->min_num_seeds;
  unsigned int max_num_seeds =  (unsigned int)options->max_num_seeds;
  unsigned int max_intron_length =  (unsigned int)options->max_intron_length;
  unsigned int flank_length =  (unsigned int)options->flank_length;
  unsigned int pair_mode =  (unsigned int)options->pair_mode;
  unsigned int pair_min_distance =  (unsigned int)options->pair_min_distance;
  unsigned int pair_max_distance =  (unsigned int)options->pair_max_distance;
  unsigned int min_intron_length =  (unsigned int)options->min_intron_length;
  unsigned int gpu_process = (unsigned int)options->gpu_process;
  float min_score =  (float)options->min_score;
  float match =   (float)options->match;
  float mismatch =   (float)options->mismatch;
  float gap_open =   (float)options->gap_open;
  float gap_extend =   (float)options->gap_extend;
  
  printf("PARAMETERS CONFIGURATION\n");
  printf("=================================================\n");
  printf("Num gpu threads %d\n", num_gpu_threads);
  printf("GPU Process: %s\n",  gpu_process == 0 ? "Disable":"Enable");
  printf("Num cpu threads %d\n",  num_cpu_threads);
  printf("RNA Server: %s\n",  rna_seq == 0 ? "Disable":"Enable");
  printf("Report all hits: %s\n",  report_all == 0 ? "Disable":"Enable");
  printf("Report best hits: %d\n",  report_best);
  printf("Report n hits: %d\n",  report_n_hits);
  printf("CAL seeker errors: %d\n",  cal_seeker_errors);
  printf("Batch size: %dBytes\n",  batch_size);
  printf("Write size: %dBytes\n",  write_size);
  printf("BWT Threads: %d\n",  bwt_threads);
  printf("Region Threads: %d\n",  region_threads);
  printf("Num CAL seekers: %d\n", num_cal_seekers);
  printf("Num SW servers: %d\n",  num_sw_servers);
  printf("SEEDING and CAL PARAMETERS\n");
  printf("\tMin. number of seeds: %d\n",  min_num_seeds);
  printf("\tMax. number of seeds: %d\n",  max_num_seeds);
  printf("\tSeed size: %d\n",  seed_size);
  printf("\tMin seed size: %d\n",  min_seed_size);
  printf("\tMin CAL size: %d\n",  min_cal_size);
  printf("\tSeeds max distance: %d\n",  seeds_max_distance);
  printf("\tFlank length: %d\n", flank_length);
  printf("RNA PARAMETERS\n");
  printf("\tMax intron length: %d\n", max_intron_length);
  printf("\tMin intron length: %d\n", min_intron_length);
  printf("PAIR-MODE PARAMETERS\n");
  printf("\tPair mode: %d\n", pair_mode);
  printf("\tMin. distance: %d\n", pair_min_distance);
  printf("\tMax. distance: %d\n", pair_max_distance);
  printf("SMITH-WATERMAN PARAMETERS\n");
  printf("\tMin score  : %0.4f\n", min_score);
  printf("\tMatch      : %0.4f\n", match);
  printf("\tMismatch   : %0.4f\n", mismatch);
  printf("\tGap open   : %0.4f\n", gap_open);
  printf("\tGap extend : %0.4f\n", gap_extend);
  printf("=================================================\n");

  free(in_filename);
  if (in_filename2 != NULL) free(in_filename2);
  free(bwt_dirname);
  //free(genome_filename);
  free(output_name);
}

//--------------------------------------------------------------------

void** argtable_options_new(void) {
	void **argtable = (void**)malloc((NUM_OPTIONS + 1) * sizeof(void*));	// NUM_OPTIONS +1 to allocate end structure
	// NOTICE that order cannot be changed as is accessed by index in other functions
	argtable[0] = arg_file0("f", "fq,fastq", NULL, "Reads file input");
	argtable[1] = arg_file0("i", "bwt-index", NULL, "BWT directory name");
	argtable[2] = arg_int0("l", "log-level", NULL, "Log debug level");
	argtable[3] = arg_lit0(NULL, "report-all", "Report all alignments");
	argtable[4] = arg_file0("o", "outdir", NULL, "Output directory");
	argtable[5] = arg_int0(NULL, "gpu-threads", NULL, "Number of GPU Threads");
	argtable[6] = arg_int0(NULL, "cpu-threads", NULL, "Number of CPU Threads");
	argtable[7] = arg_int0(NULL, "xxxx", NULL, "xxxxx");
	argtable[8] = arg_int0(NULL, "cal-seeker-errors", NULL, "Number of errors in CAL Seeker");
	argtable[9] = arg_int0(NULL, "min-cal-size", NULL, "Minimum CAL size");
	argtable[10] = arg_int0(NULL, "max-distance-seeds", NULL, "Maximum distance between seeds");
	argtable[11] = arg_int0(NULL, "read-batch-size", NULL, "Batch Size");
	argtable[12] = arg_int0(NULL, "write-batch-size", NULL, "Write Size");
	argtable[13] = arg_int0(NULL, "num-cal-seekers", NULL, "Number of CAL Seekers");
	argtable[14] = arg_int0(NULL, "num-sw-servers", NULL, "Number of Smith-Waterman servers");
	argtable[15] = arg_int0(NULL, "num-bwt-threads", NULL, "Number of BWT threads");
	argtable[16] = arg_int0(NULL, "num-region-threads", NULL, "Number of region threads");
	argtable[17] = arg_int0(NULL, "seed-size", NULL, "Number of nucleotides in a seed");
	argtable[18] = arg_int0(NULL, "min-seed-size", NULL, "Minimum number of nucleotides in a seed");
	argtable[19] = arg_int0(NULL, "cal-flank-size", NULL, "Flank length for CALs");
	argtable[20] = arg_dbl0(NULL, "sw-match", NULL, "Match value for Smith-Waterman algorithm");
	argtable[21] = arg_dbl0(NULL, "sw-mismatch", NULL, "Mismatch value for Smith-Waterman algorithm");
	argtable[22] = arg_dbl0(NULL, "sw-gap-open", NULL, "Gap open penalty for Smith-Waterman algorithm");
	argtable[23] = arg_dbl0(NULL, "sw-gap-extend", NULL, "Gap extend penalty for Smith-Waterman algorithm");
	argtable[24] = arg_dbl0(NULL, "sw-min-score", NULL, "Minimum score for valid mappings");
	argtable[25] = arg_int0(NULL, "max-intron-size", NULL, "Maximum intron size");
	argtable[26] = arg_int0(NULL, "min-intron-size", NULL, "Minimum intron size");
	argtable[27] = arg_lit0("t", "time", "Timming mode active");
	argtable[28] = arg_lit0("s", "stats", "Statistics mode active");
	argtable[29] = arg_lit0("h", "help", "Help option");
	argtable[30] = arg_file0(NULL, "xxxxx", NULL, "xxxxxxxx");
	argtable[31] = arg_file0(NULL, "xxxxx", NULL, "xxxxxxxx");
	argtable[32] = arg_file0("j", "fq2,fastq2", NULL, "Reads file input #2 (for paired mode)");
	argtable[33] = arg_int0(NULL, "paired-mode", NULL, "Pair mode: 0 = single-end, 1 = paired-end, 2 = mate-pair [Default 0]");
	argtable[34] = arg_int0(NULL, "paired-min-distance", NULL, "Minimum distance between pairs");
	argtable[35] = arg_int0(NULL, "paired-max-distance", NULL, "Maximum distance between pairs");
	argtable[36] = arg_int0(NULL, "report-n-best", NULL, "Report the <n> best alignments");
	argtable[37] = arg_int0(NULL, "report-n-hits", NULL, "Report <n> hits");
	argtable[38] = arg_int0(NULL, "min-num-seeds", NULL, "Minimum number of seeds per read");
	argtable[39] = arg_int0(NULL, "max-num-seeds", NULL, "Maximum number of seeds per read");
	argtable[40] = arg_lit0(NULL, "gpu-enable", "Enable GPU Process");
	argtable[41] = arg_end(20);

	return argtable;
}


void argtable_options_free(void **argtable) {
	if(argtable != NULL) {
		arg_freetable(argtable, NUM_OPTIONS + 1);	// struct end must also be freed
		free(argtable);
	}
}



int read_config_file(const char *filename, options_t *options) {
	if (filename == NULL || options == NULL) {
		return -1;
	}

	config_t *config = (config_t*) calloc (1, sizeof(config_t));
	int ret_code = config_read_file(config, filename);
	if (ret_code == CONFIG_FALSE) {
		LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
		return -1;
	}

	const char *tmp_string;
	long tmp_int;

	/*if(config_lookup_string(config, "app.outdir", &tmp_string)) { options->output_directory = strdup(tmp_string); }
	if(config_lookup_int(config, "app.cpu-num-threads", &tmp_int)) { options->cpu_num_threads = (int)tmp_int; }
	*/

	config_destroy(config);
	free(config);
//	free(tmp_string);

	return ret_code;
}


/**
 * @brief Initializes an options_t structure from argtable parsed CLI with default values. Notice that options are order dependent.
 * @return A new options_t structure initialized with default values.
 *
 * Initializes the only default options from options_t.
 */
options_t *read_CLI_options(void **argtable, options_t *options) {	
  if (((struct arg_file*)argtable[0])->count) { options->in_filename = strdup(*(((struct arg_file*)argtable[0])->filename)); }
  if (((struct arg_file*)argtable[1])->count) { options->bwt_dirname = strdup(*(((struct arg_file*)argtable[1])->filename)); }
  if (((struct arg_file*)argtable[2])->count) { options->log_level = *(((struct arg_int*)argtable[2])->ival); }
  if (((struct arg_file*)argtable[3])->count) { options->report_all = (((struct arg_int *)argtable[3])->count); }
  if (((struct arg_file*)argtable[4])->count) { free(options->output_name); options->output_name = strdup(*(((struct arg_file*)argtable[4])->filename)); }  
  if (((struct arg_int*)argtable[5])->count) { options->num_gpu_threads = *(((struct arg_int*)argtable[5])->ival); }
  if (((struct arg_int*)argtable[6])->count) { options->num_cpu_threads = *(((struct arg_int*)argtable[6])->ival); }
  if (((struct arg_int*)argtable[8])->count) { options->cal_seeker_errors = *(((struct arg_int*)argtable[8])->ival); }
  if (((struct arg_int*)argtable[9])->count) { options->min_cal_size = *(((struct arg_int*)argtable[9])->ival); }
  if (((struct arg_int*)argtable[10])->count) { options->seeds_max_distance = *(((struct arg_int*)argtable[10])->ival); }
  if (((struct arg_int*)argtable[11])->count) { options->batch_size = *(((struct arg_int*)argtable[11])->ival); }
  if (((struct arg_int*)argtable[12])->count) { options->write_size = *(((struct arg_int*)argtable[12])->ival); }
  if (((struct arg_int*)argtable[13])->count) { options->cal_set = 1; options->num_cal_seekers = *(((struct arg_int*)argtable[13])->ival); }
  if (((struct arg_int*)argtable[14])->count) { options->sw_set = 1; options->num_sw_servers = *(((struct arg_int*)argtable[14])->ival); }
  if (((struct arg_int*)argtable[15])->count) { options->bwt_set = 1; options->bwt_threads = *(((struct arg_int*)argtable[15])->ival); }
  if (((struct arg_int*)argtable[16])->count) { options->reg_set = 1; options->region_threads = *(((struct arg_int*)argtable[16])->ival); }
  if (((struct arg_int*)argtable[17])->count) { options->seed_size = *(((struct arg_int*)argtable[17])->ival); }
  if (((struct arg_int*)argtable[18])->count) { options->min_seed_size = *(((struct arg_int*)argtable[18])->ival); }
  if (((struct arg_int*)argtable[19])->count) { options->flank_length = *((struct arg_int*)argtable[19])->ival; }
  if (((struct arg_dbl*)argtable[20])->count) { options->match = *((struct arg_dbl*)argtable[20])->dval; }
  if (((struct arg_dbl*)argtable[21])->count) { options->mismatch = *(((struct arg_dbl*)argtable[21])->dval); }
  if (((struct arg_dbl*)argtable[22])->count) { options->gap_open = *(((struct arg_dbl*)argtable[22])->dval); }
  if (((struct arg_dbl*)argtable[23])->count) { options->gap_extend = *(((struct arg_dbl*)argtable[23])->dval); }
  if (((struct arg_dbl*)argtable[24])->count) { options->min_score = *(((struct arg_dbl*)argtable[24])->dval); }
  if (((struct arg_int*)argtable[25])->count) { options->max_intron_length = *(((struct arg_int*)argtable[25])->ival); }
  if (((struct arg_int*)argtable[26])->count) { options->min_intron_length = *(((struct arg_int*)argtable[26])->ival); }
  if (((struct arg_int*)argtable[27])->count) { options->timming = ((struct arg_int*)argtable[27])->count; }
  if (((struct arg_int*)argtable[28])->count) { options->statistics = ((struct arg_int*)argtable[28])->count; }
  if (((struct arg_int*)argtable[29])->count) { options->help = ((struct arg_int*)argtable[29])->count; }
  if (((struct arg_file*)argtable[32])->count) { options->in_filename2 = strdup(*(((struct arg_file*)argtable[32])->filename)); }
  if (((struct arg_int*)argtable[33])->count) { options->pair_mode = *(((struct arg_int*)argtable[33])->ival); }
  if (((struct arg_int*)argtable[34])->count) { options->pair_min_distance = *(((struct arg_int*)argtable[34])->ival); }
  if (((struct arg_int*)argtable[35])->count) { options->pair_max_distance = *(((struct arg_int*)argtable[35])->ival); }
  if (((struct arg_int*)argtable[36])->count) { options->report_best = *(((struct arg_int*)argtable[36])->ival); }
  if (((struct arg_int*)argtable[37])->count) { options->report_n_hits = *(((struct arg_int*)argtable[37])->ival); }
  if (((struct arg_int*)argtable[38])->count) { options->min_num_seeds = *(((struct arg_int*)argtable[38])->ival); }
  if (((struct arg_int*)argtable[39])->count) { options->max_num_seeds = *(((struct arg_int*)argtable[39])->ival); }
  if (((struct arg_int*)argtable[40])->count) { 
    #ifdef HPG_GPU
       options->gpu_process = (((struct arg_int *)argtable[40])->count); 
    #else
       options->gpu_process = 0; 
    #endif
  }

  return options;
}


options_t *parse_options(int argc, char **argv) {
  void **argtable = argtable_options_new();
  //	struct arg_end *end = arg_end(10);
  //	void **argtable = argtable_options_get(argtable_options, end);
  
  options_t *options = options_new();
  if (argc < 2) {
    usage(argtable);
    exit(-1);
  }else {

    int num_errors = arg_parse(argc, argv, argtable);

    if (((struct arg_int*)argtable[29])->count) {
      usage(argtable);
	argtable_options_free(argtable);
	options_free(options);
	exit(0);
    }
        
    if (num_errors > 0) {
      arg_print_errors(stdout, argtable[NUM_OPTIONS], "hpg-aligner");	// struct end is always allocated in the last position
      usage(argtable);
      exit(-1);
    }else {
      options = read_CLI_options(argtable, options);
      if(options->help) {
	usage(argtable);
	argtable_options_free(argtable);
	options_free(options);
	exit(0);
      }
      // Check if 'help' option has been provided.
    }
    
  }
  //	exit:
  argtable_options_free(argtable);
  //	free(end);
  //	free(argtable_options);

  // in previous versions, min. flank length was 20
  //  if (options->flank_length < 5) {
  //    options->flank_length = 5;
  //  }

  return options;
}

void usage(void **argtable) {
  printf("Usage:\n./hpg-aligner {dna | rna | bs | build-index}");
  arg_print_syntaxv(stdout, argtable, "\n");
  arg_print_glossary(stdout, argtable, "%-50s\t%s\n");
}

void usage_cli() {
  void **argtable = argtable_options_new();
  usage(argtable);
  exit(0);
}
