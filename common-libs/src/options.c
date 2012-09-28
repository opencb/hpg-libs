#include "options.h"

//============================ DEFAULT VALUES ============================

#define DEFAULT_GPU_THREADS		32
#define DEFAULT_CPU_THREADS		6
#define DEFAULT_CAL_SEEKER_ERROS	0
#define DEFAULT_MIN_CAL_SIZE		20
#define DEFAULT_SEEDS_MAX_DISTANCE	60
#define DEFAULT_BWT_THREADS		1
#define DEFAULT_BATCH_SIZE		200000
#define DEFAULT_WRITE_SIZE		500000
#define DEFAULT_NUM_CAL_SEEKERS		1
#define DEFAULT_REGION_THREADS		1
#define DEFAULT_NUM_SW_THREADS		1
#define DEFAULT_MIN_SEED_SIZE		16
#define DEFAULT_SEED_SIZE		18
#define DEFAULT_MAX_INTRON_LENGTH	1000000
#define DEFAULT_FLANK_LENGTH		20
#define DEFAULT_SW_MIN_SCORE		300
#define DEFAULT_SW_MATCH		5
#define DEFAULT_SW_MISMATCH		-4
#define DEFAULT_SW_GAP_OPEN		10
#define DEFAULT_SW_GAP_EXTEND		0.5
#define DEFAULT_MIN_INTRON_LENGTH	40


const char DEFAULT_OUTPUT_FILENAME[30] = "reads_results.bam";
const char SPLICE_EXACT_FILENAME[30]   = "exact_junctions.bed";
const char SPLICE_EXTEND_FILENAME[30]  = "extend_junctions.bed";

//========================================================================


options_t *options_new(void) {
	options_t *options = (options_t*) calloc (1, sizeof(options_t));

	/*options->in_filename = "";
	options->bwt_dirname = "";
	options->genome_filename = "";
	options->chromosome_filename = "";*/
	options->output_filename = strdup(DEFAULT_OUTPUT_FILENAME);
	options->splice_exact_filename = strdup(SPLICE_EXACT_FILENAME);
	options->splice_extend_filename = strdup(SPLICE_EXTEND_FILENAME);
	options->num_gpu_threads = DEFAULT_GPU_THREADS;
	options->num_cpu_threads = DEFAULT_CPU_THREADS;
	options->rna_seq = 0; 
	options->min_cal_size = DEFAULT_MIN_CAL_SIZE; 
	options->cal_seeker_errors = DEFAULT_CAL_SEEKER_ERROS;
	options->seeds_max_distance = DEFAULT_SEEDS_MAX_DISTANCE;
	options->bwt_threads = DEFAULT_BWT_THREADS;
	options->batch_size = DEFAULT_BATCH_SIZE;
	options->write_size = DEFAULT_WRITE_SIZE;
	options->num_cal_seekers = DEFAULT_NUM_CAL_SEEKERS;
	options->region_threads = DEFAULT_REGION_THREADS;
	options->num_sw_servers = DEFAULT_NUM_SW_THREADS;
	options->min_seed_size = DEFAULT_MIN_SEED_SIZE;
	options->seed_size = DEFAULT_SEED_SIZE;
	options->max_intron_length = DEFAULT_MAX_INTRON_LENGTH;
	options->flank_length = DEFAULT_FLANK_LENGTH;
	options->min_score = DEFAULT_SW_MIN_SCORE;
	options->match = DEFAULT_SW_MATCH;
	options->mismatch = DEFAULT_SW_MISMATCH;
	options->gap_open = DEFAULT_SW_GAP_OPEN;
	options->gap_extend = DEFAULT_SW_GAP_EXTEND;
	options->min_intron_length = DEFAULT_MIN_INTRON_LENGTH;
	options->timming = 0;
	options->statistics = 0;
	//	options->help = DEFAULT_HELP;
	
	return options;
}


void options_free(options_t *options) {
	if(options == NULL) {
		return;
	}

	if (options->splice_exact_filename != NULL)	{ free(options->splice_exact_filename); }
	if (options->splice_extend_filename  != NULL)	{ free(options->splice_extend_filename); }
	if (options->in_filename  != NULL)		{ free(options->in_filename); }
	if (options->bwt_dirname  != NULL)		{ free(options->bwt_dirname); }
	if (options->genome_filename  != NULL)		{ free(options->genome_filename); }
	if (options->chromosome_filename  != NULL)	{ free(options->chromosome_filename); }
	if (options->output_filename  != NULL)		{ free(options->output_filename); }
	
	free(options);
}


void options_display(options_t *options) {

  char* in_filename = strdup(options->in_filename);
  char* bwt_dirname =  strdup(options->bwt_dirname);
  char* genome_filename =  strdup(options->genome_filename);
  char* chromosome_filename =  strdup(options->chromosome_filename);
  char* output_filename =  strdup(options->output_filename);
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
  unsigned int max_intron_length =  (unsigned int)options->max_intron_length;
  unsigned int flank_length =  (unsigned int)options->flank_length;
  float min_score =  (float)options->min_score;
  float match =   (float)options->match;
  float mismatch =   (float)options->mismatch;
  float gap_open =   (float)options->gap_open;
  float gap_extend =   (float)options->gap_extend;
  char *splice_exact_filename =  strdup(options->splice_exact_filename);
  char *splice_extend_filename =  strdup(options->splice_extend_filename);
  unsigned int min_intron_length =  options->min_intron_length;
    
  printf("PARAMETERS CONFIGURATION\n");
  printf("=================================================\n");
  printf("Num gpu threads %d\n", num_gpu_threads);
  printf("Num cpu threads %d\n",  num_cpu_threads);
  printf("RNA Server: %s\n",  rna_seq == 0 ? "Disable":"Enable");
  printf("CAL seeker errors: %d\n",  cal_seeker_errors);
  printf("Min CAL size: %d\n",  min_cal_size);
  printf("Seeds max distance: %d\n",  seeds_max_distance);
  printf("Batch size: %dBytes\n",  batch_size);
  printf("Write size: %dBytes\n",  write_size);
  printf("BWT Threads: %d\n",  bwt_threads);
  printf("Region Threads: %d\n",  region_threads);
  printf("Num CAL seekers: %d\n", num_cal_seekers);
  printf("Num SW servers: %d\n",  num_sw_servers);
  printf("Min seed size: %d\n",  min_seed_size);
  printf("Seed size: %d\n",  seed_size);
  printf("Max intron length: %d\n", max_intron_length);
  printf("Min intron length: %d\n", min_intron_length);
  printf("Flank length: %d\n", flank_length);
  printf("SMITH-WATERMAN PARAMETERS\n");
  printf("\tMin score  : %0.4f\n", min_score);
  printf("\tMatch      : %0.4f\n", match);
  printf("\tMismatch   : %0.4f\n", mismatch);
  printf("\tGap open   : %0.4f\n", gap_open);
  printf("\tGap extend : %0.4f\n", gap_extend);
  printf("=================================================\n");

  free(in_filename);
  free(bwt_dirname);
  free(genome_filename);
  free(chromosome_filename);
  free(output_filename);
  free(splice_exact_filename);
  free(splice_extend_filename);
}

//--------------------------------------------------------------------

void** argtable_options_new(void) {
	void **argtable = (void**)malloc((NUM_OPTIONS + 1) * sizeof(void*));	// NUM_OPTIONS +1 to allocate end structure

	// NOTICE that order cannot be changed as is accessed by index in other functions
	argtable[0] = arg_file1("i", "fq,fastq", NULL, "Reads file input");
	argtable[1] = arg_file1("b", "bwt", NULL, "BWT directory name");
	argtable[2] = arg_file1("g", "genome", NULL, "Genome filename");
	argtable[3] = arg_file1("c", "chromosome", NULL, "Chromosome filename");
	argtable[4] = arg_file0("m", "match-output", NULL, "Match output filename");
	argtable[5] = arg_int0(NULL, "gpu-threads", NULL, "Number of GPU Threads");
	argtable[6] = arg_int0(NULL, "cpu-threads", NULL, "Number of CPU Threads");
	argtable[7] = arg_lit0(NULL, "rna-seq", "Active RNA Seq");
	argtable[8] = arg_int0(NULL, "cal-seeker-errors", NULL, "Number of errors in CAL Seeker");
	argtable[9] = arg_int0(NULL, "min-cal-size", NULL, "Minimum CAL size");
	argtable[10] = arg_int0(NULL, "max-distance-seeds", NULL, "Maximum distance between seeds");
	argtable[11] = arg_int0(NULL, "batch-size", NULL, "Batch Size");
	argtable[12] = arg_int0(NULL, "write-size", NULL, "Write Size");
	argtable[13] = arg_int0(NULL, "num-cal-seekers", NULL, "Number of CAL Seekers");
	argtable[14] = arg_int0(NULL, "num-sw-servers", NULL, "Number of Smith-Waterman servers");
	argtable[15] = arg_int0(NULL, "num-bwt-threads", NULL, "Number of BWT threads");
	argtable[16] = arg_int0(NULL, "num-region-threads", NULL, "Number of region threads");
	argtable[17] = arg_int0(NULL, "seed-size", NULL, "Number of nucleotides in a seed");
	argtable[18] = arg_int0(NULL, "min-seed-size", NULL, "Minimum number of nucleotides in a seed");
	argtable[19] = arg_int0(NULL, "cal-flank-length", NULL, "Flank length for CALs");
	argtable[20] = arg_dbl0(NULL, "match", NULL, "Match value for Smith-Waterman algorithm");
	argtable[21] = arg_dbl0(NULL, "mismatch", NULL, "Mismatch value for Smith-Waterman algorithm");
	argtable[22] = arg_dbl0(NULL, "gap-open", NULL, "Gap open penalty for Smith-Waterman algorithm");
	argtable[23] = arg_dbl0(NULL, "gap-extend", NULL, "Gap extend penalty for Smith-Waterman algorithm");
	argtable[24] = arg_dbl0(NULL, "min-sw-score", NULL, "Minimum score for valid mappings");
	argtable[25] = arg_int0(NULL, "max-intron-length", NULL, "Maximum intron length");
	argtable[26] = arg_int0(NULL, "min-itron-length", NULL, "Minimum intron length");
	argtable[27] = arg_lit0("t", "timing", "Timming mode active");
	argtable[28] = arg_lit0("s", "statistics", "Statistics mode active");
	argtable[29] = arg_lit0("h", "help", "Help option");
	argtable[30] = arg_file0(NULL, "splice-exact", NULL, "Splice Junctions exact filename");
	argtable[31] = arg_file0(NULL, "splice-extend", NULL, "Splice Junctions extend filename");
	
	argtable[32] = arg_end(20);

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
  //	options_t *options = (options_t*) calloc (1, sizeof(options_t));
	
  if (((struct arg_file*)argtable[0])->count) { options->in_filename = strdup(*(((struct arg_file*)argtable[0])->filename)); }
  if (((struct arg_file*)argtable[1])->count) { options->bwt_dirname = strdup(*(((struct arg_file*)argtable[1])->filename)); }
  if (((struct arg_file*)argtable[2])->count) { options->genome_filename = strdup(*(((struct arg_file*)argtable[2])->filename)); }
  if (((struct arg_file*)argtable[3])->count) { options->chromosome_filename = strdup(*(((struct arg_file*)argtable[3])->filename)); }
  if (((struct arg_file*)argtable[4])->count) { free(options->output_filename); options->output_filename = strdup(*(((struct arg_file*)argtable[4])->filename)); }
  if (((struct arg_int*)argtable[5])->count) { options->num_gpu_threads = *(((struct arg_int*)argtable[5])->ival); }
  if (((struct arg_int*)argtable[6])->count) { options->num_cpu_threads = *(((struct arg_int*)argtable[6])->ival); }
  if (((struct arg_int*)argtable[7])->count) { options->rna_seq = (((struct arg_int *)argtable[7])->count); }
  if (((struct arg_int*)argtable[8])->count) { options->cal_seeker_errors = *(((struct arg_int*)argtable[8])->ival); }
  if (((struct arg_int*)argtable[9])->count) { options->min_cal_size = *(((struct arg_int*)argtable[9])->ival); }
  if (((struct arg_int*)argtable[10])->count) { options->seeds_max_distance = *(((struct arg_int*)argtable[10])->ival); }
  if (((struct arg_int*)argtable[11])->count) { options->batch_size = *(((struct arg_int*)argtable[11])->ival); }
  if (((struct arg_int*)argtable[12])->count) { options->write_size = *(((struct arg_int*)argtable[12])->ival); }
  if (((struct arg_int*)argtable[13])->count) { options->num_cal_seekers = *(((struct arg_int*)argtable[13])->ival); }
  if (((struct arg_int*)argtable[14])->count) { options->num_sw_servers = *(((struct arg_int*)argtable[14])->ival); }
  if (((struct arg_int*)argtable[15])->count) { options->bwt_threads = *(((struct arg_int*)argtable[15])->ival); }
  if (((struct arg_int*)argtable[16])->count) { options->region_threads = *(((struct arg_int*)argtable[16])->ival); }
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
  if (((struct arg_file*)argtable[30])->count) { free(options->splice_exact_filename); options->splice_exact_filename = strdup(*(((struct arg_file*)argtable[30])->filename)); }
  if (((struct arg_file*)argtable[31])->count) { free(options->splice_extend_filename); options->splice_extend_filename = strdup(*(((struct arg_file*)argtable[31])->filename)); }
  
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
  return options;
}

void usage(void **argtable) {
  printf("Usage:\n./main-cpu {qc | filter | prepro}");
  arg_print_syntaxv(stdout, argtable, "\n");
  arg_print_glossary(stdout, argtable, "%-50s\t%s\n");
}
