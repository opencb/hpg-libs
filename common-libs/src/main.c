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
//#include "rna_server.h"
#include "batch_aligner.h"
#include "statistics.h"

//-------------------------------------------------------
// constants
//-------------------------------------------------------

#define OPTIONS 			30
#define MIN_ARGC  			5
#define NUM_SECTIONS_TIME 		10
#define NUM_SECTIONS_STATISTICS 	6
#define NUM_SECTIONS_STATISTICS_SB	23
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
unsigned int id_splice = 1;
timing_t* timing_p = NULL;
int num_of_chromosomes;
statistics_t *statistics_p;

//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------

enum option_req { NO_APPEAR=0, APPEAR=1, NO_REQUIRED=2, REQUIRED=3 };
enum type_option { FILE_T=0, INTEGER_T=1, FLOAT_T=2, NONE_T=3, OFILE_T=4, STRING_T=5 };

typedef struct option_item {
  char *option;
  char *value;
  int id;
  int required;
  int type;
} option_item_t;

enum options_id{  READS_INPUT        = 0,  	BWT_DIRECTORY     = 1, 		GENOME_FILE        = 2, 	CHROMOSOME_FILE        = 3, 
		  OUTPUT_FILE        = 4,  	MISMATCH_FILE     = 5, 		GPU_THREADS        = 6, 	CPU_THREADS            = 7, 
		  RNA_ALIGNER        = 8,  	SEED_SIZE         = 9, 		MIN_SEED_SIZE      = 10, 	SEEDS_MAX_DISTANCE     = 11,
		  BATCH_SIZE         = 12, 	CAL_SEEKERS_NUM   = 13,		SMITH_WATERMAN_NUM = 14,	BWT_NUM                = 15,
		  MIN_CAL_SIZE       = 16,	CAL_SEEKER_ERRORS = 17,		FLANK_LENGTH  	   = 18,	MATCH_VALUE            = 19,
		  MISMATCH_VALUE     = 20,   	GAP_OPEN_PENALTY  = 21,		GAP_EXTEND_PENALTY = 22,	MINIMUM_SCORE          = 23,
		  TIMING_MODE        = 24,	STATISTICS_MODE   = 25,		PAIR_INPUT         = 26,        PAIR_MODE              = 27,
                  PAIR_MIN_DISTANCE  = 28,      PAIR_MAX_DISTANCE = 29,          UNKNOWN 	   = -1
	      };

option_item_t long_options[OPTIONS] = {
  {"<Unknown>",			"<Unknown>",		                          		UNKNOWN,  	   	-1,    	    	-1},
  {"-i", 			"Reads File Input",	                           		READS_INPUT,  	 	REQUIRED,    	FILE_T},
  {"-p", 			"Reads Pair File Input",	                           	PAIR_INPUT,  	 	NO_REQUIRED,   	FILE_T},
  {"-b", 			"Bwt Directory Name",	                           		BWT_DIRECTORY, 	 	REQUIRED,    	FILE_T},
  {"-g", 			"Genome Filename",	                           		GENOME_FILE, 	 	NO_REQUIRED, 	FILE_T},
  {"-c", 			"Chromosome Filename",	                           		CHROMOSOME_FILE, 	NO_REQUIRED, 	FILE_T},
  {"-o", 			"Output Filename",                             		        OUTPUT_FILE, 	 	NO_REQUIRED, 	OFILE_T},
  {"--gpu-threads", 		"Number Gpu Threads",	                           		GPU_THREADS,     	NO_REQUIRED, 	INTEGER_T},
  {"--omp-threads", 		"Number CPU Threads",	                           		CPU_THREADS,     	NO_REQUIRED, 	INTEGER_T},
  {"--rna-seq", 		"",			                           		RNA_ALIGNER, 	 	NO_REQUIRED, 	NONE_T},
  {"--cal-seeker-errors", 	"Number of errors in cal seeker",                       	CAL_SEEKER_ERRORS, 	NO_REQUIRED, 	INTEGER_T},
  {"--min-cal-size", 	        "Minimum cal size",   	                          		MIN_CAL_SIZE,       	NO_REQUIRED, 	INTEGER_T},
  {"--max-distance-seeds", 	"Max Distance Segment",                           		SEEDS_MAX_DISTANCE,     NO_REQUIRED, 	INTEGER_T},
  {"--batch-size", 		"Batch Size",		                          		BATCH_SIZE, 		NO_REQUIRED, 	INTEGER_T},
  {"--num-cal-seekers", 	"Number of CAL seekers",                          		CAL_SEEKERS_NUM, 	NO_REQUIRED, 	INTEGER_T},
  {"--num-sw-servers", 		"Number of Smith-Waterman servers",               		SMITH_WATERMAN_NUM, 	NO_REQUIRED, 	INTEGER_T},
  {"--num-bwt-threads",	        "Number of bwt threads",               		                BWT_NUM, 	        NO_REQUIRED, 	INTEGER_T},
  {"--seed-size", 		"Seed size",                                      		SEED_SIZE, 		NO_REQUIRED, 	INTEGER_T},
  {"--min-seed-size",	        "Minimum number of nucleotides in a seed",              	MIN_SEED_SIZE, 		NO_REQUIRED, 	INTEGER_T},
  {"--cal-flank-length",	"Flank length for CALs",                          		FLANK_LENGTH, 		NO_REQUIRED, 	INTEGER_T},
  {"--match",	                "Match value for Smith-Waterman algorithm",       		MATCH_VALUE, 		NO_REQUIRED, 	FLOAT_T},
  {"--mismatch",	        "Mismatch value for Smith-Waterman algorithm",    		MISMATCH_VALUE, 	NO_REQUIRED, 	FLOAT_T},
  {"--gap-open",	        "Gap open penalty for Smith-Waterman algorithm",  		GAP_OPEN_PENALTY, 	NO_REQUIRED, 	FLOAT_T},
  {"--gap-extend",	        "Gap extend penalty for Smith-Waterman algorithm",		GAP_EXTEND_PENALTY, 	NO_REQUIRED, 	FLOAT_T},
  {"--min-score", 		"Minimum score for valid mappings (0..1)",        		MINIMUM_SCORE, 		NO_REQUIRED, 	FLOAT_T},
  {"--pair-mode", 		"Pair mode: single-end, paired-end or mate-pair",	        PAIR_MODE, 		NO_REQUIRED, 	STRING_T},
  {"--pair-min-distance", 	"Minimum distance between pairs (400)",        		        PAIR_MIN_DISTANCE, 	NO_REQUIRED, 	INTEGER_T},
  {"--pair-max-distance", 	"Maximum distance between pairs (1000)",       		        PAIR_MAX_DISTANCE, 	NO_REQUIRED, 	INTEGER_T},
  {"-t", 		        "",		                               			TIMING_MODE, 		NO_REQUIRED, 	NONE_T},
  {"-s", 		        "",		                               			STATISTICS_MODE,	NO_REQUIRED, 	NONE_T}  
}; 

option_item_t searchOptions(char *str, int *appearOptions) {
  register int i;
  
  for(i=0; i< OPTIONS; i++){
    if(strcmp(str, long_options[i].option) == 0){
      appearOptions[i]=APPEAR;
      return long_options[i];
    }
  }
  
  return long_options[0];
}

void usagePrint() {	
  register int i;
  printf("\nUsage:\n");
  printf("  ./%s [Options]", __FILE__);
  
  for(i=1; i<OPTIONS; i++){
    if( long_options[i].required == REQUIRED){
      if(strcmp("", long_options[i].value) != 0){
	printf(" %s <%s>", long_options[i].option, long_options[i].value);
      }else{
	printf(" %s", long_options[i].option);
      }
    }
  }
  
  printf("\n\nDescription \n");
  for(i=1; i<OPTIONS; i++){
    if(strcmp("", long_options[i].value) != 0){
      printf("  %-25s <%s>\n", long_options[i].option, long_options[i].value);
    }else{
      printf("  %-25s\n", long_options[i].option);
    }
  }
}

//===============================================================
//  main
//===============================================================
extern int unmapped_by_max_cals_counter[100];
extern int unmapped_by_score_counter[100];

int main(int argc, char* argv[]) {
  char *in_filename1 = NULL;
  char *in_filename2 = NULL;
  char *bwt_dirname;
  char *genome_filename;
  char *chromosome_filename;
  char *output_filename = "reads_results.bam";
  //char* mismatch_filename = "MismatchsResults.fq";
  unsigned int num_gpu_threads = 32;
  unsigned int num_cpu_threads = 6;
  unsigned int rna_seq = 0; 
  unsigned int cal_seeker_errors = 0; //select value default
  unsigned int min_cal_size = 20; //select value default
  unsigned int seeds_max_distance = 60; //select value default
  unsigned int bwt_threads = 4; //select value default
  unsigned int batch_size = 200000; //  0.2 MB
  unsigned int write_size = 500000;  //  0.5 MB
  unsigned int num_cal_seekers = 1;
  unsigned int num_sw_servers = 1;
  unsigned int min_seed_size = 16;
  unsigned int seed_size = 18;
  //unsigned int num_seeds_per_cal = 2;
  unsigned int flank_length = 20;
  float min_score = 0.60f;
  float match = 5.0f;
  float mismatch = -4.0f;
  float gap_open = 10.0f;
  float gap_extend = 0.5f;
  unsigned int version;
  char *splice_filename = "SpliceJunctions.bed";

  int is_pair = 0;
  int flags = SINGLE_END_FLAG;
  char *pair_mode = "single-end";
  unsigned int pair_min_distance = 400, pair_max_distance = 1000;
  pair_mng_t *pair_mng = NULL;

  //printf("A=%i C=%i G=%i N=%i T=%i\n", 'A'&7, 'C'&7, 'G'&7, 'N'&7, 'T'&7);
  //exit(-1);

  register int i;
  option_item_t optionData;
  /*
  if(argc < MIN_ARGC)  
  {
    printf("ERROR: No readInput file or BWT directory specified\n");
    usagePrint();
    exit(-1);
  }*/
  int *appearOptions;
  int valueAtoi;
  float valueAtof = 0.0f;
  appearOptions=(int *)calloc(OPTIONS, sizeof(int));
  
  for(i=1; i<argc; i++){
    optionData=searchOptions(argv[i], appearOptions);
    
    if(optionData.id == UNKNOWN){
	printf("ERROR: Unrecognized option %s\n", argv[i]);
	usagePrint();
	exit(-1);
    }
    
    else if( (strcmp("", optionData.value) != 0) ){

      if(i+1 >= argc){
	printf("ERROR: Usage %s <%s> \n", optionData.option, optionData.value);
	usagePrint();
	exit(-1);
      }
      
      i++;
      
      if (optionData.type == FILE_T) {
	if (fopen(argv[i],"r") == 0) {
	  printf("ERROR: File %s not found\n", argv[i]);
	  usagePrint();
	  exit(-1);
	}
      } else if (optionData.type == INTEGER_T) {
	sscanf(argv[i], "%i", &valueAtoi);
	if(valueAtoi < 0) {
	  printf("ERROR: All values are positive\n");
	  usagePrint();
	  exit(-1);
	}
      } else if (optionData.type == FLOAT_T) {
	//sscanf(argv[i], "%.2f", &valueAtof);
	valueAtof = atof(argv[i]);
      }
    }
	      
    switch(optionData.id) {
        case READS_INPUT:
	     in_filename1 = argv[i];
	     break;
        case PAIR_INPUT:
	     in_filename2 = argv[i];
	     break;
        case BWT_DIRECTORY:
	     bwt_dirname = argv[i];
	     break;
        case GENOME_FILE:
	     genome_filename = argv[i];
	     break;
        case CHROMOSOME_FILE:
	     chromosome_filename = argv[i];
	     break;
        case OUTPUT_FILE:
	     output_filename = argv[i];
	     break;
        case GPU_THREADS:
	     num_gpu_threads = valueAtoi;
	     break;
        case CPU_THREADS:
	     num_cpu_threads = valueAtoi;
	     break;
        case RNA_ALIGNER:
	     rna_seq = 1;
	     break;
        case CAL_SEEKER_ERRORS:
	     cal_seeker_errors = valueAtoi;
	     break;
        case MIN_CAL_SIZE:
	     min_cal_size = valueAtoi;
	     break;
        case MIN_SEED_SIZE:
	     min_seed_size = valueAtoi;
	     break;
        case BATCH_SIZE:
	     batch_size = valueAtoi;
	     break;
        case CAL_SEEKERS_NUM:
	     num_cal_seekers = valueAtoi;
	     break;
        case SMITH_WATERMAN_NUM:
	     num_sw_servers = valueAtoi;
	     break;
        case SEEDS_MAX_DISTANCE:
	     seeds_max_distance = valueAtoi;
	     break;
        case SEED_SIZE:
	     seed_size = valueAtoi;
	     break;
        case BWT_NUM:
	     bwt_threads = valueAtoi;
	     break;
        case FLANK_LENGTH:
	     flank_length = valueAtoi;
	     break;
        case MATCH_VALUE:
	     match = valueAtof;
	     break;
        case MISMATCH_VALUE:
	     mismatch = valueAtof;
	     break;
        case GAP_OPEN_PENALTY:
	     gap_open = valueAtof;
	     break;
        case GAP_EXTEND_PENALTY:
	     gap_extend = valueAtof;
	     break;
        case MINIMUM_SCORE:
	     min_score = valueAtof;
	     break;
        case PAIR_MODE:
	     pair_mode = argv[i];
	     if (strcmp(pair_mode, "paired-end") == 0) {
	       flags = PAIRED_END_FLAG;
	       is_pair = 1;
	     } else if (strcmp(pair_mode, "mate-pair") == 0) {
	       flags = MATE_PAIR_FLAG;
	       is_pair = 1;
	     } else {
	       flags = SINGLE_END_FLAG;
	       is_pair = 0;
	     }
	     break;
        case PAIR_MIN_DISTANCE:
	     pair_min_distance = valueAtoi;
	     break;
        case PAIR_MAX_DISTANCE:
	     pair_max_distance = valueAtoi;
	     break;
        case TIMING_MODE:
	     time_on = 1;
	     break;
	case STATISTICS_MODE:
	     statistics_on = 1;
	     break;
      }
  }
  
  for (i = 0; i < OPTIONS; i++) {
       if( (appearOptions[i] == NO_APPEAR) && (long_options[i].required == REQUIRED)) {
	    printf("ERROR: Parameters");
      
	    for (i = 0; i < OPTIONS; i++) {
		 if(long_options[i].required == REQUIRED) {
		      if(strcmp("", long_options[i].value) != 0) {
			   printf(" %s <%s>", long_options[i].option, long_options[i].value);
		      } else {
			   printf(" %s", long_options[i].option);
		      }
		 }
	    }
	    printf(" are required\n");
	    usagePrint();
	    exit(-1);
       }
  }
  free(appearOptions);
  
  int cuda;
  
  /*#ifdef CUDA_VERSION
     cuda = 1;
  #else
     cuda = 0;
  #endif*/

  printf("PARAMETERS CONFIGURATION\n");
  printf("=================================================\n");
  printf("Num gpu threads %d\n", num_gpu_threads);
  printf("Num cpu threads %d\n", num_cpu_threads);
  printf("RNA Server: %s\n", rna_seq == 0 ? "Disable":"Enable");
  printf("CAL seeker errors: %d\n", cal_seeker_errors);
  printf("Min CAL size: %d\n", min_cal_size);
  printf("Seeds max distance: %d\n", seeds_max_distance);
  printf("Batch size: %dBytes\n", batch_size);
  printf("Write size: %dBytes\n", write_size);
  printf("BWT Threads: %d\n", bwt_threads);
  printf("Num CAL seekers: %d\n", num_cal_seekers);
  printf("Num SW servers: %d\n", num_sw_servers);
  printf("Min seed size: %d\n", min_seed_size);
  printf("Seed size: %d\n", seed_size);
  printf("Flank length: %d\n", flank_length);
  printf("Pair mode: %s\n", pair_mode);
  printf("Min. distance between pairs: %d\n", pair_min_distance);
  printf("Max. distance between pairs: %d\n", pair_max_distance);
  printf("SMITH-WATERMAN PARAMETERS\n");
  printf("\tMin score  : %0.4f\n", min_score);
  printf("\tMatch      : %0.4f\n", match);
  printf("\tMismatch   : %0.4f\n", mismatch);
  printf("\tGap open   : %0.4f\n", gap_open);
  printf("\tGap extend : %0.4f\n", gap_extend);
  printf("=================================================\n");

  // timing
  if (time_on) { 
    char* labels_time[NUM_SECTIONS_TIME] = {"Initialization BWT index   ", 
					    "Initialization genome index", 
					    "FastqQ reader              ", 
					    "BWT server                 ", 
					    "Region Seeker              ", 
					    "CAL seeker                 ", 
					    "Rna server                 ", 
					    "Batch writer               ", 
					    "Free main memory           ", 
					    "Total time                 "};
    
    int num_threads[NUM_SECTIONS_TIME] = {1, 1, 1, 1, 1, num_cal_seekers, num_sw_servers, 1, 1, 1};
    timing_p = timing_new((char**) labels_time, (int*) num_threads, NUM_SECTIONS_TIME);
    
    timing_start(MAIN_INDEX, 0, timing_p);
    //timing_start(INIT_INDEX, 0, timing_p);
    
  }
  if(statistics_on) {
    /*
     * FastqQ reader      : Num Batches %d 
     * BWT server         : Num Reads process %d, Num reads mapped %d, Num reads unmapped %d, Num mappings report %d 
     * Region Seeker      : Total num regions %d
     * CAL seeker         : Total cals %d, Num Reads unmapped %d
     * Rna server         : Num reads mapped %d, Num reads unmapped %d 
     * Batch writer       : Total num mappings report %d
     * Total 		  : Total reads %d, Total reads mapped %d, Total reads unmapped %d, Total splice junctions %d 	
     */				
    char* labels_statistics[NUM_SECTIONS_STATISTICS] = {"FastqQ reader              ", 
							"BWT server                 ", 
							"Region seeker              ", 
							"CAL seeker                 ", 
							"Rna server                 ",
							"Total 	Statistics	    "};
    
    char *sub_labels_statistics[NUM_SECTIONS_STATISTICS_SB] = { "Num batches             ",
								"Total reads             ",
								
								"Num batches process     ",
								"Reads process           ",
								"   Reads mapped         ", 
								"   Reads unmapped       ",
								"Mappings report         ",
								
								"Num batches process     ",
								"Num regions             ",
								
								"Num batches process     ",
								"Total cals              ",
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
				      
    unsigned int num_values[NUM_SECTIONS_STATISTICS] = {2, 5, 2, 3, 7, 4};
    statistics_p = statistics_new((char **)labels_statistics, (char **)sub_labels_statistics, (unsigned int *)num_values, NUM_SECTIONS_STATISTICS, NUM_SECTIONS_STATISTICS_SB);
    
  }
  // initialize some structures: Burrow-Wheeler objects, genome, nucletotide table...
  // (all these initializations could be performed in parallel)
  //

  printf("Reading bwt index...\n");
  if(time_on){ timing_start(INIT_BWT_INDEX, 0, timing_p); }
  bwt_index_t *bwt_index_p = bwt_index_new(bwt_dirname);
  if(time_on){ timing_stop(INIT_BWT_INDEX, 0, timing_p); }
  printf("Reading bwt index done !!\n");
  // (errors, threads, max aligns) 
  bwt_optarg_t *bwt_optarg_p = bwt_optarg_new(1, bwt_threads, 200);
  //  bwt_optarg_t *bwt_optarg_p = bwt_optarg_new(1, bwt_threads, 500);
  
  
  //GOOD LUCK(20, 60, 18, 16, 0)
  //                  cal_optarg_new(min_cal_size, max_cal_distance, seed_size, min_seed_size, num_error)
  cal_optarg_t *cal_optarg_p = cal_optarg_new(min_cal_size, seeds_max_distance, seed_size, min_seed_size, cal_seeker_errors);
  
  printf("reading genome...\n");
  if(time_on){ timing_start(INIT_GENOME_INDEX, 0, timing_p); }
  genome_t* genome_p = genome_new(genome_filename, chromosome_filename);
  if(time_on){ timing_stop(INIT_GENOME_INDEX, 0, timing_p); }
  printf("reading genome done !!\n");

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

  void *context_p = NULL;


 /* if (cuda) {
    int num_gpus = 2;
    context_p = (void *) gpu_context_new(num_gpus, num_gpu_threads);
  }*/
  
  // lists to communicate/synchronize the different threads
  //
  list_t read_list;
  list_t unmapped_reads_list;
  list_t regions_list;
  list_t sw_list; 
  list_t pair_list;
  list_t write_list;

  list_init("read", 1, 5, &read_list);
  list_init("unmmaped reads", 1, 5, &unmapped_reads_list);
  list_init("regions", 1, 5, &regions_list);
  
  if (is_pair) {
    list_init("sw", num_cal_seekers, 5, &sw_list);
    list_init("pair", 1 + num_cal_seekers, 5, &pair_list);
    list_init("write", 1 + num_sw_servers, 5, &write_list);

    pair_mng = pair_mng_new(pair_mode, pair_min_distance,
			    pair_max_distance);
  } else {
    list_init("sw", num_cal_seekers, 5, &sw_list); 
    list_init("pair", 0, 5, &pair_list);
    list_init("write", 1 + num_sw_servers + num_cal_seekers, 
	      5, &write_list);
  }

  // in sw_list, the threads insert are the cal_seekers and
  // in case of pair mode, the bwt_server too
  /*
  printf("writers to read_list = %d\n", list_get_writers(&read_list));
  printf("writers to unmapped_reads_list = %d\n", list_get_writers(&unmapped_reads_list));
  printf("writers to regions_list = %d\n", list_get_writers(&regions_list));
  printf("writers to sw_list = %d\n", list_get_writers(&sw_list));
  printf("writers to pair_list = %d\n", list_get_writers(&pair_list));
  printf("writers to write_list = %d\n", list_get_writers(&write_list));
  */

  //start delete ...
  //list_decr_writers(&unmapped_reads_list);
  //list_decr_writers(&regions_list);
  //list_decr_writers(&sw_list);
  //list_decr_writers(&write_list);
  //list_decr_writers(&write_list);
  //list_decr_writers(&write_list);
  //end delete

  // launch threads in parallel
  //
  omp_set_nested(1);
  //omp_set_num_threads(num_cpu_threads);
  
  //Rna Version
  version = 2;


  for (int i = 0; i < num_cpu_threads; i++) {
    unmapped_by_max_cals_counter[i] = 0;
    unmapped_by_score_counter[i] = 0;
  }

  //*********************************************************************
  //*********************************************************************
  //*********************************************************************
  {

    list_t read_list, write_list;
    list_init("read", 1, 24, &read_list);
    list_init("write", num_cpu_threads, 24, &write_list);

    printf("writers to read_list = %d\n", list_get_writers(&read_list));
    printf("writers to write_list = %d\n", list_get_writers(&write_list));
    
    bwt_server_input_t bwt_input;
    bwt_server_input_init(NULL, 0, bwt_optarg_p, bwt_index_p, 
			  NULL, NULL, 0, NULL, &bwt_input);
    
    region_seeker_input_t region_input;
    region_seeker_input_init(NULL, cal_optarg_p, bwt_optarg_p, 
			     bwt_index_p, NULL, &region_input);
    
    cal_seeker_input_t cal_input;
    cal_seeker_input_init(NULL, cal_optarg_p, NULL, 0, 
			  NULL, NULL, &cal_input);
    
    pair_server_input_t pair_input;
    pair_server_input_init(pair_mng, NULL, NULL, NULL, &pair_input);
    
    sw_server_input_t sw_input;
    sw_server_input_init(NULL, NULL, 0, match, mismatch, 
			 gap_open, gap_extend, min_score, flank_length, genome_p, 
			 &sw_input);


    
    double t_total;
    struct timeval t1, t2;

    #pragma omp parallel sections num_threads(3 + num_cpu_threads)
    {
      printf("Principal Sections %d threads\n", omp_get_num_threads());

      // fastq batch reader
      #pragma omp section
      {
	fastq_batch_reader_input_t input;
	fastq_batch_reader_input_init(in_filename1, in_filename2, flags, batch_size, 
				      &read_list, &input);
	fastq_batch_reader(&input);
      }

      // batch aligner
      #pragma omp section
      {
	gettimeofday(&t1, NULL);
        #pragma omp parallel num_threads(num_cpu_threads)
	{
	  batch_aligner_input_t input;
	  batch_aligner_input_init(&read_list, &write_list,
				   &bwt_input, &region_input, &cal_input,
				   (is_pair ? &pair_input : NULL), &sw_input,
				   &input);
	  
	  batch_aligner(&input);
	}
	gettimeofday(&t2, NULL);
	printf("\n\naligner time (no IO): %0.5f sec\n\n", (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6);
      }
      // batch writer
      #pragma omp section
      {
	batch_writer_input_t input;
	batch_writer_input_init("new.sam", splice_filename, &write_list, &input);
	batch_writer2(&input);
      }
    }
  }

//  printf("main.c: aborting by debugging\n");
//  abort();

  size_t total_item = 0;
  double max_time = 0, total_throughput = 0;
  printf("\nBWT time:\n");
  for (int i = 0; i < num_cpu_threads; i++) {
    printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f BWT/s (reads)\n", 
	   i, bwt_time[i] / 1e6, thr_batches[i], thr_bwt_items[i], 1e6 * thr_bwt_items[i] / bwt_time[i]);
    total_item += thr_bwt_items[i];
    total_throughput += (1e6 * thr_bwt_items[i] / bwt_time[i]);
    if (max_time < bwt_time[i]) max_time = bwt_time[i];
  }
  printf("\n\tTotal BWTs: %d, Max time = %0.4f, Throughput = %0.2f BWT/s\n", total_item, max_time / 1e6, total_throughput);

  total_item = 0; max_time = 0; total_throughput = 0;
  printf("\nSeeding time:\n");
  for (int i = 0; i < num_cpu_threads; i++) {
    printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f BWT/s (seeds)\n", 
	   i, seeding_time[i] / 1e6, thr_batches[i], thr_seeding_items[i], 1e6 * thr_seeding_items[i] / seeding_time[i]);
    total_item += thr_seeding_items[i];
    total_throughput += (1e6 * thr_seeding_items[i] / seeding_time[i]);
    if (max_time < seeding_time[i]) max_time = seeding_time[i];
  }
  printf("\n\tTotal BWTs: %d, Max time = %0.4f, Throughput = %0.2f BWT/s\n", total_item, max_time / 1e6, total_throughput);

  total_item = 0; max_time = 0; total_throughput = 0;
  printf("\nCAL time:\n");
  for (int i = 0; i < num_cpu_threads; i++) {
    printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f CAL/s)\n", 
	   i, cal_time[i] / 1e6, thr_batches[i], thr_cal_items[i], 1e6 * thr_cal_items[i] / cal_time[i]);
    total_item += thr_cal_items[i];
    total_throughput += (1e6 * thr_cal_items[i] / cal_time[i]);
    if (max_time < cal_time[i]) max_time = cal_time[i];
  }
  printf("\n\tTotal CALs: %d, Max time = %0.4f, Throughput = %0.2f CAL/s\n", total_item, max_time / 1e6, total_throughput);

  total_item = 0; max_time = 0; total_throughput = 0;
  printf("\nSW time:\n");
  for (int i = 0; i < num_cpu_threads; i++) {
    printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f SW/s)\n", 
	   i, sw_time[i] / 1e6, thr_batches[i], thr_sw_items[i], 1e6 * thr_sw_items[i] / sw_time[i]);
    total_item += thr_sw_items[i];
    total_throughput += (1e6 * thr_sw_items[i] / sw_time[i]);
    if (max_time < sw_time[i]) max_time = sw_time[i];
  }
  printf("\n\tTotal SWs: %d, Max time = %0.4f, Throughput = %0.2f SW/s\n", total_item, max_time / 1e6, total_throughput);

  //*********************************************************************
  //*********************************************************************
  //*********************************************************************

    // free memory
    bwt_optarg_free(bwt_optarg_p);
    cal_optarg_free(cal_optarg_p);
    bwt_index_free(bwt_index_p);

    if (pair_mng != NULL) pair_mng_free(pair_mng);
    if (time_on) { timing_start(FREE_MAIN, 0, timing_p); }
    genome_free(genome_p);
    if (time_on) { timing_stop(FREE_MAIN, 0, timing_p); }
    
    /*if (cuda) {
      gpu_context_free((gpu_context_t*) context_p);
      }
      if (time_on) { timing_stop(FREE_INDEX, 0, timing_p); }*/
    
    if (time_on) { 
      timing_stop(MAIN_INDEX, 0, timing_p);
      timing_display(timing_p);
      timing_free(timing_p);
    }
    
    if(statistics_on) {
      statistics_display(statistics_p);
      statistics_free(statistics_p);
    }


    int by_cals = 0, by_score = 0;
    for (int i = 0; i < num_cpu_threads; i++) {
      by_cals += unmapped_by_max_cals_counter[i];
      by_score += unmapped_by_score_counter[i];
    }
    printf("unmapped by MAX_CALS = %d\n", by_cals);
    printf("unmapped by SW score = %d\n", by_score);
    
    return 0;
}


//-----------------------------------------------------
//-----------------------------------------------------
