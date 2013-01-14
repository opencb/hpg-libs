#include "rna_aligner.h"
#define NUM_SECTIONS_TIME 		7

void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     options_t *options) {

  int path_length = strlen(options->output_name);
  int extend_length = 0;
  if (options->extend_name) {
    extend_length = strlen(options->extend_name);
  }

  char *reads_results = (char *)calloc((60 + extend_length), sizeof(char));
  char *extend_junctions = (char *)calloc((60 + extend_length), sizeof(char));
  char *exact_junctions = (char *)calloc((60 + extend_length), sizeof(char));


  char *output_filename = (char *)calloc((path_length + extend_length + 60), sizeof(char));
  char *extend_filename = (char *)calloc((path_length + extend_length + 60), sizeof(char));
  char *exact_filename = (char *)calloc((path_length + extend_length + 60), sizeof(char));

  if (options->extend_name) {
    strcat(reads_results, "/");
    strcat(reads_results, options->extend_name);
    strcat(reads_results, "_alignments.bam");  

    strcat(extend_junctions, "/");
    strcat(extend_junctions, options->extend_name);
    strcat(extend_junctions, "_extend_junctions.bed");

    strcat(exact_junctions, "/");
    strcat(exact_junctions, options->extend_name);
    strcat(exact_junctions, "_exact_junctions.bed");
 
  } else {
    strcat(reads_results, "/alignments.bam");
    strcat(extend_junctions, "/extend_junctions.bed");
    strcat(exact_junctions, "/exact_junctions.bed");
  } 

  strcat(output_filename, options->output_name);
  strcat(output_filename, reads_results);
  free(reads_results);

  strcat(extend_filename, options->output_name);
  strcat(extend_filename, extend_junctions);
  free(extend_junctions);

  strcat(exact_filename, options->output_name);
  strcat(exact_filename, exact_junctions);
  free(exact_junctions);


  //************** Set Threads to sections **************//
  size_t cpu_threads = options->num_cpu_threads;

  if (!options->bwt_set &&
      !options->reg_set && 
      !options->cal_set &&
      !options->sw_set &&
      cpu_threads > 4) {
    LOG_DEBUG("Auto Thread configuration ...");
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
  }
  //****************************************************//
  LOG_DEBUG("Auto Thread Configuration Done !");

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
    
  }

  // display selected options
  LOG_DEBUG("Displaying options...\n");
  options_display(options);

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
  omp_set_nested(1);
  
  if (time_on) { 
    timing_start(MAIN_INDEX, 0, timing_p);
  }

  #pragma omp parallel sections num_threads(7)
  {
    
    #pragma omp section
    {
      fastq_batch_reader_input_t input;
      fastq_batch_reader_input_init(options->in_filename, options->in_filename2, 
				    options->pair_mode, options->batch_size, 
				    &read_list, &input);

      fastq_batch_reader_aligner(&input);

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

      #pragma omp parallel num_threads( options->num_sw_servers)
      {
	rna_server_omp_smith_waterman(&input, chromosome_avls);
      }

      write_chromosome_avls(chromosome_avls, &write_list,  
			    extend_filename, 
			    exact_filename,
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
      batch_writer_input_init( output_filename,  
			       exact_filename,  
			       extend_filename,  
			       &write_list, genome, &input);
      batch_writer2(&input);
      }
  }

  if (time_on) { 
    timing_stop(MAIN_INDEX, 0, timing_p);
  }


  free(output_filename);
  free(exact_filename);
  free(extend_filename);

}

//--------------------------------------------------------------------
