#include "rna_aligner.h"

//--------------------------------------------------------------------

void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     options_t *options) {
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
  #pragma omp parallel sections num_threads(7)
  {
    printf("Principal Sections %d threads\n", omp_get_num_threads());
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
