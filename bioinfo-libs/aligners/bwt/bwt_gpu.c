#include "bwt_gpu.h"
#include "bwt.h"
#include "containers/array_list.h" 

//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------
double kl_time = 0.0, s_time = 0.0;



size_t bwt_map_exact_seed_batch_gpu(fastq_batch_t *batch,
				    bwt_optarg_t *bwt_optarg,
				    cal_optarg_t *cal_optarg,
				    bwt_index_t *index,
				    gpu_context_t *gpu_context,
				    array_list_t **mapping_list) {
  
  struct timeval start_time, end_time;

  

  size_t num_mappings = 0;
  alignment_t *alignment;
  region_t *region;

  size_t num_reads = batch->num_reads;
  size_t max_seeds = 0, num_seeds, min_seed_size;
  char *code_seqs = (char *) malloc(batch->data_size);


  char *cigar, *seq_dup, *header, *quality;
  size_t idx, key, error, pos;

  size_t l_aux, k_aux, len, offset, start;

  int found = 0;
  size_t num_mapped_reads = 0, num_unmapped_reads = 0;
  
  /*printf("Real solutions\n");
  array_list_t *test_list = array_list_new(num_reads + 1, 1.25f, 
			     COLLECTION_MODE_SYNCHRONIZED);
  */
  /* char *seq = (char *)malloc(sizeof(char)*300);
  size_t seq_len = batch->data_indices[1] -  batch->data_indices[0] - 1;

  memcpy(seq, batch->seq, seq_len);
  seq[seq_len] = '\0';
  num_mappings = bwt_map_exact_seeds_seq(seq, 
					 cal_optarg->seed_size, 
					 cal_optarg->min_seed_size, 
					 bwt_optarg, index, 
					 test_list);
  
  //printf("===================================\n");
  //printf("num_mappings = %d\n", num_mappings);
  
  //printf("seq: %s\n", seq);
      
  for (size_t i = 0; i < num_mappings; i++) {
    printf("mapping: %i\n", i);
    region_t *alignment = array_list_get(i, test_list);
    printf("@DNA chromosome:%i | strand:%i | position:%i\n", 
	   alignment->chromosome_id,
	   alignment->strand,
	   alignment->start);
  }
  printf("===================================\n\n");
  
  printf("---------------\n");
  */

  replaceBases(batch->seq, code_seqs, batch->data_size);

  
  //printf("Reads %i\n", num_reads);
  for (int r = 0; r < num_reads; r++) {
    len = batch->data_indices[r+1] - batch->data_indices[r] - 1;
    num_seeds = len / cal_optarg->seed_size;
    min_seed_size = len % cal_optarg->seed_size;
    
    if (min_seed_size > 0) {
      num_seeds++;
    }
    //FOR EXACT PROCESS NEED PROCESS NEXT 'PASADA :)'
    //if ( min_seed_size > cal_optarg->min_seed_size) {
    //num_seeds++;
    //}

    if (num_seeds > max_seeds) {
      max_seeds = num_seeds;
    }
  }

  size_t num_kls = 2 * (num_reads * (2 * max_seeds));
  //printf("kl length = %i\n", (2 * num_reads)*(2 * max_seeds));
  size_t *k_values = (size_t*) calloc(num_kls, sizeof(size_t));
  size_t *l_values = (size_t*) calloc(num_kls, sizeof(size_t));
  unsigned int offset_strand = 0;
  unsigned int offset_pasada = 0;
  unsigned int offset_read = 0;

  /*
  // for debugging
  printf("num_reads = %d\n", num_reads);
  for (int r = 0; r < num_reads; r++) {
    printf("\nread #%d\n", r);
    for (int i =  batch->data_indices[r]; i < batch->data_indices[r+1]; i++) {
      printf("%c", batch->seq[i]);
    }
    printf("\n");
    for (int i =  batch->data_indices[r]; i < batch->data_indices[r+1]; i++) {
      printf("%i", code_seqs[i]);
    }  
  }

  printf("\n");

  printf("(k, l) : before\n");
  for (int i = 0; i < num_kls; i++) {
    printf("%i :(k, l) : (%lu, %lu)\n", i, k_values[i], l_values[i]);
    if (i == num_reads) printf("---\n");
  }
  */


  start_timer(start_time);
  //printf("Calculating kl values in gpu...\n");
  gpu_get_kl_values(cal_optarg->seed_size, cal_optarg->min_seed_size, 
		    max_seeds, num_reads, batch->data_size,
		    code_seqs, batch->data_indices,
		    gpu_context, k_values, l_values, SEED_MODE);
  
  stop_timer(start_time, end_time, kl_time);
  num_mappings = 0;
  start_timer(start_time);
  //printf("Calculating kl values in gpu END\n");

  //printf("Process kl in cpu...\n");
  // for debugging
  //#pragma omp parallel for private(offset_strand, offset_read, offset_pasada, len, k_aux, l_aux, key, start, region, idx) reduction(+:unm_mappings)
  for (unsigned int r = 0; r < num_reads; r++) {
    //printf("Read %i :: (k, l) : after searching\n", r);
    offset_strand = 0;
    offset_read = r * (2*max_seeds);
    for (unsigned int s = 0; s < 2; s++) {
      //printf(" Strand %i\n", s);
      offset_pasada = 0;
      for(unsigned int p = 0; p < 2; p++) {
	//printf("\t\tPasada %i\n", p);
	for(unsigned int v = 0; v < max_seeds; v++) { 
	  /*printf("\t\tPos: %i :(k, l) : (%lu, %lu)\n", 
		 offset_read + v + offset_strand + offset_pasada,
	  	 k_values[offset_read + v + offset_strand + offset_pasada], 
	  	 l_values[offset_read + v + offset_strand + offset_pasada]);
	  */
	  /**************************** Process KL Values ***************************/
	  
	  len = batch->data_indices[r + 1] - batch->data_indices[r] - 1;
	  k_aux =  k_values[offset_read + v + offset_strand + offset_pasada];
	  l_aux =  l_values[offset_read + v + offset_strand + offset_pasada];
	  
	  if (l_aux - k_aux + 1 < bwt_optarg->max_alignments_per_read) {
	    for (size_t j = k_aux; j <= l_aux; j++) {
	      if (index->S.ratio == 1) {
		key = index->S.vector[j];
	      } else {
		key = getScompValue(j, &index->S, &index->h_C, &index->h_O);
	      }
	      //printf("----> key value: %d\n", key);
	      
	      idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
	      //printf("----> idx value: %d\n", idx);
	      //chromosome = index->karyotype.chromosome + (idx-1) * IDMAX;
	      
	      if(key + len <= index->karyotype.offset[idx]) {
		start = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
		
		//printf("\tStrand:%i\tchromosome:%d\tStart:%u\tLength:%u\n",
		//       (r > num_reads), idx, start, len);
		
		// save all into one alignment structure and insert to the list
		num_mappings++;
		region = region_new(idx, s, start, start + cal_optarg->seed_size - 1, 0, 0, 0);
		//printf("Insert region\n");
		if (!array_list_insert((void*) region, mapping_list[r])) {
		  printf("Error to insert item into array list\n");
		}
		//printf("Insert ok!\n");
	      }
	    }	    
	    /**************************************************************************/
	  }
	  
	}
	offset_pasada = max_seeds;	
      }
      offset_strand =  num_reads * (2 * max_seeds);
    }
  }

  stop_timer(start_time, end_time, time_search_seed);  

  //  printf("Process kl END\n");
  free(code_seqs);
  free(k_values);
  free(l_values);

  return num_mappings;
}



size_t bwt_map_exact_batch_gpu(fastq_batch_t *batch,
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       gpu_context_t *gpu_context,
			       fastq_batch_t *unmapped_batch,
			       array_list_t *mapping_list) {

  struct timeval start_time, end_time;

  size_t num_mappings = 0;
  alignment_t *alignment;

  size_t num_reads = batch->num_reads;
  size_t num_kls = 2 * num_reads;
  size_t *k_values = (size_t*) calloc(num_kls, sizeof(size_t));
  size_t *l_values = (size_t*) calloc(num_kls, sizeof(size_t));

  char *code_seqs = (char *) malloc(batch->data_size);
  replaceBases(batch->seq, code_seqs, batch->data_size);

  unmapped_batch->num_reads = 0;

  /*
  // for debugging
  printf("num_reads = %d\n", num_reads);
  for (int r = 0; r < num_reads; r++) {
    printf("\nread #%d\n", r);
    for (int i =  batch->data_indices[r]; i < batch->data_indices[r+1]; i++) {
      printf("%c", batch->seq[i]);
    }
    printf("\n");
    for (int i =  batch->data_indices[r]; i < batch->data_indices[r+1]; i++) {
      printf("%i", code_seqs[i]);
    }  
  }

  printf("\n");

  printf("(k, l) : before\n");
  for (int i = 0; i < num_kls; i++) {
    printf("%i :(k, l) : (%lu, %lu)\n", i, k_values[i], l_values[i]);
    if (i == num_reads) printf("---\n");
  }
  */

  start_timer(start_time);

  gpu_get_kl_values(0, 0, 0, num_reads, batch->data_size,
		    code_seqs, batch->data_indices,
		    gpu_context, k_values, l_values, NORMAL_MODE);

  stop_timer(start_time, end_time, kl_time);

  /*
  // for debugging
  printf("(k, l) : after searching\n");
  for (int i = 0; i < num_kls; i++) {
    printf("%i :(k, l) : (%lu, %lu)\n", i, k_values[i], l_values[i]);
    if (i == num_reads) printf("---\n");
  }
  */

  char *cigar, *seq_dup, *header, *quality;
  size_t idx, key, error, pos;

  size_t l_aux, k_aux, len, offset, start;

  int found = 0;
  size_t num_mapped_reads = 0, num_unmapped_reads = 0;

  start_timer(start_time);

  for (int r = 0; r < num_reads; r++) {
    found = 0;
    len = batch->data_indices[r+1] - batch->data_indices[r] - 1;

    offset = 0;
    for (int i = 0; i < 2; i++) {
      
      k_aux = k_values[r + offset];
      l_aux = l_values[r + offset];
      offset = num_reads;

      if (l_aux - k_aux + 1 < bwt_optarg->max_alignments_per_read) {
	
	for (size_t j = k_aux; j <= l_aux; j++) {
	  if (index->S.ratio == 1) {
	    key = index->S.vector[j];
	  } else {
	    key = getScompValue(j, &index->S, &index->h_C, &index->h_O);
	  }
	  //printf("----> key value: %d\n", key);
	  
	  idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
	  //printf("----> idx value: %d\n", idx);
	  //chromosome = index->karyotype.chromosome + (idx-1) * IDMAX;
	  
	  if(key + len <= index->karyotype.offset[idx]) {
	    found = 1;
	    start = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	    
	    //	    printf("\tStrand:%i\tchromosome:%d\tStart:%u\tLength:%u\n",
	    //		   (r > num_reads), idx, start, len);
	    
	    cigar = (char *) malloc (sizeof(char) * len);
	    sprintf(cigar, "%d=\0", len);

	    seq_dup = (char *) malloc(sizeof(char) * (len + 1));
	    memcpy(seq_dup, &batch->seq[batch->data_indices[r]], len + 1);
	    quality = strndup(&batch->quality[batch->data_indices[r]], len + 1);
	    header = strndup(&batch->header[batch->header_indices[r]], batch->header_indices[r+1] - batch->header_indices[r] + 1);

	    // save all into one alignment structure and insert to the list
	    alignment = alignment_new();
	    alignment_init_single_end(header, seq_dup, quality, (r > num_reads), 
				      idx - 1,				      
				      start, 
				      cigar, 1, 255, 1, (num_mappings > 0), alignment);	    
	    
	    if (!array_list_insert((void*) alignment, mapping_list)) {
	      printf("Error to insert item into array list\n");
	    }
	    
	    num_mappings++;
	  }
	}
      }
    }
    if (!found) {
      //printf("%s\n", &batch->seq[batch->data_indices[r]]);
      unmapped_batch->num_reads++;
      //      exit(-1);
    }
  }
  stop_timer(start_time, end_time, s_time);


  free(code_seqs);
  free(k_values);
  free(l_values);

  return num_mappings;
}

//-----------------------------------------------------------------------------

/*
size_t bwt_map_exact_batch(fastq_batch_t *batch,
			   bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, 
			   fastq_batch_t *unmapped_batch,
			   array_list_t *mapping_list) {
  
  //size_t length_header, length_seq;
  //fastq_read_t *fq_read;
  //char header[1024], seq[1024], quality[1024]; 
  size_t num_mappings, num_mappings_tot; //num_unmapped;
  size_t read_pos;
  size_t num_threads = bwt_optarg->num_threads;
  //unsigned int th_id;
  size_t num_reads = batch->num_reads;
  //size_t batch_individual_size = batch->data_size / num_threads;
  size_t chunk = MAX(1, num_reads/(num_threads*10));

  array_list_t **individual_mapping_list_p   = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  //array_list_t **individual_unmapping_list_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  
  alignment_t *mapping;
  
  //struct timeval start_time, end_time;
  //double timer, parallel_t, sequential_t;
  unsigned int j, header_len, quality_len;
  size_t read_id = 0;
  for (unsigned int i = 0; i < num_reads; i++) {
    individual_mapping_list_p[i] = array_list_new(bwt_optarg->max_alignments_per_read, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
    //individual_unmapping_list_p[i] = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }  
  
  //printf("%d Threads %d chunk\n", num_threads, chunk);
  //num_mappings = 0;
  omp_set_num_threads(num_threads);
  
  //  start_timer(start_time);
  //printf("Process %d reads\n", num_reads);
  #pragma omp parallel for schedule(dynamic, chunk)
  for (unsigned int i = 0; i < num_reads; i++) {
 
    bwt_map_exact_seq(&(batch->seq[batch->data_indices[i]]), 
		      bwt_optarg, index, 
		      individual_mapping_list_p[i]);
  }

  //printf("Out mappings = %d\n", num_mappings_tot);

  read_pos = 0;
  unmapped_batch->header_indices[read_pos] = 0;
  unmapped_batch->data_indices[read_pos] = 0;
  for (unsigned int i = 0; i < num_reads; i++) {
    num_mappings = array_list_size(individual_mapping_list_p[i]);
    if(num_mappings){
      for (unsigned int j = 0; j < num_mappings; j++) {
	mapping = (alignment_t *) array_list_get(j, individual_mapping_list_p[i]);
	if(mapping != NULL){
	  header_len = batch->header_indices[i + 1] - batch->header_indices[i];
	  mapping->query_name = (char *)malloc(sizeof(char)*header_len);
	  get_to_first_blank(&(batch->header[batch->header_indices[i]]), header_len, mapping->query_name);
	  //printf("--->header id %s to %s\n",read->id, header_id);
	  
	  quality_len = batch->data_indices[i + 1] - batch->data_indices[i];
	  mapping->quality = (char *)malloc(sizeof(char)*quality_len);
	  memcpy(mapping->quality, &(batch->quality[batch->data_indices[i]]), quality_len);
	  //array_list_insert( mapping, mapping_list);
	  array_list_insert( mapping, mapping_list);
	}else{
	  printf("Error to extract item\n");
	}
      }
    }else{
      //fq_read = (fastq_read_t *)array_list_get(j, individual_unmapping_list_p[i]);
      read_pos++;
      header_len = batch->header_indices[i + 1] - batch->header_indices[i];
      quality_len = batch->data_indices[i + 1] - batch->data_indices[i];
      
      memcpy(&(unmapped_batch->header[unmapped_batch->header_indices[read_pos - 1]]), &(batch->header[batch->header_indices[i]]), header_len);
      memcpy(&(unmapped_batch->seq[unmapped_batch->data_indices[read_pos - 1]]),      &(batch->seq[batch->data_indices[i]]),      quality_len);
      memcpy(&(unmapped_batch->quality[unmapped_batch->data_indices[read_pos - 1]]),  &(batch->quality[batch->data_indices[i]]),  quality_len);
      
      unmapped_batch->data_indices[read_pos]   = unmapped_batch->data_indices[read_pos - 1]   + quality_len;
      unmapped_batch->header_indices[read_pos] = unmapped_batch->header_indices[read_pos - 1] + header_len;
      
    }
  }
  //printf("\tReads unmapped %d\n", read_pos);  
  unmapped_batch->num_reads = read_pos;

  for (unsigned int i = 0; i < num_reads; i++) {
    array_list_free(individual_mapping_list_p[i], NULL);
    //individual_unmapping_list_p[i] = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }  

  free(individual_mapping_list_p);
  //free(individual_unmapping_list_p);
  
  return array_list_size(mapping_list);


  /*size_t length_header, length_seq;
  fastq_read_t *fq_read;
  char header[1024], seq[1024], quality[1024]; 
  size_t num_mappings, num_unmapped;
  size_t read_pos;
  size_t num_threads = bwt_optarg->num_threads;
  short int th_id;
  size_t num_reads = batch->num_reads;
  size_t batch_individual_size = batch->data_size / num_threads;
  const size_t chunk = MAX(1, num_reads/(num_threads*10));
  
  array_list_t **individual_mapping_list_p   = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  array_list_t **individual_unmapping_list_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  
  alignment_t *mapping;
  
  struct timeval start_time, end_time;
  double timer, parallel_t, sequential_t;
  
  //#pragma omp parallel for  
  for (int i = 0; i < num_threads; i++) {
   individual_mapping_list_p[i]   = array_list_new(50000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
   individual_unmapping_list_p[i] = array_list_new(50000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }
  
  //printf("%d Threads %d chunk\n", num_threads, chunk);
  num_mappings = 0;
  omp_set_num_threads(num_threads);
  
  //start_timer(start_time);
  #pragma omp parallel for private(th_id, fq_read, num_mappings) schedule(dynamic, chunk)
  for (int i = 0; i < num_reads; i++) {
    th_id = omp_get_thread_num();
    
    fq_read = fastq_read_new(&(batch->header[batch->header_indices[i]]), 
			     &(batch->seq[batch->data_indices[i]]),
			     &(batch->quality[batch->data_indices[i]]));
    //fq_read = fastq_read_new(&header, &seq, &quality);
    //printf("Soy Thread %d and I process read %i : %s\n", th_id, i, fq_read->sequence);
    num_mappings = bwt_map_exact_read(fq_read, bwt_optarg, index, individual_mapping_list_p[th_id]);
    
    if (!num_mappings) {
      //Read not mapped. Store in unmapped_batch
      //printf("Insert Error Read\n");
      if (!array_list_insert((void *)fq_read, individual_unmapping_list_p[th_id])) {
	printf("Error to insert item into array list\n");
      }
    }else{
      fastq_read_free(fq_read);
    }
  }
  //stop_timer(start_time, end_time, timer);
  //parallel_t = timer;
  
  //start_timer(start_time);
  //printf("Join Results\n");
  
  //Join results
  for (int i = 0; i < num_threads; i++) {
    num_mappings = array_list_size(individual_mapping_list_p[i]);
    for (int j = 0; j < num_mappings; j++) {
      mapping = (alignment_t *)array_list_get(j, individual_mapping_list_p[i]);
      array_list_insert(mapping, mapping_list);
    }
  }

  //printf("Generate Batch\n");
  
  //Generate batch unmapped reads
  read_pos = 0;
  unmapped_batch->header_indices[read_pos] = 0;
  unmapped_batch->data_indices[read_pos] = 0;
  for (int i = 0; i < num_threads; i++) {
    num_unmapped = array_list_size(individual_unmapping_list_p[i]);
    for (int j = 0; j < num_unmapped; j++) {
      fq_read = (fastq_read_t *)array_list_get(j, individual_unmapping_list_p[i]);
      read_pos++;
      length_header = strlen(fq_read->id) + 1;
      length_seq = strlen(fq_read->sequence) + 1;
      memcpy(&(unmapped_batch->header[unmapped_batch->header_indices[read_pos - 1]]), fq_read->id,       length_header);
      memcpy(&(unmapped_batch->seq[unmapped_batch->data_indices[read_pos - 1]]),      fq_read->sequence, length_seq);
      memcpy(&(unmapped_batch->quality[unmapped_batch->data_indices[read_pos - 1]]),  fq_read->quality,  length_seq);
      
      unmapped_batch->data_indices[read_pos]   = unmapped_batch->data_indices[read_pos - 1]   + length_seq;
      unmapped_batch->header_indices[read_pos] = unmapped_batch->header_indices[read_pos - 1] + length_header;
    
      //fastq_read_free(fq_read);
    }
  }
  
  unmapped_batch->num_reads = read_pos;
  
  //printf("Free\n");
  
  for (int i = 0; i < num_threads; i++) {
    array_list_free(individual_unmapping_list_p[i], fastq_read_free);
    //array_list_free(individual_mapping_list_p[i], alignment_free);
  }

  free(individual_mapping_list_p);
  free(individual_unmapping_list_p);
  
  //stop_timer(start_time, end_time, timer);
  
  //sequential_t = timer - parallel_t;
  
  //global_parallel += parallel_t;
  //global_sequential += sequential_t;
  
  //printf("Parallel %d\n", parallel_t);
  
  //printf("%4.06f\t%4.06f\t%4.06f\n", parallel_t / 1000000, sequential_t / 1000000, timer / 1000000);
  
  return array_list_size(mapping_list);*/
/*
}
*/
//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------
/*
size_t bwt_map_exact_seed(char *seq, 
			  size_t seq_start, size_t seq_end,
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index, 
			  array_list_t *mapping_list){
  
  //printf("Process New Seeds\n");

  region_t *region;
  size_t start = 0;
  size_t end = seq_end - seq_start;
  result result;
  char *code_seq = &seq[seq_start];
  //(char *)malloc(sizeof(char)*(seq_end - seq_start + 1)); //= &seq[seq_start];
  //memcpy(code_seq, &seq[seq_start], seq_end - seq_start);
  //code_seq[seq_end - seq_start] = '\0';
  
  //size_t start = 0;
  //size_t end = seq_end - seq_start;
  size_t len = seq_end - seq_start;
  unsigned int num_mappings = 0;
  char plusminus[2] = "-+";
  size_t idx, key, direction, error, pos;
  results_list *r_list;
  //  result *r;
  size_t l_aux, k_aux;
  alignment_t *alignment;
  //  size_t len = strlen(seq);
  //  size_t start = 0;
  //  size_t end = len - 1;
  size_t start_mapping;
  
  
  for (short int type = 1; type >= 0; type--) {
       result.k = 0;
       result.l = index->h_O.siz - 2;
       result.start = start;
       result.end = end;
      if (type == 1) {
	result.pos = end;
	BWExactSearchBackward(code_seq, &index->h_C, &index->h_C1, &index->h_O, &result);
	//BWExactSearchBackward(code_seq, start, end, &index->h_C, &index->h_C1, &index->h_O, result_p);
      }else{
	result.pos = start;
	BWExactSearchForward(code_seq, &index->h_rC, &index->h_rC1, &index->h_rO, &result);
	//BWExactSearchForward(code_seq, start, end, &index->h_rC, &index->h_rC1, &index->h_rO, result_p);
      }
      k_aux = result.k;
      l_aux = result.l;
      if (l_aux - k_aux + 1 < bwt_optarg->max_alignments_per_read) {
	for (size_t j = k_aux; j <= l_aux; j++) {
	  if (index->S.ratio == 1) {
	    key = index->S.vector[j];
	  } else {
	    key = getScompValue(j, &index->S, &index->h_C, &index->h_O);
	  }
	  //printf("----> key value: %d\n", key);

	  idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
	  //printf("----> idx value: %d\n", idx);
	  //chromosome = index->karyotype.chromosome + (idx-1) * IDMAX;
 
	  if(key + len <= index->karyotype.offset[idx]) {
	    //start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	    //printf("Strand:%c\tchromosome:%s\tStart:%u\tend:%u\n",plusminus[type],
	    //	    index->karyotype.chromosome + (idx-1) * IDMAX,
	    //	    start_mapping, start_mapping + len);
	    //
	    start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	    // save all into one alignment structure and insert to the list
	    region = region_new(idx, !type, start_mapping, start_mapping + len);

	    if(!array_list_insert((void*) region, mapping_list)){
		  printf("Error to insert item into array list\n");
	    }
	    
	    num_mappings++;
	  }
	}
	//      } else {
	//	printf("****** too many mappings (%d) for this sequence (max. %d)!!\n", l_aux - k_aux + 1, bwt_optarg->max_alignments_per_read);
      }
  }
  //  free(result_p);
	
  return num_mappings;
  
}

size_t bwt_map_inexact_seed(char *seq, 
			    size_t seq_start, size_t seq_end,
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    array_list_t *mapping_list) {
  region_t *region; 

  char *code_seq = &seq[seq_start];
  size_t start = 0;
  size_t end = seq_end - seq_start;
  size_t len = end + 1;

  //  size_t len = strlen(seq);
  //  size_t start = 0;
  //  size_t end = len - 1;
  size_t start_mapping;
  //  char *code_seq = (char *) calloc(len, sizeof(char));
  //  replaceBases(seq, code_seq, len);

  // calculate vectors k and l
  size_t *k0 = (size_t *) malloc(len * sizeof(size_t));
  size_t *l0 = (size_t *) malloc(len * sizeof(size_t));
  size_t *k1 = (size_t *) malloc(len * sizeof(size_t));
  size_t *l1 = (size_t *) malloc(len * sizeof(size_t));
  size_t *ki0 = (size_t *) malloc(len * sizeof(size_t));
  size_t *li0 = (size_t *) malloc(len * sizeof(size_t));
  size_t *ki1 = (size_t *) malloc(len * sizeof(size_t));
  size_t *li1 = (size_t *) malloc(len * sizeof(size_t));
  
  BWExactSearchVectorBackward(code_seq, start, end, 0, index->h_O.m_count - 2,
			      k1, l1, &index->h_C, &index->h_C1, &index->h_O);

  BWExactSearchVectorForward(code_seq, start, end, 0, index->h_Oi.m_count - 2,
			     ki1, li1, &index->h_C, &index->h_C1, &index->h_Oi);

  BWExactSearchVectorForward(code_seq, start, end, 0, index->h_rO.m_count - 2,
			     k0, l0, &index->h_rC, &index->h_rC1, &index->h_rO);

  BWExactSearchVectorBackward(code_seq, start, end, 0, index->h_rOi.m_count - 2,
			      ki0, li0, &index->h_rC, &index->h_rC1, &index->h_rOi);
  

  // compare the vectors k and l to get mappings in the genome
  unsigned int num_mappings = 0;
  char plusminus[2] = "-+";
  int idx, key, direction, error, pos;
  //results_list *r_list;
  results_list r_list;
  result *r;
  alignment_t *alignment;

  for (int type = 1; type >= 0; type--) {

    //    r_list = (results_list *) calloc(1, sizeof(results_list));
    //new_results_list(r_list, 2000);
    new_results_list(&r_list, 2000);
    r_list.num_results = 0;
    //r_list->n = 0;
    r_list.read_index = 0;

    if (type == 1) {
      //      BWIterativeSearch1(code_seq, start, end, k1, l1, ki1, li1, 
      //			 &index->h_C, &index->h_C1, &index->h_O, &index->h_Oi, r_list);
      BWSearch1(code_seq, start, end, k1, l1, ki1, li1, 
		&index->h_C, &index->h_C1, &index->h_O, &index->h_Oi, &r_list);
    } else {
	//      BWIterativeSearch1(code_seq, start, end, ki0, li0, k0, l0, 
	//     			 &index->h_rC, &index->h_rC1, &index->h_rOi, &index->h_rO, r_list);
      BWSearch1(code_seq, start, end, ki0, li0, k0, l0, 
		&index->h_rC, &index->h_rC1, &index->h_rOi, &index->h_rO, &r_list);
    }

    for (size_t ii = 0; ii < r_list.num_results; ii++) {
      //for (size_t ii = 0; ii < r_list->n; ii++) {
      r = &r_list.list[ii];
      if (r->l - r->k + 1 < bwt_optarg->max_alignments_per_read) {
	for (unsigned int j = r->k; j <= r->l; j++) {
	  if (type) {
	    direction = r->dir;
	  } else {
	    direction = !r->dir;
	  }
	  if (index->S.ratio == 1) {
	    key = (direction)
	      ? index->Si.n - index->Si.vector[j] - len - 1
	      : index->S.vector[j];
	  } else {
	    key = (direction)
	      ? index->Si.n - getScompValue(j, &index->Si, &index->h_C,
					    &index->h_Oi) - len - 1
	      : getScompValue(j, &index->S, &index->h_C, &index->h_O);
	  }
	  idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
	  if(key + len <= index->karyotype.offset[idx]) {
	    //printf("\tvalue idx=%d\n", idx);
	    //printf("\t%s\t%c\t%s %u %s error: %s, pos: %i, base: %i ",
	   //   "nothing", plusminus[type],
	  //    index->karyotype.chromosome + (idx-1) * IDMAX,
	   //   index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]),
	   //   seq, bwt_error_type(r->err_kind[0]), r->position[0], r->base[0]);
	   //	  
	    start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	  // save all into one alignment structure and insert to the list
	    region = region_new(idx, !type, start_mapping, start_mapping + end);
	    
	    if(!array_list_insert((void*) region, mapping_list)){
	      printf("Error to insert item into array list\n");
	    }
	    num_mappings++;
	  }
	}
      }//end if max solution
    }
    //free_results_list(r_list);
    free(r_list.list);
    //    free(r_list);
  } // end for type 

  //free(code_seq);
  free(k0);
  free(l0);
  free(k1);
  free(l1);
  free(ki0);
  free(li0);
  free(ki1);
  free(li1);
  return num_mappings;

}
*/
//-----------------------------------------------------------------------------
/*
size_t bwt_map_inexact_seq(char *seq, 
			   bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, 
			   array_list_t *mapping_list) {
  

  //printf("\tIn function...\n");
  char *seq_dup;

  //return bwt_map_inexact_seq_by_pos(seq, start, end,
  //					bwt_optarg, index, mapping_list);
  
  size_t len = strlen(seq);
  size_t start = 0;
  size_t end = len - 1;

  char *code_seq = (char *) malloc(len * sizeof(char));

  replaceBases(seq, code_seq, len);

  // calculate vectors k and l
  size_t *k0 = (size_t *) malloc(len * sizeof(size_t));
  size_t *l0 = (size_t *) malloc(len * sizeof(size_t));
  size_t *k1 = (size_t *) malloc(len * sizeof(size_t));
  size_t *l1 = (size_t *) malloc(len * sizeof(size_t));
  size_t *ki0 = (size_t *) malloc(len * sizeof(size_t));
  size_t *li0 = (size_t *) malloc(len * sizeof(size_t));
  size_t *ki1 = (size_t *) malloc(len * sizeof(size_t));
  size_t *li1 = (size_t *) malloc(len * sizeof(size_t));

  BWExactSearchVectorBackward(code_seq, start, end, 0, index->h_O.siz - 2,
			      k1, l1, &index->h_C, &index->h_C1, &index->h_O);
  
  BWExactSearchVectorForward(code_seq, start, end, 0, index->h_Oi.siz - 2,
			     ki1, li1, &index->h_C, &index->h_C1, &index->h_Oi);
  
  BWExactSearchVectorForward(code_seq, start, end, 0, index->h_rO.siz - 2,
			     k0, l0, &index->h_rC, &index->h_rC1, &index->h_rO);
  
  BWExactSearchVectorBackward(code_seq, start, end, 0, index->h_rOi.siz - 2,
			      ki0, li0, &index->h_rC, &index->h_rC1, &index->h_rOi);
  

  // compare the vectors k and l to get mappings in the genome
  size_t num_mappings = 0;
  char plusminus[2] = "-+";
  size_t idx, key, direction;
  char error;
  int pos;
  //results_list *r_list;
  results_list r_list;
  result *r;
  alignment_t *alignment;
  char error_debug = 0;
  char *cigar_dup;
  char cigar[1024];
  size_t cigar_len, num_cigar_ops;
  error = MISMATCH;

  new_results_list(&r_list, bwt_optarg->max_alignments_per_read);

  for (int type = 1; type >= 0; type--) {

    r_list.num_results = 0;
    r_list.read_index = 0;

    //    printf("*** bwt.c: calling BWSearch1 with type = %d...\n", type);
    if (type == 1) {
      BWSearch1(code_seq, start, end, k1, l1, ki1, li1, 
		&index->h_C, &index->h_C1, &index->h_O, &index->h_Oi, &r_list);
    } else {
      BWSearch1(code_seq, start, end, ki0, li0, k0, l0, 
		&index->h_rC, &index->h_rC1, &index->h_rOi, &index->h_rO, &r_list);
    }
    //    printf("*** bwt.c: calling BWSearch1 with type = %d (num_results = %d). Done !!\n", type, r_list.num_results);
    
    for (size_t ii = 0; ii < r_list.num_results; ii++) {
      //for (size_t ii = 0; ii < r_list->n; ii++) {
     
      r = &r_list.list[ii];
      //printf("Errors number %d\n", r->num_mismatches);
      //if(r == NULL){printf("ERROR: Result list position null");exit(-1);}
      if(!r->num_mismatches)
	error = 0;
      else
	error = r->err_kind[0];

      pos = r->err_pos[0];
      //pos = r->position[0];
	  
      if (r->l - r->k + 1 < bwt_optarg->max_alignments_per_read) {
	//printf("\tk=%d - l=%d\n", r->k, r->l);      
	for (unsigned int j = r->k; j <= r->l; j++) {
	  if (type) {
	    direction = r->dir;
	  } else {
	    direction = !r->dir;
	  }
	  if (index->S.ratio == 1) {
	    key = (direction)
	      ? index->Si.siz - index->Si.vector[j] - len - 1
	      : index->S.vector[j];
	  } else {
	    key = (direction)
	      ? index->Si.siz - getScompValue(j, &index->Si, &index->h_C,
					    &index->h_Oi) - len - 1
	      : getScompValue(j, &index->S, &index->h_C, &index->h_O);
	  }
	  idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
	  if(key + len <= index->karyotype.offset[idx]) {
*/	    
	    /*alignment_p->flags = CODED_SEQ_FLAG;
	      alignment_p->header_p = &read_batch_p->header[read_batch_p->header_indices[i]];
	      alignment_p->read_p = seq;
	      alignment_p->quality_p = &read_batch_p->quality[read_batch_p->data_indices[i]];
	      
	      mapping_info_p = mapping_info_new(genome_p->chromosome + (index-1) * IDMAX,
	      plusminus[type],
	      genome_p->start[index-1] + (key - genome_p->offset[index-1]));
	      memcpy(&mapping_info_p->errors, r, sizeof(result));
	      mapping_item_p = list_item_new(0, 0, mapping_info_p);
	      list_insert_item(mapping_item_p, &alignment_p->mapping_list);*/
/*	    
	    // cigar processing
	    
	    //printf("Errors? %d\n", error);
	    if (error == 0) {
	      sprintf(cigar, "%i=\0", len);
	      num_cigar_ops = 1;
	      
	    } else if (error == MISMATCH) {
	      
	      if (pos == 0) {
		sprintf(cigar, "1S%iM\0", len-1);
		num_cigar_ops = 2;
	      } else if (pos == len - 1) {
		sprintf(cigar, "%iM1S\0", len-1);
		num_cigar_ops = 2;
	      } else {
		sprintf(cigar, "%iM\0", len);
		num_cigar_ops = 1;
	      }
	      
	    } else if (error == INSERTION) {
	      
	      if (pos == 0) {
		sprintf(cigar, "1H%iM\0", len-1);
		num_cigar_ops = 2;
	      } else if (pos == len - 1) {
		sprintf(cigar, "%iM1H\0", len-1);
		num_cigar_ops = 2;
	      } else {
		sprintf(cigar, "%iM1I%iM\0", pos, len - pos - 1);
		num_cigar_ops = 3;
	    }
	      
	    } else if (error == DELETION) {
	      
	      if (pos == 0) {
		sprintf(cigar, "1I%iM\0", len);
		num_cigar_ops = 2;
	      } else if (pos == len - 1) {
		sprintf(cigar, "%iM1I\0", len);
		num_cigar_ops = 2;
	      } else {
		sprintf(cigar, "%iM1D%iM\0", pos, len - pos);
		num_cigar_ops = 3;
	      }

	    }else{
	      printf("NUM MAPPINGS %d -> POS %d -> ERROR %d -> (%d):%s", num_mappings, pos, error, len, seq);
	      continue;
	      //exit(-1);
	      //error_debug = 1;
	    }
	    //printf("Value idx = %d\n", idx);
	    //  printf("\t%s\t%c\t%s %u %s error: %s, pos: %i, base: %i cigar: %s\n",
	    //  "nothing", plusminus[type],
	    //  index->karyotype.chromosome + (idx-1) * IDMAX,
	    //  index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]),
	    //  seq, bwt_error_type(r->err_kind[0]), r->position[0], r->base[0], cigar);
	      
	    if(!error_debug){
	      cigar_len = strlen(cigar) + 1;
	      cigar_dup = (char *)calloc(cigar_len, sizeof(char));
	      memcpy(cigar_dup, cigar, cigar_len);
	      
	      seq_dup = (char *)malloc(sizeof(char)*(len + 1));
	      memcpy(seq_dup ,seq, len + 1);
	    
	      // save all into one alignment structure and insert to the list
	      alignment = alignment_new();
	      alignment_init_single_end(NULL, seq_dup, NULL, !type, 
					idx - 1, //index->karyotype.chromosome + (idx-1) * IDMAX,
					index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]), 
					cigar_dup, num_cigar_ops, 255, 1, (num_mappings > 0), alignment);
	      
	      array_list_insert((void*) alignment, mapping_list);
	    
	      num_mappings++;
            }
	  }
	}//end for k and l
      }//end if max solutions
      //free(r);
    }
    //free_results_list(r_list);
    //    free(r_list);
  } // end for type 
  
  free(r_list.list);
  free(code_seq);

  free(k0);
  free(l0);
  free(k1);
  free(l1);
  free(ki0);
  free(li0);
  free(ki1);
  free(li1);
  //printf("\tOut function\n");

  return num_mappings;
}
*/
//-----------------------------------------------------------------------------
/*
size_t bwt_map_inexact_seqs_by_pos(char *seqs, 
				   unsigned int num_reads,
				   unsigned int starts,
				   unsigned int ends,
				   bwt_optarg_t *bwt_optarg, 
				   bwt_index_t *index, 
				   array_list_t *mapping_list) {

  unsigned int num_threads = bwt_optarg->num_threads;
  array_list_t **th_mapping_list = (array_list_t**) calloc(num_threads, sizeof(array_list_t*));
  
  const unsigned int chunk = MAX(1, num_reads / (num_threads * 10));
  unsigned int num_mappings = 0; 
  unsigned int num_mappings_tot = 0;
  short int th_id;
  alignment_t *mapping;
  
  for (int i = 0; i < num_threads; i++) {
    th_mapping_list[i] = array_list_new(100000, 1.25f, 
					COLLECTION_MODE_SYNCHRONIZED);
  }
  
  omp_set_num_threads(num_threads);
  
  //printf("Initialization ok!Process %d reads x %d threads\n", num_reads, num_theads);
  
  #pragma omp parallel for private(th_id) schedule(dynamic, chunk) reduction(+:num_mappings_tot)
  for (int i = 0; i < num_reads; i++) {
    th_id = omp_get_thread_num();
    printf("I'm thread %d and process read %d\n", th_id,  i);
    num_mappings_tot += bwt_map_inexact_seq(seqs[i], bwt_optarg, index, 
						th_mapping_list[th_id]);
  }
  
  printf("Process reads ok! Total mappings %d\n\n", num_mappings_tot);
  
  // mergeresults
  for (int i = 0; i < num_threads; i++) {
    num_mappings = array_list_size(th_mapping_list[i]);
   
    for (int j = 0; j < num_mappings; j++) {
      mapping = (alignment_t *)array_list_get(j, th_mapping_list[i]);
      array_list_insert( mapping, mapping_list);
    }
  }
  
  return num_mappings_tot;
}
*/	  
//-----------------------------------------------------------------------------

size_t bwt_map_exact_seeds_seq_gpue(char *seq, size_t seed_size, size_t min_seed_size,
				    bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				    array_list_t *mapping_list) {
  
/*
  size_t len = strlen(seq);
  size_t offset, num_seeds = len / seed_size;

  char *code_seq = (char *) calloc(len, sizeof(char));

  replaceBases(seq, code_seq, len);

  // first 'pasada'
  offset = 0;
  for (size_t i = 0; i < num_seeds; i++) {
    //printf("%s\n", &seq[offset]);
    //printf("1, seed %d: start = %d, end = %d\n", i, offset, offset + seed_size - 1);
    bwt_map_exact_seed(code_seq, offset, offset + seed_size - 1,
			 bwt_optarg, index, mapping_list);
    offset += seed_size;
  }

  // special processing for the last seed !!
  if (len % seed_size >= min_seed_size) {
    //printf("%s\n", &seq[offset]);
    //printf("1', : start = %d, end = %d\n", offset, len - 1);
    bwt_map_exact_seed(code_seq, offset, len - 1,
			 bwt_optarg, index, mapping_list);
  }

  // second 'pasada', shifting seeds by (seed_size / 2)
  offset = seed_size / 2;
  num_seeds = (len - seed_size / 2) / seed_size;
  for (size_t i = 0; i < num_seeds; i++) {
    //printf("%s\n", &seq[offset]);
    //printf("2, seed %d: start = %d, end = %d\n", i, offset, offset + seed_size - 1);
    bwt_map_exact_seed(code_seq, offset, offset + seed_size - 1,
			 bwt_optarg, index, mapping_list);
    offset += seed_size;
  }

  // again, special processing for the last seed !!
  if ((len - seed_size / 2) % seed_size >= min_seed_size) {
    //printf("%s\n", &seq[offset]);
    //printf("2',: start = %d, end = %d\n", offset, len - 1);
    bwt_map_exact_seed(code_seq, offset, len - 1,
			 bwt_optarg, index, mapping_list);
  }

  //printf("mapping list size = %d\n", array_list_size(mapping_list));
  free(code_seq);
  //  exit(-1);
*/

  return array_list_size(mapping_list);
  
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

