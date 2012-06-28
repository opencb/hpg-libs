#include "bwt.h"


//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_seq(char *seq, 
				   bwt_optarg_t *bwt_optarg, 
				   bwt_index_t *index, 
				   array_list_t *mapping_list);

unsigned int bwt_map_exact_read(fastq_read_t *read, 
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index, 
				    array_list_t *mapping_list);

unsigned int bwt_map_exact_seqs(char **seqs, 
				    unsigned int num_reads,
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index, 
				    array_list_t *mapping_list);

unsigned int bwt_map_exact_batch(fastq_batch_t *batch,
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index, 
				     fastq_batch_t *unmapped_batch,
				     array_list_t *mapping_list);

//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seeds_seq(char *seq, seed_t *seeds, unsigned int num_seeds,
					bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
					array_list_t *mapping_list);

unsigned int bwt_map_inexact_seed(char *seq,
                                     bwt_optarg_t *bwt_optarg,
                                     bwt_index_t *index,
                                     array_list_t *mapping_list);

unsigned int bwt_map_inexact_seq(char *seq, 
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index, 
				     array_list_t *mapping_list);

unsigned int bwt_map_inexact_read(fastq_read_t *read, 
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list);

unsigned int bwt_map_inexact_seqs(char **seqs, 
				      unsigned int num_reads,
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list);

unsigned int bwt_map_inexact_batch(fastq_batch_t *batch,
				       bwt_optarg_t *bwt_optarg, 
				       bwt_index_t *index, 
				       fastq_batch_t *unmapped_batch,
				       array_list_t *mapping_list);


//------------------------------------------------------------------------------

char *bwt_error_type(char error_kind);

//------------------------------------------------------------------------------
// Parameters for seed
//------------------------------------------------------------------------------

void seed_init(unsigned int starts, unsigned int end, seed_t *seed_p){
  seed_p->starts = starts;
  seed_p->end = end;
}

//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
// Paratemers for the candidate alignment localizations (CALs)
//------------------------------------------------------------------------------

cal_optarg_t *cal_optarg_new(const unsigned int min_cal_size, 
			     const unsigned int max_cal_distance, 
			     const unsigned int seed_size,
			     const unsigned int min_seed_size){
			       
  cal_optarg_t *cal_optarg_p = (cal_optarg_t *)malloc(sizeof(cal_optarg_t));
  cal_optarg_p->min_cal_size = min_cal_size;
  cal_optarg_p->max_cal_distance = max_cal_distance;
  cal_optarg_p->min_seed_size = min_seed_size;
  cal_optarg_p->seed_size = seed_size;
  
  return cal_optarg_p;
}
    
   
void cal_optarg_free(cal_optarg_t *optarg){
  free(optarg);
}

//------------------------------------------------------------------------------

cal_t *cal_new(const unsigned int chromosome_id, 
	       const short int strand,
	       const size_t start, 
	       const size_t end){
		 
  cal_t *cal_p = (cal_t *)malloc(sizeof(cal_t));  

  cal_p->chromosome_id = chromosome_id;
  cal_p->end = end;
  cal_p->start = start;
  cal_p->strand = strand;
  
  return cal_p;
}



void cal_free(cal_t *cal){
  free(cal);
}

//------------------------------------------------------------------------------

region_t *region_new(const unsigned int chromosome_id, 
	       const short int strand,
	       const size_t start, 
	       const size_t end){
		 
  region_t *region_p = (region_t *)malloc(sizeof(region_t));  
  region_p->chromosome_id = chromosome_id;
  region_p->end = end;
  region_p->start = start;
  region_p->strand = strand;
  
  return region_p;
}


void region_free(region_t *region_p){
  free(region_p);
}

//------------------------------------------------------------------------------

read_cals_t *read_cals_new(const fastq_read_t *read){
    read_cals_t *read_p = (read_cals_t *)malloc(sizeof(read_cals_t));
    read_p->cal_list = array_list_new(100000, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
    read_p->read = read;
}

void read_cals_free(read_cals_t *read_cals){
  fastq_read_free(read_cals->read);
  array_list_free(read_cals->cal_list, cal_free);
  free(read_cals);
}


//------------------------------------------------------------------------------
//Options settings Burrows-Wheeler Transform 
//------------------------------------------------------------------------------
bwt_optarg_t *bwt_optarg_new(const unsigned int num_errors,
			     const unsigned int num_threads,
			     const unsigned int max_alginments_per_read){
  
  bwt_optarg_t *bwt_optarg_p = (bwt_optarg_t *)malloc(sizeof(bwt_optarg_t));
  
  bwt_optarg_p->num_errors = num_errors;
  bwt_optarg_p->num_threads = num_threads;
  bwt_optarg_p->max_alignments_per_read = max_alginments_per_read;
  
  return bwt_optarg_p;
}

void bwt_optarg_free(bwt_optarg_t *optarg){
  free(optarg);
}
//-------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Index for the Burrows-Wheeler transform
//-----------------------------------------------------------------------------
char * bwt_error_type(char error_kind){
  char *msg_eror;
  switch(error_kind){
    case DELETION:
      msg_eror = "Deletion";
      break;
    case INSERTION:
      msg_eror = "Insertion";
      break;
    case MISMATCH:
      msg_eror = "Mismatch";
      break;
    default:
      msg_eror = "Unknown";
      break;
  }
  return msg_eror;
}

bwt_index_t *bwt_index_new(const char *dirname) {
  bwt_index_t *index = (bwt_index_t*) calloc(1, sizeof(bwt_index_t));

  index->dirname = strdup(dirname);

  readUIntVector(&index->h_C, dirname, "C");
  readUIntVector(&index->h_C1, dirname, "C1");
  readCompMatrix(&index->h_O, dirname, "O");
  readUIntCompVector(&index->S, dirname, "Scomp");

  reverseStrandC(&index->h_rC, &index->h_C,
  		 &index->h_rC1, &index->h_C1);

  reverseStrandO(&index->h_rO, &index->h_O);

  // 1-error handling 
  readCompMatrix(&index->h_Oi, dirname, "Oi");
  readUIntCompVector(&index->Si, dirname, "Scompi");

  reverseStrandO(&index->h_rOi, &index->h_Oi);

  // load karyotype
  char path[4096];
  sprintf(path, "%s/index.txt", dirname);
  load_exome_file(&index->karyotype, path);

  return index;
}

//-----------------------------------------------------------------------------

void bwt_index_free(bwt_index_t *index) {
  if (index != NULL) return;

  free(index->dirname);

  free(index->h_C.vector);
  free(index->h_C1.vector);
  free(index->h_rC.vector);
  free(index->h_rC1.vector);

  freeCompMatrix(&index->h_O);
  free(index->S.vector);

  // 1-error handling
  freeCompMatrix(&index->h_Oi);
  free(index->Si.vector);

  free(index);
}

//-----------------------------------------------------------------------------
// general functions
//-----------------------------------------------------------------------------

unsigned int bwt_map_seq(char *seq, 
			     bwt_optarg_t *bwt_optarg, 
			     bwt_index_t *index, 
			     array_list_t *mapping_list) {

  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_seq(seq, bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_seq(seq, bwt_optarg, index, mapping_list);  
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_read(fastq_read_t *read, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list) {

  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_read(read, bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_read(read, bwt_optarg, index, mapping_list);  
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_seqs(char **seqs, 
			      unsigned int num_reads,
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list) {

  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_seqs(seqs, num_reads, 
				  bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_seqs(seqs, num_reads, 
				  bwt_optarg, index, mapping_list);  
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_batch(fastq_batch_t *batch,
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       fastq_batch_t *unmapped_batch,
			       array_list_t *mapping_list) {

  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_batch(batch, bwt_optarg, 
				   index, unmapped_batch, mapping_list);
  }

  return bwt_map_inexact_batch(batch, bwt_optarg, 
				   index, unmapped_batch, mapping_list);  
}

//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_seq(char *seq, 
				   bwt_optarg_t *bwt_optarg, 
				   bwt_index_t *index, 
				   array_list_t *mapping_list) {
  
  
  
  unsigned int len = strlen(seq);
  int start = 0;
  int end = len - 1;
  char *code_seq = (char *) calloc(len, sizeof(char));
  result *result_p = (result *) calloc(1, sizeof(result));
  unsigned int l_aux, k_aux;
  unsigned int num_mappings = 0;
  char plusminus[2] = "-+";
  int idx, key, error, pos;
  alignment_t *alignment;
  char *cigar_p;
  unsigned int start_mapping;
  char *seq_dup;
  replaceBases(seq, code_seq, len);
  
  //printf("Search Read (%d): %s\n", len, seq);
  
  for (short int type = 1; type >= 0; type--) {
      if (type == 1) {
	result_p->k = 0;
	result_p->l = index->h_Oi.m_count - 2;
	BWExactSearchForward(code_seq, start, end, &index->h_C, &index->h_C1, &index->h_Oi, result_p);
      }else{
	result_p->k = 0;
	result_p->l = index->h_rO.m_count - 2;
	BWExactSearchForward(code_seq, start, end, &index->h_rC, &index->h_rC1, &index->h_rO, result_p);
      }
      
      k_aux = result_p->k;
      l_aux = result_p->l;
      
      if (l_aux - k_aux + 1 < bwt_optarg->max_alignments_per_read) {
	for (int j = k_aux; j <= l_aux; j++) {
	  if (index->S.ratio == 1) {
	    key = (type)
	      ? index->Si.n - index->Si.vector[j] - len - 1
	      : index->S.vector[j];
	  } else {
	    key = (type)
	      ? index->Si.n - getScompValue(j, &index->Si, &index->h_C,
					    &index->h_Oi) - len - 1
	      : getScompValue(j, &index->S, &index->h_C, &index->h_O);
	  }
	  idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
	  if(key + len <= index->karyotype.offset[idx]) {
	    start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	    /*printf("Strand:%c\tchromosome:%s\tStart:%u\n",plusminus[type],
		    index->karyotype.chromosome + (idx-1) * IDMAX,
		    start_mapping);
	    */
	    cigar_p = (char *)malloc(sizeof(char)*len);
	    sprintf(cigar_p, "%d=\0", len);
	    
	    seq_dup = (char *)malloc(sizeof(char)*len);
	    memcpy(seq_dup, seq, len);
	    // save all into one alignment structure and insert to the list
	    alignment = alignment_new();
	    alignment_init_single_end(NULL, seq, NULL, type, 
				      idx, //index->karyotype.chromosome + (idx-1) * IDMAX,
				      start_mapping, 
				      cigar_p, 1, 255, 1, (num_mappings > 0), alignment);

	    if(!array_list_insert((void*) alignment, mapping_list)){
         	printf("Error to insert item into array list\n");
	    }
	    
	    num_mappings++;
	  }
	}
      }
  }
  
  free(result_p);
  free(code_seq);
  return num_mappings;
  
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_read(fastq_read_t *read, 
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index, 
				    array_list_t *mapping_list) {
  
  unsigned int padding = array_list_size(mapping_list);
  
  unsigned int num_mappings = bwt_map_exact_seq(read->sequence, 
						    bwt_optarg, index, 
						    mapping_list);
  alignment_t *mapping;
  
  
  for (int i = 0; i < num_mappings ; i++) {
    //printf("access to pos:%d-padding:%d\n", i + padding, padding);
    mapping = (alignment_t *) array_list_get(i + padding, mapping_list);
    if(mapping != NULL){
      mapping->query_name = read->id;
      mapping->quality = read->quality;
    }else{
      printf("Error to extract item\n");
    }
    //printf("Store Ok!\n"); 
    //printf("\tSEQ:%s-ID:%s\n", mapping->sequence, mapping->query_name);
  }
  
  return num_mappings;
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_seqs(char **seqs, 
				    unsigned int num_reads,
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index, 
				    array_list_t *mapping_list) {
  
  unsigned int num_threads = bwt_optarg->num_threads;
  array_list_t **individual_mapping_list_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  const unsigned int chunk = MAX(1, num_reads/( num_threads*10));
  unsigned int num_mappings = 0; 
  unsigned int num_mappings_tot = 0;
  short int th_id;
  alignment_t *mapping;
  
  
  
  for(int i = 0; i < num_threads; i++){
   individual_mapping_list_p[i]  = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }
  
  omp_set_num_threads(num_threads);
  
  //printf("Initialization ok!Process %d reads x %d threads\n", num_reads, num_threads);
  
  #pragma omp parallel for private(th_id) schedule(dynamic, chunk) reduction(+:num_mappings_tot)
  for(int i = 0; i < num_reads; i++){
    th_id = omp_get_thread_num();
    //printf("I'm thread %d and process read %d\n", th_id,  i);
    num_mappings_tot += bwt_map_exact_seq(seqs[i], bwt_optarg, index, individual_mapping_list_p[th_id]);
  }
  
  //printf("Process reads ok! Total mappings %d\n\n", num_mappings_tot);
  
  //Join results
  for(int i = 0; i < num_threads; i++ ){
    num_mappings = array_list_size(individual_mapping_list_p[i]);
   
    for(int j = 0; j < num_mappings; j++){
      mapping = (alignment_t *)array_list_get(j, individual_mapping_list_p[i]);
      array_list_insert( mapping, mapping_list);
    }
  }

  free(individual_mapping_list_p);

  return num_mappings_tot;
}


//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_batch(fastq_batch_t *batch,
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index, 
				     fastq_batch_t *unmapped_batch,
				     array_list_t *mapping_list) {
  
  unsigned int length_header, length_seq;
  fastq_read_t *fq_read;
  char header[1024], seq[1024], quality[1024]; 
  unsigned int num_mappings, num_mappings_tot, num_unmapped;
  unsigned int read_pos;
  unsigned int num_threads = bwt_optarg->num_threads;
  unsigned int th_id;
  unsigned int num_reads = batch->num_reads;
  unsigned int batch_individual_size = batch->data_size / num_threads;
  const unsigned int chunk = MAX(1, num_reads/(num_threads*10));
  
  array_list_t **individual_mapping_list_p   = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  array_list_t **individual_unmapping_list_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  
  alignment_t *mapping;
  
  struct timeval start_time, end_time;
  double timer, parallel_t, sequential_t;
  
  //#pragma omp parallel for  
  for(int i = 0; i < num_threads; i++){
   individual_mapping_list_p[i]   = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
   individual_unmapping_list_p[i] = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }
  
  //printf("%d Threads %d chunk\n", num_threads, chunk);
  num_mappings = 0;
  omp_set_num_threads(num_threads);
  
  start_timer(start_time);
  #pragma omp parallel for private(th_id, fq_read, num_mappings) schedule(dynamic, chunk)
  for(int i = 0; i < num_reads; i++){
    th_id = omp_get_thread_num();
    
    fq_read = fastq_read_new(&(batch->header[batch->header_indices[i]]),  &(batch->seq[batch->data_indices[i]]), &(batch->quality[batch->data_indices[i]]));
    //fq_read = fastq_read_new(&header, &seq, &quality);
    //printf("Soy Thread %d and I process read %i : %s\n", th_id, i, fq_read->sequence);
    num_mappings = bwt_map_exact_read(fq_read, bwt_optarg, index, individual_mapping_list_p[th_id]);
    
    if(!num_mappings){
      //Read not mapped. Store in unmapped_batch
      //printf("Insert Error Read\n");
      if(!array_list_insert((void*)fq_read, individual_unmapping_list_p[th_id])){
         printf("Error to insert item into array list\n");
      }
    }
  }
  stop_timer(start_time, end_time, timer);
  parallel_t = timer;
  
  
  start_timer(start_time);
  //printf("Join Results\n");
  num_mappings_tot = 0;
  
  //Join results
  for(int i = 0; i < num_threads; i++ ){
    num_mappings = array_list_size(individual_mapping_list_p[i]);
    num_mappings_tot += num_mappings;
    for(int j = 0; j < num_mappings; j++){
      mapping = (alignment_t *)array_list_get(j, individual_mapping_list_p[i]);
      array_list_insert( mapping, mapping_list);
    }
  }

  //printf("Generate Batch\n");
  
  //Generate batch unmapped reads
  read_pos = 0;
  unmapped_batch->header_indices[read_pos] = 0;
  unmapped_batch->data_indices[read_pos] = 0;
  for(int i = 0; i < num_threads; i++){
    num_unmapped = array_list_size(individual_unmapping_list_p[i]);
    for(int j = 0; j < num_unmapped; j++){
      fq_read = (fastq_read_t *)array_list_get(j, individual_unmapping_list_p[i]);
      read_pos++;
      length_header = strlen(fq_read->id) + 1;
      length_seq = strlen(fq_read->sequence) + 1;
      memcpy(&(unmapped_batch->header[unmapped_batch->header_indices[read_pos - 1]]), fq_read->id,       length_header);
      memcpy(&(unmapped_batch->seq[unmapped_batch->data_indices[read_pos - 1]]),      fq_read->sequence, length_seq);
      memcpy(&(unmapped_batch->quality[unmapped_batch->data_indices[read_pos - 1]]),  fq_read->quality,  length_seq);
      
      unmapped_batch->data_indices[read_pos]   = unmapped_batch->data_indices[read_pos - 1]   + length_seq;
      unmapped_batch->header_indices[read_pos] = unmapped_batch->header_indices[read_pos - 1] + length_header;
    }
  }
  
  unmapped_batch->num_reads = read_pos;
  
  //printf("Free\n");
  
  for(int i = 0; i < num_threads; i++){
    array_list_free(individual_unmapping_list_p[i], fastq_read_free);
    //array_list_free(individual_mapping_list_p[i], alignment_free);
  }

  free(individual_mapping_list_p);
  free(individual_unmapping_list_p);
  
  stop_timer(start_time, end_time, timer);
  
  sequential_t = timer - parallel_t;
  
  global_parallel += parallel_t;
  global_sequential += sequential_t;
  
  //printf("Parallel %d\n", parallel_t);
  
  printf("%4.06f\t%4.06f\t%4.06f\n", parallel_t / 1000000, sequential_t / 1000000, timer / 1000000);
  
  return num_mappings_tot;
}

//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seed(char *seq, 
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list) {
    
  region_t *region_p; 
  unsigned int len = strlen(seq);
  int start = 0;
  int end = len - 1;
  size_t start_mapping;
  char *code_seq = (char *) calloc(len, sizeof(char));

  replaceBases(seq, code_seq, len);

  // calculate vectors k and l
  size_t *k0, *l0, *k1, *l1;
  size_t *ki0, *li0, *ki1, *li1;
  
  BWExactSearchBackwardVector(code_seq, start, end, 0, index->h_O.m_count - 2,
			      &k1, &l1, &index->h_C, &index->h_C1, &index->h_O);

  BWExactSearchForwardVector(code_seq, start, end, 0, index->h_Oi.m_count - 2,
			     &ki1, &li1, &index->h_C, &index->h_C1, &index->h_Oi);

  BWExactSearchForwardVector(code_seq, start, end, 0, index->h_rO.m_count - 2,
			     &k0, &l0, &index->h_rC, &index->h_rC1, &index->h_rO);

  BWExactSearchBackwardVector(code_seq, start, end, 0, index->h_rOi.m_count - 2,
			      &ki0, &li0, &index->h_rC, &index->h_rC1, &index->h_rOi);
  

  // compare the vectors k and l to get mappings in the genome
  unsigned int num_mappings = 0;
  char plusminus[2] = "-+";
  int idx, key, direction, error, pos;
  results_list *r_list;
  result *r;
  alignment_t *alignment;

  for (int type = 1; type >= 0; type--) {

    r_list = (results_list *) calloc(1, sizeof(results_list));
    new_results_list(r_list, 2000);
    r_list->n = 0;
    r_list->read_index = 0;

    if (type == 1) {
      BWIterativeSearch1(code_seq, start, end, k1, l1, ki1, li1, 
			 &index->h_C, &index->h_C1, &index->h_O, &index->h_Oi, r_list);
    } else {
      BWIterativeSearch1(code_seq, start, end, ki0, li0, k0, l0, 
     			 &index->h_rC, &index->h_rC1, &index->h_rOi, &index->h_rO, r_list);
    }

    for (size_t ii = 0; ii < r_list->n; ii++) {
      r = &r_list->list[ii];

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
	  
	  /*printf("%s\t%c\t%s %u %s error: %s, pos: %i, base: %i cigar: %s\n",
		 "nothing", plusminus[type],
		 index->karyotype.chromosome + (idx-1) * IDMAX,
		 index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]),
		 seq, bwt_error_type(r->err_kind[0]), r->position[0], r->base[0], cigar);*/
		  
	  start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	  // save all into one alignment structure and insert to the list
	  region_p = region_new(idx, type, start_mapping, start_mapping + end);
	  
	  if(!array_list_insert((void*) region_p, mapping_list)){
	    printf("Error to insert item into array list\n");
	  }
	  num_mappings++;
	}
      }
    }
    free_results_list(r_list);
  } // end for type 

  free(code_seq);
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

//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seq(char *seq, 
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index, 
				     array_list_t *mapping_list) {


  /*unsigned int len = strlen(seq);
  unsigned int start = 0;
  unsigned int end = len - 1;*/
  char *seq_dup;
  /*char *codeSeq = (char *) calloc(len, sizeof(char));
  replaceBases(seq, codeSeq, len);*/

  //return bwt_map_inexact_seq_by_pos(seq, start, end,
  //					bwt_optarg, index, mapping_list);
  
  unsigned int len = strlen(seq);
  int start = 0;
  int end = len - 1;

  char *code_seq = (char *) calloc(len, sizeof(char));

  replaceBases(seq, code_seq, len);

  // calculate vectors k and l
  size_t *k0, *l0, *k1, *l1;
  size_t *ki0, *li0, *ki1, *li1;
  
  BWExactSearchBackwardVector(code_seq, start, end, 0, index->h_O.m_count - 2,
			      &k1, &l1, &index->h_C, &index->h_C1, &index->h_O);

  BWExactSearchForwardVector(code_seq, start, end, 0, index->h_Oi.m_count - 2,
			     &ki1, &li1, &index->h_C, &index->h_C1, &index->h_Oi);

  BWExactSearchForwardVector(code_seq, start, end, 0, index->h_rO.m_count - 2,
			     &k0, &l0, &index->h_rC, &index->h_rC1, &index->h_rO);

  BWExactSearchBackwardVector(code_seq, start, end, 0, index->h_rOi.m_count - 2,
			      &ki0, &li0, &index->h_rC, &index->h_rC1, &index->h_rOi);
  

  // compare the vectors k and l to get mappings in the genome
  unsigned int num_mappings = 0;
  char plusminus[2] = "-+";
  int idx, key, direction, error, pos;
  results_list *r_list;
  result *r;
  alignment_t *alignment;
  
  unsigned int cigar_len;
  char *cigar_dup;
  char cigar[1024];
  unsigned int num_cigar_ops;

  for (int type = 1; type >= 0; type--) {

    r_list = (results_list *) calloc(1, sizeof(results_list));
    new_results_list(r_list, 2000);
    r_list->n = 0;
    r_list->read_index = 0;

    if (type == 1) {
      BWIterativeSearch1(code_seq, start, end, k1, l1, ki1, li1, 
			 &index->h_C, &index->h_C1, &index->h_O, &index->h_Oi, r_list);
    } else {
      BWIterativeSearch1(code_seq, start, end, ki0, li0, k0, l0, 
     			 &index->h_rC, &index->h_rC1, &index->h_rOi, &index->h_rO, r_list);
    }

    for (size_t ii = 0; ii < r_list->n; ii++) {
      r = &r_list->list[ii];
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
	  


	  // cigar processing
	  error = r->err_kind[0];
	  pos = r->position[0];

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

	  } else if (error == DELETION) {

	    if (pos == 0) {
	      sprintf(cigar, "1H%iM\0", len-1);
	      num_cigar_ops = 2;
	    } else if (pos == len - 1) {
	      sprintf(cigar, "%iM1H\0", len-1);
	      num_cigar_ops = 2;
	    } else {
	      sprintf(cigar, "%iM1D%iM\0", pos, len - pos - 1);
	      num_cigar_ops = 3;
	    }

	  } else if (error == INSERTION) {

	    if (pos == 0) {
	      sprintf(cigar, "1I%iM\0", len);
	      num_cigar_ops = 2;
	    } else if (pos == len - 1) {
	      sprintf(cigar, "%iM1I\0", len);
	      num_cigar_ops = 2;
	    } else {
	      sprintf(cigar, "%iM1I%iM\0", pos, len - pos);
	      num_cigar_ops = 3;
	    }

	  }
	  
	  /*printf("%s\t%c\t%s %u %s error: %s, pos: %i, base: %i cigar: %s\n",
		 "nothing", plusminus[type],
		 index->karyotype.chromosome + (idx-1) * IDMAX,
		 index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]),
		 seq, bwt_error_type(r->err_kind[0]), r->position[0], r->base[0], cigar);*/
		  
	  cigar_len = strlen(cigar);
	  cigar_dup = (char *)malloc(sizeof(char)*cigar_len);
	  memcpy(cigar_dup, &cigar, cigar_len + 1);
	  
	  seq_dup = (char *)malloc(sizeof(char)*len);
	  memcpy(seq_dup ,seq, len);
	  // save all into one alignment structure and insert to the list
	  alignment = alignment_new();
	  alignment_init_single_end(NULL, seq_dup, NULL, type, 
				    idx, //index->karyotype.chromosome + (idx-1) * IDMAX,
				    index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]), 
				    cigar_dup, num_cigar_ops, 255, 1, (num_mappings > 0), alignment);

	  array_list_insert((void*) alignment, mapping_list);

	  num_mappings++;
	}
      }
    }

    free_results_list(r_list);
  } // end for type 

  free(code_seq);
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

//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seqs_by_pos(char *seqs, 
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
	  
//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_read(fastq_read_t *read, 
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list) {
  
  unsigned int padding = array_list_size(mapping_list);
  
  unsigned int num_mappings = bwt_map_inexact_seq(read->sequence, 
						      bwt_optarg, index, 
						      mapping_list);
  //printf("padding:%d - num_mappings:%d\n", padding, num_mappings);
  alignment_t *mapping;
  for (int i = 0; i < num_mappings; i++) {
    mapping = (alignment_t *) array_list_get(i + padding, mapping_list);
    
    if(mapping != NULL){
      mapping->query_name = read->id;
      mapping->quality = read->quality;
    }else{
      printf("Error to extract item\n");
    }
  }
  
  return num_mappings;
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seqs(char **seqs, 
				      unsigned int num_reads,
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list) {
  
  unsigned int num_threads = bwt_optarg->num_threads;
  array_list_t **individual_mapping_list_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  const unsigned int chunk = MAX(1, num_reads/( num_threads*10));
  unsigned int num_mappings = 0; 
  unsigned int num_mappings_tot = 0;
  short int th_id;
  alignment_t *mapping;
  
  
  
  for(int i = 0; i < num_threads; i++){
   individual_mapping_list_p[i]  = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }
  
  omp_set_num_threads(num_threads);
  
  //printf("Initialization ok!Process %d reads x %d threads\n", num_reads, num_threads);
  
  #pragma omp parallel for private(th_id) schedule(dynamic, chunk) reduction(+:num_mappings_tot)
  for(int i = 0; i < num_reads; i++){
    th_id = omp_get_thread_num();
    //printf("I'm thread %d and process read %d\n", th_id,  i);
    num_mappings_tot += bwt_map_inexact_seq(seqs[i], bwt_optarg, index, individual_mapping_list_p[th_id]);
  }
  
  //printf("Process reads ok! Total mappings %d\n\n", num_mappings_tot);
  
  //Join results
  for(int i = 0; i < num_threads; i++ ){
    num_mappings = array_list_size(individual_mapping_list_p[i]);
   
    for(int j = 0; j < num_mappings; j++){
      mapping = (alignment_t *)array_list_get(j, individual_mapping_list_p[i]);
      array_list_insert( mapping, mapping_list);
    }
  }
  
  free(individual_mapping_list_p);
  
  return num_mappings_tot;
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_batch(fastq_batch_t *batch,
				       bwt_optarg_t *bwt_optarg, 
				       bwt_index_t *index, 
				       fastq_batch_t *unmapped_batch,
				       array_list_t *mapping_list) {
  
  unsigned int length_header, length_seq;
  fastq_read_t *fq_read;
  char header[1024], seq[1024], quality[1024]; 
  unsigned int num_mappings, num_mappings_tot, num_unmapped;
  unsigned int read_pos;
  unsigned int num_threads = bwt_optarg->num_threads;
  unsigned int th_id;
  unsigned int num_reads = batch->num_reads;
  unsigned int batch_individual_size = batch->data_size / num_threads;
  unsigned int chunk = MAX(1, num_reads/(num_threads*10));
  array_list_t **individual_mapping_list_p   = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  array_list_t **individual_unmapping_list_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  
  alignment_t *mapping;
  
  struct timeval start_time, end_time;
  double timer, parallel_t, sequential_t;
  
  for(int i = 0; i < num_threads; i++){
   individual_mapping_list_p[i]   = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
   individual_unmapping_list_p[i] = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }  
  
  //printf("%d Threads %d chunk\n", num_threads, chunk);
  num_mappings = 0;
  omp_set_num_threads(num_threads);
  
  start_timer(start_time);
  #pragma omp parallel for private(th_id, fq_read, num_mappings) schedule(dynamic, chunk)
  for(int i = 0; i < num_reads; i++){
    th_id = omp_get_thread_num();
    //th_id = 0;
    //printf("I'm Thread %d of %d\n", th_id, omp_get_thread_num());
    fq_read = fastq_read_new(&(batch->header[batch->header_indices[i]]),  &(batch->seq[batch->data_indices[i]]), &(batch->quality[batch->data_indices[i]]));
    num_mappings = bwt_map_inexact_read(fq_read, bwt_optarg, index, individual_mapping_list_p[th_id]);
    //num_mappings = 0;
    if(!num_mappings){
      //Read not mapped. Store in unmapped_batch
      //printf("\tNo mapped \n");
      array_list_insert((void*)fq_read, individual_unmapping_list_p[th_id]);
    }else{
      //printf("\tMapped \n");
    }
  }
  stop_timer(start_time, end_time, timer);
  parallel_t = timer;
  
  start_timer(start_time);
  //printf("Join Results\n");
  num_mappings_tot = 0;
  //Join results
  for(int i = 0; i < num_threads; i++ ){
    num_mappings = array_list_size(individual_mapping_list_p[i]);
    num_mappings_tot += num_mappings;
    for(int j = 0; j < num_mappings; j++){
      mapping = (alignment_t *)array_list_get(j, individual_mapping_list_p[i]);
      array_list_insert( mapping, mapping_list);
    }
  }
  //printf("Generate Batch\n");
  //Generate batch unmapped reads
  read_pos = 0;
  unmapped_batch->header_indices[read_pos] = 0;
  unmapped_batch->data_indices[read_pos] = 0;
  for(int i = 0; i < num_threads; i++){
    num_unmapped = array_list_size(individual_unmapping_list_p[i]);
    for(int j = 0; j < num_unmapped; j++){
      fq_read = (fastq_read_t *)array_list_get(j, individual_unmapping_list_p[i]);
      read_pos++;
      length_header = strlen(fq_read->id) + 1;
      length_seq = strlen(fq_read->sequence) + 1;
      memcpy(&(unmapped_batch->header[unmapped_batch->header_indices[read_pos - 1]]), fq_read->id,       length_header);
      memcpy(&(unmapped_batch->seq[unmapped_batch->data_indices[read_pos - 1]]),      fq_read->sequence, length_seq);
      memcpy(&(unmapped_batch->quality[unmapped_batch->data_indices[read_pos - 1]]),  fq_read->quality,  length_seq);
      
      unmapped_batch->data_indices[read_pos]   = unmapped_batch->data_indices[read_pos - 1]   + length_seq;
      unmapped_batch->header_indices[read_pos] = unmapped_batch->header_indices[read_pos - 1] + length_header;
    }
  }
  
  unmapped_batch->num_reads = read_pos;
  //printf("Free\n");
  
  for(int i = 0; i < num_threads; i++){
    array_list_free(individual_unmapping_list_p[i], fastq_read_free);
    //array_list_free(individual_mapping_list_p[i], alignment_free);
  }

  free(individual_mapping_list_p);
  free(individual_unmapping_list_p);
  
  stop_timer(start_time, end_time, timer);
  
  sequential_t = timer - parallel_t;
  
  global_parallel += parallel_t;
  global_sequential += sequential_t;
  
  //printf("Parallel %d\n", parallel_t);
  
  printf("%4.06f\t%4.06f\t%4.06f\n", parallel_t / 1000000, sequential_t / 1000000, timer / 1000000);
  return num_mappings_tot;

}

//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seeds_seq2(char *seq, bwt_optarg_t *bwt_optarg, 
					bwt_index_t *index, array_list_t *mapping_list) {

  size_t len = strlen(seq);
  size_t offset, num_seeds = len / seed_size;

  char *code_seq = (char *) calloc(len, sizeof(char));

  replaceBases(seq, code_seq, len);

  // first 'pasada'
  offset = 0;
  for (size_t i = 0; i < num_seeds; i++) {
    bwt_map_inexact_seed(code_seq, offset, offset + seed_size - 1,
			 bwt_optarg, index, mapping_list);
    offset += seed_size;
  }

  // special processing for the last seed !!
  if (len % seed_size >= bwt_optarg->min_seed_size) {
    bwt_map_inexact_seed(code_seq, offset, len - 1,
			 bwt_optarg, index, mapping_list);
  }

  // second 'pasada', shifting seeds by (seed_size / 2)
  offset = seed_size / 2;
  num_seeds = (len - seed_size / 2) / seed_size;
  for (size_t i = 0; i < num_seeds; i++) {
    bwt_map_inexact_seed(code_seq, offset, offset + seed_size - 1,
			 bwt_optarg, index, mapping_list);
    offset += seed_size;
  }

  // again, special processing for the last seed !!
  if ((len - seed_size / 2) % seed_size >= bwt_optarg->min_seed_size) {
    bwt_map_inexact_seed(code_seq, offset, len - 1,
			 bwt_optarg, index, mapping_list);
  }

  return array_list_size(mapping_list);
  
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seeds_seq(char *seq, seed_t *seeds, unsigned int num_seeds,
					bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
					array_list_t *mapping_list){
  
  printf("Process Read(%d): %s\n", strlen(seq), seq);
  char *seq_seed = (char *) calloc(sizeof(char), strlen(seq));
  
  unsigned int num_mappings = 0;
  unsigned int mappings_tot = 0;
  unsigned int len, end;
  
  for(int i = 0; i < num_seeds; i++){
    printf("Process read_seed([%d-%d])\n", seeds[i].starts, seeds[i].end);
    len = seeds[i].end - seeds[i].starts + 1;
    memcpy(seq_seed, seq + seeds[i].starts, len);
    seq_seed[len] = '\0';
    printf("\tSeed: %s (len = %d)\n",   seq_seed, len);
    num_mappings = bwt_map_inexact_seed(seq_seed, bwt_optarg, index, mapping_list);
    mappings_tot += num_mappings;
    printf("\tMappings : %d (total = %d, mapping list size = %d)\n", num_mappings, mappings_tot, 
	   array_list_size(mapping_list));
  }
  
  printf("End Mapping seeds\n");
  return mappings_tot;
  
}

unsigned int bwt_map_inexact_seeds_read(fastq_read_t *read, seed_t *seeds, unsigned int num_seeds,
					    bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
					    array_list_t *mapping_list){
  
	return bwt_map_inexact_seeds_seq(read->sequence, seeds, num_seeds, bwt_optarg, index, mapping_list);				  
}

//-----------------------------------------------------------------------------
// Seed functions
//-----------------------------------------------------------------------------

seed_t *create_seeds_from_seq(char *seq, 
				  unsigned int seed_size,
				  unsigned int min_seed_size,
				  unsigned int *num_seeds) {

  unsigned int end_acum = 0;
  unsigned int nb_seeds = 0;
  unsigned int len = strlen(seq);
  double res = ceil((double)len / (double)seed_size);
  unsigned int max_num_seeds = (unsigned int)res * 2;
  
  
  printf("seq(%d):%s \n", len, seq);
  printf("\t seed_size = %d :: len = %d :: max_num_seeds = %d :: min_seed_size = %d\n", seed_size, len,  max_num_seeds, min_seed_size);
  
  seed_t *seed, *seeds = (seed_t*) calloc(max_num_seeds, sizeof(seed_t));
    
  // generate seeds by splitting read into seeds
  // first set of seeds starting at position 0
  for (unsigned int i = 0; i < len; i += seed_size) {
    seed_init(i, i + seed_size - 1, &seeds[nb_seeds++]);
   // printf("1º-[%d-%d]\n", i, i + seed_size);
  }
  
  if (nb_seeds > 0) {
    seed = &seeds[nb_seeds - 1];
    if ( seed->end >= len) {
      if (len - seed->starts >= min_seed_size) {
	seed->end = len - 1;
      } else {
	nb_seeds--;
      }
    }
  }
  
  // second set of seeds starting at position (seed_size / 2)
  for (unsigned int i = (seed_size / 2); i < len; i += seed_size) {
    seed_init(i, i + seed_size - 1, &seeds[nb_seeds++]);
    //printf("2º-[%d-%d]\n", i, i + seed_size);
  }
  
  if (nb_seeds > 0) {
    seed = &seeds[nb_seeds - 1];
    if ( seed->end >= len) {
      if ((len - seed->starts) >= min_seed_size) {
	seed->end = len - 1;
      } else {
	nb_seeds--;
      }
    }
  }

  *num_seeds = nb_seeds;
  
  return seeds;
}

//-----------------------------------------------------------------------------

seed_t *create_seeds_from_read(fastq_read_t *read, 
				   unsigned int seed_size,
				   unsigned int min_seed_size,
				   unsigned int *num_seeds) {
  
    return create_seeds_from_seq(read->sequence,  
				     seed_size,
				     min_seed_size,
				     num_seeds);
}

//-----------------------------------------------------------------------------

seed_t **create_seeds_from_seqs(char **seqs, 
				    unsigned int seed_size,
				    unsigned int min_seed_size,
				    unsigned int *num_seeds,
				    unsigned int num_reads,
				    bwt_optarg_t *bwt_optarg){
  
    unsigned int num_threads = bwt_optarg->num_threads;
    const unsigned int chunk = MAX(1, num_reads / (num_threads * 10));
    seed_t **seeds = (seed_t **) calloc(num_reads, sizeof(seed_t *));
    
    omp_set_num_threads(num_threads);
    
    //printf("Initialization ok!Process %d reads x %d threads\n", num_reads, num_threads);
    
    #pragma omp parallel for private(num_seeds) schedule(dynamic, chunk)
    for (int i = 0; i < num_reads; i++) {
      seeds[i] = create_seeds_from_seq(seqs[i], 
					   seed_size,
					   min_seed_size,
					   &num_seeds[i]);
    }
    
    return seeds;
}

//-----------------------------------------------------------------------------

//TODO: Generate batch of seeds? complete implementation, seed_batch_t??
seed_t **create_seeds_from_batch(fastq_batch_t *batch, 
				     unsigned int seed_size,
				     unsigned int min_seed_size,
				     unsigned int *num_seeds,
				     bwt_optarg_t *bwt_optarg){
  
    unsigned int num_threads = bwt_optarg->num_threads;
    unsigned int num_reads = batch->num_reads;
    const unsigned int chunk = MAX(1, num_reads / (num_threads * 10));
    seed_t **seeds = (seed_t **) calloc(num_reads, sizeof(seed_t *));
    
    //fastq_read_t fq_read;
    omp_set_num_threads(num_threads);
    
    //printf("Initialization ok!Process %d reads x %d threads\n", num_reads, num_threads);
    
    #pragma omp parallel for private(num_seeds) schedule(dynamic, chunk)
    for(int i = 0; i < num_reads; i++){
       //fq_read = fastq_read_new(&(batch->header[batch->header_indices[i]]),  &(batch->seq[batch->data_indices[i]]), &(batch->quality[batch->data_indices[i]]));
       seeds[i] = create_seeds_from_seq(&(batch->seq[batch->data_indices[i]]), 
					    seed_size,
					    min_seed_size,
					    &num_seeds[i]);
       if(num_seeds[i] != 0) {
          //Insert to split reads batch
       }
    }
    
    return seeds;
  
}

//-----------------------------------------------------------------------------
// CAL functions
//-----------------------------------------------------------------------------

unsigned int bwt_generate_cal_list(array_list_t *mapping_list,
				       cal_optarg_t *cal_optarg,
				       array_list_t *cal_list) {
  cal_t *cal, *cal1;
  region_t *region_p;
  unsigned int extend, added, cal_idx;
  unsigned int min_cal_size = cal_optarg->min_cal_size;
  unsigned int max_cal_distance  = cal_optarg->max_cal_distance;
  unsigned int num_mappings = array_list_size(mapping_list);
  unsigned int num_cals;
  unsigned int chromosome_id;
  short int strand;
  size_t start, end;
  
  int print = 0;
  
  printf("Start generate CALs: %d mappings\n", num_mappings);
  // from the results mappings, generates CALs
  for (unsigned int m = 0; m < num_mappings; m++) {
    region_p = array_list_get(m, mapping_list);

    chromosome_id = region_p->chromosome_id;
    strand = region_p->strand;
    start = region_p->start;
    end = region_p->end;

      if ( (strand == 1 && chromosome_id == 1 && (start >= 21010025 && start <= 21210025)) ||
	   (strand == 1 && chromosome_id == 1 && (end >= 21010078 && end <= 21210078)) ) {
	printf("----> mapping, start = %d, end = %d\n", start, end);
      }

    extend = 0;
    added = 0;
    num_cals = array_list_size(cal_list);
    //printf("Process region strand %d chromosome %d [%d-%d]\n", strand, chromosome_id, start, end);
    for (unsigned int c = 0; c < num_cals; c++) {
      cal = array_list_get(c, cal_list);

      print = 0;
      if ( (strand == 1 && cal->chromosome_id == 1 && (cal->start >= 21010025 && cal->start <= 21210025)) ||
	   (strand == 1 && cal->chromosome_id == 1 && (cal->end >= 21010078 && cal->end <= 21210078)) ) {
	print = 1;
      }

      // locate the alignment among the CALs
      if (chromosome_id == cal->chromosome_id && 
	  strand == cal->strand) {
	// chromosome and strand match
	if (end + max_cal_distance < cal->start ||
	    cal->end + max_cal_distance < start) {
	  // too far, next CAL...
	  continue;
	} else {
	  added = 1;
	  if (start >= cal->start && end <= cal->end) {
	    // the aligmnent is INTO the CAL
	    break;
	  } else {
	    cal_idx = c;
	    extend = 1;
	    
	    if (print) 
	      printf("mapping: start = %d end = %d\n", start, end);

	    if (start < cal->start) {
	      if (print) 
		printf("\tChange start: chrm %d CAL %d [%d-%d] -> new start = %d\n", cal->chromosome_id, c, cal->start, cal->end, start);
	      cal->start = start; 
	    }
	    if (end > cal->end) {
	      if (print) 
		printf("\tChange end: chrm %d CAL %d [%d-%d] -> new end = %d\n", cal->chromosome_id, c, cal->start, cal->end, end);
	      cal->end = end; 
	    }
	    if (print) 
	      printf("\t\tExtend CAL %d Result :: %d-%d\n", c, cal->start, cal->end);
	    break;
	  }
	}
      }
    } // end for cal_list->size

    if (added) {
      if (extend && num_cals > 1) {
	//printf("Extend\n");
	// are there CALS overlapping after extending a CAL?
	start = cal->start;
	end = cal->end;
	chromosome_id = cal->chromosome_id;
	strand = cal->strand;
	extend = 0;
	num_cals = array_list_size(cal_list);
	for (unsigned int c = 0; c < num_cals; c++) {
	  if (c != cal_idx) {
	    cal1 = array_list_get(c, cal_list);
	    if((cal1->chromosome_id == chromosome_id) && (cal->strand == strand) ){
	      //printf("\t1.%d + %d = %d >= %d ??\n", end, max_cal_distance, (size_t)(end + max_cal_distance), cal1->start);
	      if(end + max_cal_distance < cal1->start ||
		 cal1->end + max_cal_distance < start) {
		  // too far, next CAL...
		  continue;
	      }
	      else{
		if (print) 
		  printf("\tFound CAL in range. Start End Actualization cal:[%d-%d] extend cal1:[%d-%d]\n", cal->start, cal->end, cal1->start, cal1->end);
		if(start > cal1->start){
		  cal->start = cal1->start;
		  extend = 1;
		}
		if(end < cal1->end){
		  cal->end = cal1->end;
		  extend = 1;
		}
		if (print) 
		  printf("\t\tResult cal [%d-%d] and Delete %d\n", cal->start, cal->end, extend);
	      }
	    }
	    if (extend) {
	      array_list_remove_at(c, cal_list);
	      break;
	    }
	  }
	}
      }
    } else {
      array_list_insert(cal_new(chromosome_id, strand, start, end),
			cal_list);
    }
  }
  printf("End\n");
  
  num_cals = array_list_size(cal_list);
  // removing 'small' CALs from list
  for (int i = num_cals - 1; i >= 0 ; i--) {
    cal = array_list_get(i, cal_list);
    //printf("(%d - %d + 1= %d) < %d \n", cal->end, cal->start, cal->end - cal->start + 1, min_cal_size);
    if (cal->end - cal->start + 1 < min_cal_size) {
      array_list_remove_at(i, cal_list);
    } else{
      printf("Not Remove strand %d chromosome %d CAL %d-%d\n", cal->strand, cal->chromosome_id, cal->start, cal->end);
    }
  }
  
  return array_list_size(cal_list);
}

//

//-----------------------------------------------------------------------------

unsigned int bwt_find_cals_from_seq(char *seq, 
					bwt_optarg_t *bwt_optarg, 
					bwt_index_t *index, 
					cal_optarg_t *cal_optarg, 
					array_list_t *cal_list) {
  // generate the seeds
  unsigned int num_seeds = 0;
  seed_t *seeds = create_seeds_from_seq(seq,
					cal_optarg->seed_size,
					cal_optarg->min_seed_size,
					&num_seeds);
  
  for(int i = 0; i < num_seeds; i++){
    printf("Seed %d: [%d-%d]\n", i, seeds[i].starts, seeds[i].end);
  }

  // create list where to store the seed mappings
  array_list_t *mapping_list = array_list_new(100000, 
					      1.25f, 
					      COLLECTION_MODE_SYNCHRONIZED); 
 printf("Array list complete\n");

  // map seeds
  unsigned int num_mappings = bwt_map_inexact_seeds_seq(seq, seeds, num_seeds,
							bwt_optarg, index, 
							mapping_list);
  

  printf("num_mappings = %d, mapping_list_size = %d\n", num_mappings, array_list_size(mapping_list));

  // generate cals from mapped seeds
  unsigned int num_cals = bwt_generate_cal_list(mapping_list,
						cal_optarg, cal_list);
  // free memory for mapping list
  array_list_free(mapping_list, NULL);

  return num_cals;
}

//-----------------------------------------------------------------------------

unsigned int bwt_find_cals_from_read(fastq_read_t *read, 
					 bwt_optarg_t *bwt_optarg, 
					 bwt_index_t *index, 
					 cal_optarg_t *cal_optarg, 
					 array_list_t *cal_list) {

  return bwt_find_cals_from_seq(read->sequence, 
				    bwt_optarg, index,
				    cal_optarg, cal_list);
}

//-----------------------------------------------------------------------------

unsigned int bwt_find_cals_from_seqs(char **seqs, 
					 unsigned int num_reads,
					 bwt_optarg_t *bwt_optarg, 
					 bwt_index_t *index, 
					 unsigned int *num_unmapped,
					 int *unmapped_array,
					 cal_optarg_t *cal_optarg, 
					 array_list_t *cal_list) {
  
    unsigned int num_seeds = 0;
    unsigned int num_threads = bwt_optarg->num_threads;
    const unsigned int chunk = MAX(1, num_reads/( num_threads*10));
    seed_t **seeds_p = (seed_t **) calloc(num_reads, sizeof(seed_t *));
    unsigned int num_mappings;
    fastq_read_t fq_read;
    seed_t *seeds;
    unsigned int th_id;
    unsigned int num_cals;
    unsigned int total_cals;
    cal_t *cal;
    // create list where to store the seed mappings
    array_list_t **mapping_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
    array_list_t **cal_individual_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
    
    for(int i = 0; i < num_threads; i++){
	mapping_list[i] = array_list_new(100000, 
		                         1.25f, 
		                         COLLECTION_MODE_SYNCHRONIZED);
	
        cal_individual_list[i] = array_list_new(100000, 
		                         1.25f, 
		                         COLLECTION_MODE_SYNCHRONIZED);
    }

    omp_set_num_threads(num_threads);
    
    //printf("Initialization ok!Process %d reads x %d threads\n", num_reads, num_threads);
    num_unmapped = 0;
    #pragma omp parallel for firstprivate(num_seeds) private(th_id, seeds, num_mappings, num_cals) schedule(dynamic, chunk)
    for(int i = 0; i < num_reads; i++){
       th_id = omp_get_thread_num();
       // generate the seeds
       seeds = create_seeds_from_seq(seqs[i],
				         cal_optarg->seed_size,
				         cal_optarg->min_seed_size,
				         &num_seeds);
       /*for(int i = 0; i < num_seeds; i++){
          printf("Seed %d: [%d-%d]\n", i, seeds[i].starts, seeds[i].end);
         }*/

        //printf("Array list complete\n");
       // map seeds
       num_mappings = bwt_map_inexact_seeds_seq(seqs[i], seeds, num_seeds,
						    bwt_optarg, index, 
						    mapping_list[th_id]);

       // generate cals from mapped seeds
       num_cals = bwt_generate_cal_list(mapping_list[th_id],
                                            cal_optarg, cal_individual_list[th_id]);
      

       if(num_seeds == 0 || num_mappings == 0 || num_cals == 0){
		unmapped_array[*num_unmapped] = i;
		*num_unmapped += 1;
       }
      // free memory for mapping list
      array_list_free(mapping_list, NULL);

    }
    
    for(int i = 0; i < num_threads; i++){
        num_cals = array_list_size(cal_individual_list[i]);
        total_cals += num_cals;
        array_list_free(mapping_list[i], region_free);
	for(int j = 0; j < num_cals; j++){
    	    cal = (cal_t *)array_list_get(j, cal_individual_list[i]);
            array_list_insert(cal, cal_list);
	}
    }

    return total_cals;
}

//-----------------------------------------------------------------------------

unsigned int bwt_find_cals_from_batch(fastq_batch_t *batch,
					  bwt_optarg_t *bwt_optarg, 
					  bwt_index_t *index, 
					  fastq_batch_t *unmapped_batch,
					  cal_optarg_t *cal_optarg, 
					  array_list_t *cal_list) {
  return 0;
}

//-----------------------------------------------------------------------------
