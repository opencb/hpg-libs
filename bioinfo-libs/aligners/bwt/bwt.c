#include "bwt.h"


//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------

size_t bwt_map_exact_seq(char *seq, 
			 bwt_optarg_t *bwt_optarg, 
			 bwt_index_t *index, 
			 array_list_t *mapping_list);

size_t bwt_map_exact_read(fastq_read_t *read, 
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index, 
			  array_list_t *mapping_list);

size_t bwt_map_exact_seqs(char **seqs, 
			  size_t num_reads,
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index, 
			  char *out_status,
			  array_list_t *mapping_list);

size_t bwt_map_exact_batch(fastq_batch_t *batch,
			   bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, 
			   fastq_batch_t *unmapped_batch,
			   array_list_t *mapping_list);

size_t bwt_map_exact_seeds_seq(char *seq, size_t seed_size, size_t min_seed_size,
				 bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				 array_list_t *mapping_list);
//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------
size_t bwt_map_exact_seeds_seq(char *seq, size_t seed_size, size_t min_seed_size,
				 bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				 array_list_t *mapping_list);


size_t bwt_map_inexact_seeds_seq(char *seq, size_t seed_size, size_t min_seed_size,
				 bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				 array_list_t *mapping_list);

size_t bwt_map_inexact_seed(char *seq, size_t start, size_t end,
			    bwt_optarg_t *bwt_optarg,
			    bwt_index_t *index,
			    array_list_t *mapping_list);

//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seq(char *seq, 
			   bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, 
			   array_list_t *mapping_list);

size_t bwt_map_inexact_read(fastq_read_t *read, 
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    array_list_t *mapping_list);

size_t bwt_map_inexact_seqs(char **seqs, 
			    size_t num_reads,
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    char *out_status,
			    array_list_t *mapping_list);

size_t bwt_map_inexact_batch(fastq_batch_t *batch,
			     bwt_optarg_t *bwt_optarg, 
			     bwt_index_t *index, 
			     fastq_batch_t *unmapped_batch,
			     array_list_t *mapping_list);


//------------------------------------------------------------------------------

char *bwt_error_type(char error_kind);

//------------------------------------------------------------------------------
// Paratemers for the candidate alignment localizations (CALs)
//------------------------------------------------------------------------------

cal_optarg_t *cal_optarg_new(const size_t min_cal_size, 
			     const size_t max_cal_distance, 
			     const size_t seed_size,
			     const size_t min_seed_size,
			     const size_t num_errors){
			       
  cal_optarg_t *cal_optarg_p = (cal_optarg_t *)malloc(sizeof(cal_optarg_t));
  cal_optarg_p->min_cal_size = min_cal_size;
  cal_optarg_p->max_cal_distance = max_cal_distance;
  cal_optarg_p->min_seed_size = min_seed_size;
  cal_optarg_p->seed_size = seed_size;
  cal_optarg_p->num_errors = num_errors;

  return cal_optarg_p;
}
    
void cal_optarg_free(cal_optarg_t *optarg){
  free(optarg);
}

//------------------------------------------------------------------------------

cal_t *cal_new(const size_t chromosome_id, 
	       const short int strand,
	       const size_t start, 
	       const size_t end) {
		 
  cal_t *cal = (cal_t *)malloc(sizeof(cal_t));  

  cal->chromosome_id = chromosome_id;
  cal->end = end;
  cal->start = start;
  cal->strand = strand;
  
  return cal;
}


void cal_free(cal_t *cal){
  free(cal);
}

//------------------------------------------------------------------------------

short_cal_t *short_cal_new(const size_t start, 
	                   const size_t end){	 
  short_cal_t *short_cal = (short_cal_t *)malloc(sizeof(short_cal_t));  
  
  short_cal->end = end;
  short_cal->start = start;
  
  return short_cal;
}


void short_cal_free(short_cal_t *short_cal){
  free(short_cal);
}

//------------------------------------------------------------------------------

region_t *region_new(const size_t chromosome_id, 
		     const short int strand,
		     const size_t start, 
		     const size_t end){
		 
  region_t *region = (region_t *) malloc(sizeof(region_t));  

  region->chromosome_id = chromosome_id;
  region->end = end;
  region->start = start;
  region->strand = strand;
  
  return region;
}


void region_free(region_t *region){
  free(region);
}

//------------------------------------------------------------------------------

read_cals_t *read_cals_new(const fastq_read_t *read) {

    read_cals_t *read_cals = (read_cals_t *) malloc(sizeof(read_cals_t));

    read_cals->cal_list = array_list_new(100000, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
    read_cals->read = read;

    return read_cals;
}

void read_cals_free(read_cals_t *read_cals){
  fastq_read_free(read_cals->read);
  array_list_free(read_cals->cal_list, cal_free);
  free(read_cals);
}

//------------------------------------------------------------------------------
//Options settings Burrows-Wheeler Transform 
//------------------------------------------------------------------------------

bwt_optarg_t *bwt_optarg_new(const size_t num_errors,
			     const size_t num_threads,
			     const size_t max_alginments_per_read) {
  
  bwt_optarg_t *bwt_optarg = (bwt_optarg_t *) calloc(1, sizeof(bwt_optarg_t));
  
  bwt_optarg->num_errors = num_errors;
  bwt_optarg->num_threads = num_threads;
  bwt_optarg->max_alignments_per_read = max_alginments_per_read;
  
  return bwt_optarg;
}

void bwt_optarg_free(bwt_optarg_t *optarg) {
  if (optarg != NULL) {
    free(optarg);
  }
}

//-----------------------------------------------------------------------------
// Index for the Burrows-Wheeler transform
//-----------------------------------------------------------------------------

char *bwt_error_type(char error_kind){

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

//-----------------------------------------------------------------------------

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
  char path[strlen(dirname) + 512];
  //  sprintf(path, "%s/index.txt", dirname);
  //load_exome_file(&index->karyotype, path);

  load_exome_file(&index->karyotype, dirname);

  return index;
}

//-----------------------------------------------------------------------------

void bwt_index_free(bwt_index_t *index) {
  if (index == NULL) return;

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

void bwt_generate_index_files(char *ref_file, char *output_dir, 
			      unsigned int s_ratio) {

  byte_vector X, B, Bi;
  vector C, C1;
  comp_vector S, Si, Scomp, Scompi;
  comp_vector R, Ri, Rcomp, Rcompi;
  comp_matrix O, Oi;

  exome ex;

  initReplaceTable();

  // Calculating BWT
  calculateBWT(&B, &S, &X, 0, &ex, ref_file);

  save_exome_file(&ex, output_dir);

  saveCharVector(&X, output_dir, "X");
  free(X.vector);

  printUIntVector(S.vector, S.n);
  printUIntVector(B.vector, B.n);

  // Calculating prefix-trie matrices C and O
  calculateC(&C, &C1, &B, 0);
  calculateO(&O, &B);

  printUIntVector(C.vector, C.n);
  printUIntVector(C1.vector, C1.n);
  printCompMatrix(O);

  saveCharVector(&B, output_dir, "B");
  free(B.vector);
  saveUIntVector(&C, output_dir, "C");
  free(C.vector);
  saveUIntVector(&C1, output_dir, "C1");
  free(C1.vector);
  saveCompMatrix(&O, output_dir, "O");
  freeCompMatrix(&O);

  // Calculating R
  calculateR(&S, &R);

  printUIntVector(R.vector, R.n);

  // Calculating Scomp Rcomp
  calculateSRcomp(&S, &Scomp, s_ratio);
  printUIntVector(Scomp.vector, Scomp.n);
  calculateSRcomp(&R, &Rcomp, s_ratio);
  printUIntVector(Rcomp.vector, Rcomp.n);


  saveUIntCompVector(&S, output_dir, "S");
  free(S.vector);
  saveUIntCompVector(&R, output_dir, "R");
  free(R.vector);
  saveUIntCompVector(&Scomp, output_dir, "Scomp");
  free(Scomp.vector);
  saveUIntCompVector(&Rcomp, output_dir, "Rcomp");
  free(Rcomp.vector);

  //Calculating BWT of reverse reference
  calculateBWT(&Bi, &Si, &X, 1, NULL, ref_file);

  saveCharVector(&X, output_dir, "Xi");
  free(X.vector);

  printUIntVector(Bi.vector, Bi.n);
  printUIntVector(Si.vector, Si.n);

  //Calculating inverted prefix-trie matrix Oi
  calculateO(&Oi, &Bi);

  printCompMatrix(Oi);

  saveCharVector(&Bi, output_dir, "Bi");
  free(Bi.vector);

  saveCompMatrix(&Oi, output_dir, "Oi");
  freeCompMatrix(&Oi);

  //Calculating Ri
  calculateR(&Si, &Ri);

  printUIntVector(Ri.vector, Ri.n);

  // Calculating Scompi Rcompi
  calculateSRcomp(&Si, &Scompi, s_ratio);
  printUIntVector(Scompi.vector, Scompi.n);
  calculateSRcomp(&Ri, &Rcompi, s_ratio);
  printUIntVector(Rcompi.vector, Rcompi.n);

  saveUIntCompVector(&Si, output_dir, "Si");
  free(Si.vector);
  saveUIntCompVector(&Ri, output_dir, "Ri");
  free(Ri.vector);
  saveUIntCompVector(&Scompi, output_dir, "Scompi");
  free(Scompi.vector);
  saveUIntCompVector(&Rcompi, output_dir, "Rcompi");
  free(Rcompi.vector);
}


//-----------------------------------------------------------------------------
// general functions
//-----------------------------------------------------------------------------

size_t bwt_map_seq(char *seq, 
		   bwt_optarg_t *bwt_optarg, 
		   bwt_index_t *index, 
		   array_list_t *mapping_list) {
  
  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_seq(seq, bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_seq(seq, bwt_optarg, index, mapping_list);  
}

//-----------------------------------------------------------------------------

size_t bwt_map_read(fastq_read_t *read, 
		    bwt_optarg_t *bwt_optarg, 
		    bwt_index_t *index, 
		    array_list_t *mapping_list) {
  
  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_read(read, bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_read(read, bwt_optarg, index, mapping_list);  
}

//-----------------------------------------------------------------------------

size_t bwt_map_seqs(char **seqs, 
		    size_t num_reads,
		    bwt_optarg_t *bwt_optarg, 
		    bwt_index_t *index, 
		    char *out_status,
		    array_list_t *mapping_list) {

  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_seqs(seqs, num_reads, bwt_optarg, index, 
			      out_status, mapping_list);
  }
  
  return bwt_map_inexact_seqs(seqs, num_reads, bwt_optarg, index, 
			      out_status, mapping_list);  
}

//-----------------------------------------------------------------------------

size_t bwt_map_batch(fastq_batch_t *batch,
		     bwt_optarg_t *bwt_optarg, 
		     bwt_index_t *index, 
		     fastq_batch_t *unmapped_batch,
		     array_list_t *mapping_list) {
  
  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_batch(batch, bwt_optarg, index, 
			       unmapped_batch, mapping_list);
  }

  return bwt_map_inexact_batch(batch, bwt_optarg, index, 
			       unmapped_batch, mapping_list);  
}

//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------

size_t bwt_map_exact_seq(char *seq, 
			 bwt_optarg_t *bwt_optarg, 
			 bwt_index_t *index, 
			 array_list_t *mapping_list) {
  //printf("\tIn function...\n");
  unsigned int len = strlen(seq);
  int start = 0;
  int end = len - 1;
  char *code_seq = (char *) malloc(len * sizeof(char));
  result result;
  size_t l_aux, k_aux;
  unsigned int num_mappings = 0;
  char plusminus[2] = "-+";
  size_t idx, key, error, pos;
  alignment_t *alignment;
  char *cigar_p, *seq_dup;
  unsigned int start_mapping;
  //char *chromosome = (char *)malloc(sizeof(char)*10);

  replaceBases(seq, code_seq, len);
  
  //  printf("---> EXACT Search Read (%d): %s\n", len, seq);
  
  for (short int type = 1; type >= 0; type--) {
    result.k = 0;
    result.l = index->h_O.siz - 2;
    result.start = start;
    result.end = end;
    if (type == 1) {
      result.pos = end;
      BWExactSearchBackward(code_seq, &index->h_C, &index->h_C1, &index->h_O, &result);
    } else {
      result.pos = start;
      BWExactSearchForward(code_seq, &index->h_rC, &index->h_rC1, &index->h_rO, &result);
    }
      
    k_aux = result.k;
    l_aux = result.l;
    //printf("\tk=%d - l=%d\n", k_aux, l_aux);      
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
	  start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	  /*printf("\tStrand:%c\tchromosome:%d\tStart:%u\n",plusminus[type],
	    idx,
	    start_mapping);*/
	  
	  cigar_p = (char *)malloc(sizeof(char)*len);
	  sprintf(cigar_p, "%d=\0", len);
	  
	  seq_dup = (char *)malloc(sizeof(char)*(len + 1));
	  memcpy(seq_dup, seq, len + 1);
	  // save all into one alignment structure and insert to the list
	  alignment = alignment_new();
	  alignment_init_single_end(NULL, seq_dup, NULL, !type, 
				    idx - 1,				      
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
  
  free(code_seq);
  //printf("\tOut function\n");
  return num_mappings;
}

//-----------------------------------------------------------------------------

size_t bwt_map_exact_read(fastq_read_t *read, 
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index, 
			  array_list_t *mapping_list) {
  
  size_t padding = array_list_size(mapping_list);
  
  size_t num_mappings = bwt_map_exact_seq(read->sequence, 
					  bwt_optarg, index, 
					  mapping_list);
  alignment_t *mapping;
  
  unsigned int header_len, quality_len;
  
  for (int i = 0; i < num_mappings ; i++) {
    //printf("access to pos:%d-padding:%d\n", i + padding, padding);
    mapping = (alignment_t *) array_list_get(i + padding, mapping_list);
    if (mapping != NULL) {
      header_len = strlen(read->id) + 1;
      mapping->query_name = (char *)malloc(sizeof(char)*header_len);
      get_to_first_blank(read->id, header_len, mapping->query_name);
      //printf("--->header id %s to %s\n",read->id, header_id);
      //free(read->id);
      
      quality_len = strlen(read->quality) + 1;
      mapping->quality = (char *)malloc(sizeof(char)*quality_len);
      memcpy(mapping->quality, read->quality, quality_len);
      //mapping->quality = read->quality;
    } else {
      printf("Error to extract item\n");
    }
    //printf("Store Ok!\n"); 
    //printf("\tSEQ:%s-ID:%s\n", mapping->sequence, mapping->query_name);
  }
  
  return num_mappings;
}

//-----------------------------------------------------------------------------

size_t bwt_map_exact_seqs(char **seqs, 
			  size_t num_reads,
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index, 
			  char *out_status,
			  array_list_t *mapping_list) {
  
  unsigned int num_threads = bwt_optarg->num_threads;
  array_list_t **individual_mapping_list_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  const size_t chunk = MAX(1, num_reads/( num_threads*10));
  size_t num_mappings = 0; 

  short int th_id;
  alignment_t *mapping;
  
  for(int i = 0; i < num_threads; i++){
   individual_mapping_list_p[i]  = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }
  
  omp_set_num_threads(num_threads);
  
  //printf("Initialization ok!Process %d reads x %d threads\n", num_reads, num_threads);
  
  #pragma omp parallel for private(th_id, num_mappings) schedule(dynamic, chunk)
  for(int i = 0; i < num_reads; i++){
    th_id = omp_get_thread_num();
    //printf("I'm thread %d and process read %d\n", th_id,  i);
    num_mappings = bwt_map_exact_seq(seqs[i], bwt_optarg, index, individual_mapping_list_p[th_id]);
    if(num_mappings){
      out_status[i] = 1;
    }else{
      out_status[i] = 0;
    }
    num_mappings = 0;
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

  return array_list_size(mapping_list);
}


//-----------------------------------------------------------------------------

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
}

//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------

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
  char *code_seq = &seq[seq_start];/*(char *)malloc(sizeof(char)*(seq_end - seq_start + 1)); //= &seq[seq_start];
  memcpy(code_seq, &seq[seq_start], seq_end - seq_start);*/
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
	    /*printf("Strand:%c\tchromosome:%s\tStart:%u\tend:%u\n",plusminus[type],
		    index->karyotype.chromosome + (idx-1) * IDMAX,
		    start_mapping, start_mapping + len);
	    */
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
	    /*printf("\t%s\t%c\t%s %u %s error: %s, pos: %i, base: %i ",
	      "nothing", plusminus[type],
	      index->karyotype.chromosome + (idx-1) * IDMAX,
	      index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]),
	      seq, bwt_error_type(r->err_kind[0]), r->position[0], r->base[0]);
	    */	  
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

//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seq(char *seq, 
			   bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, 
			   array_list_t *mapping_list) {
  

  //printf("\tIn function...\n");
  /*unsigned int len = strlen(seq);
  unsigned int start = 0;
  unsigned int end = len - 1;*/
  char *seq_dup;
  /*char *codeSeq = (char *) calloc(len, sizeof(char));
  replaceBases(seq, codeSeq, len);*/

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
	    /*printf("Value idx = %d\n", idx);
	      printf("\t%s\t%c\t%s %u %s error: %s, pos: %i, base: %i cigar: %s\n",
	      "nothing", plusminus[type],
	      index->karyotype.chromosome + (idx-1) * IDMAX,
	      index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]),
	      seq, bwt_error_type(r->err_kind[0]), r->position[0], r->base[0], cigar);
	    */  
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
	    }/*else{
	      //printf("ERROR: Cigar bad generated\n");
	      error_debug = 0;
	    }*/
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

size_t bwt_map_inexact_read(fastq_read_t *read, 
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    array_list_t *mapping_list) {
  
  size_t padding = array_list_size(mapping_list);
  
  size_t num_mappings = bwt_map_inexact_seq(read->sequence, 
					    bwt_optarg, index, 
					    mapping_list);
  //printf("padding:%d - num_mappings:%d\n", padding, num_mappings);
  alignment_t *mapping;
  unsigned int header_len, quality_len;
  
  for (int i = 0; i < num_mappings; i++) {
    mapping = (alignment_t *) array_list_get(i + padding, mapping_list);
    
    if(mapping != NULL){
      header_len = strlen(read->id) + 1;
      mapping->query_name = (char *)malloc(sizeof(char)*header_len);
      get_to_first_blank(read->id, header_len, mapping->query_name);
      //printf("--->header id %s to %s\n",read->id, header_id);
      //free(read->id);
      
      quality_len = strlen(read->quality) + 1;
      mapping->quality = (char *)malloc(sizeof(char)*quality_len);
      memcpy(mapping->quality, read->quality, quality_len);
      //mapping->quality = read->quality;
    }else{
      printf("Error to extract item\n");
    }
  }
  
  return num_mappings;
}

//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seqs(char **seqs, 
			    size_t num_reads,
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    char *out_status,
			    array_list_t *mapping_list) {
  
  unsigned int num_threads = bwt_optarg->num_threads;
  array_list_t **individual_mapping_list_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  const unsigned int chunk = MAX(1, num_reads/( num_threads*10));
  size_t num_mappings = 0; 
  size_t num_mappings_tot = 0;
  short int th_id;
  alignment_t *mapping;
  
  
  for (int i = 0; i < num_threads; i++) {
    individual_mapping_list_p[i]  = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }
  
  omp_set_num_threads(num_threads);
  
  //printf("Initialization ok!Process %d reads x %d threads\n", num_reads, num_threads);
  
  #pragma omp parallel for private(th_id, num_mappings) schedule(dynamic, chunk)
  for (int i = 0; i < num_reads; i++) {
    th_id = omp_get_thread_num();
    //printf("I'm thread %d and process read %d\n", th_id,  i);
    num_mappings = bwt_map_inexact_seq(seqs[i], bwt_optarg, index, individual_mapping_list_p[th_id]);
    if(num_mappings == 0){
      out_status[i] = 0;
    }else{
      out_status[i] = 1;
    }
  }

  
  //printf("Process reads ok! Total mappings %d\n\n", num_mappings_tot);
  
  //Join results
  for (int i = 0; i < num_threads; i++) {
    num_mappings = array_list_size(individual_mapping_list_p[i]);
   
    for (int j = 0; j < num_mappings; j++) {
      mapping = (alignment_t *)array_list_get(j, individual_mapping_list_p[i]);
      array_list_insert( mapping, mapping_list);
    }
  }
  
  free(individual_mapping_list_p);
  
  return array_list_size(mapping_list);
}

//-----------------------------------------------------------------------------

size_t bwt_map_inexact_batch(fastq_batch_t *batch,
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
 
    bwt_map_inexact_seq(&(batch->seq[batch->data_indices[i]]), 
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

}

//-----------------------------------------------------------------------------

void bwt_map_inexact_batch_by_filter(fastq_batch_t *batch,
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index,
				     size_t num_selected, size_t *selected,
				     size_t num_lists, array_list_t **lists,
				     size_t *num_mapped, size_t *mapped,
				     size_t *num_unmapped, size_t *unmapped) {
  
  size_t num_mappings = 0, total = 0, header_len;
  size_t num_reads = batch->num_reads;

  alignment_t *alignment;
  
  *num_mapped = 0;
  *num_unmapped = 0;

  if (selected == NULL || num_selected == num_reads) {

    for (size_t i = 0; i < num_reads; i++) {
      num_mappings = bwt_map_inexact_seq(&(batch->seq[batch->data_indices[i]]), 
					 bwt_optarg, index, 
					 lists[i]);

      if (num_mappings > 0) {
	mapped[(*num_mapped)++] = i;
	array_list_set_flag(1, lists[i]);

	for (size_t j = 0; j < num_mappings; j++) {
	  alignment = (alignment_t *) array_list_get(j, lists[i]);

	  header_len = batch->header_indices[i+1] - batch->header_indices[i] + 1;
	  alignment->query_name = (char *) malloc(sizeof(char) * header_len);
	  get_to_first_blank(&(batch->header[batch->header_indices[i]]), header_len, alignment->query_name);

	  //	  alignment->query_name = strdup(&(batch->header[batch->header_indices[i]]));
	  alignment->quality = strdup(&(batch->quality[batch->data_indices[i]]));
	}
      } else {
	unmapped[(*num_unmapped)++] = i;
	array_list_set_flag(0, lists[i]);
      }
    }
  } else {
    printf("not yet implemented, true filter\n");
    exit(-1);
  }
}

//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seeds_seq(char *seq, size_t seed_size, size_t min_seed_size,
				 bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				 array_list_t *mapping_list) {
  
  size_t len = strlen(seq);
  size_t offset, num_seeds = len / seed_size;

  char *code_seq = (char *) calloc(len, sizeof(char));

  replaceBases(seq, code_seq, len);

  // first 'pasada'
  offset = 0;
  for (size_t i = 0; i < num_seeds; i++) {
    //printf("1, seed %d: start = %d, end = %d\n", i, offset, offset + seed_size - 1);
    bwt_map_inexact_seed(code_seq, offset, offset + seed_size - 1,
			 bwt_optarg, index, mapping_list);
    offset += seed_size;
  }

  // special processing for the last seed !!
  if (len % seed_size >= min_seed_size) {
    //printf("1', : start = %d, end = %d\n", offset, len - 1);
    bwt_map_inexact_seed(code_seq, offset, len - 1,
			 bwt_optarg, index, mapping_list);
  }

  // second 'pasada', shifting seeds by (seed_size / 2)
  offset = seed_size / 2;
  num_seeds = (len - seed_size / 2) / seed_size;
  for (size_t i = 0; i < num_seeds; i++) {
    //printf("2, seed %d: start = %d, end = %d\n", i, offset, offset + seed_size - 1);
    bwt_map_inexact_seed(code_seq, offset, offset + seed_size - 1,
			 bwt_optarg, index, mapping_list);
    offset += seed_size;
  }

  // again, special processing for the last seed !!
  if ((len - seed_size / 2) % seed_size >= min_seed_size) {
    //printf("2',: start = %d, end = %d\n", offset, len - 1);
    bwt_map_inexact_seed(code_seq, offset, len - 1,
			 bwt_optarg, index, mapping_list);
  }

  //printf("mapping list size = %d\n", array_list_size(mapping_list));
  free(code_seq);
  //  exit(-1);
  return array_list_size(mapping_list);
  
}

size_t bwt_map_exact_seeds_seq(char *seq, size_t seed_size, size_t min_seed_size,
				 bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				 array_list_t *mapping_list) {
  
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
  return array_list_size(mapping_list);
  
}
//-----------------------------------------------------------------------------
// CAL functions
//-----------------------------------------------------------------------------

int print_item(void *item, void *dummy){
        short_cal_t *coordenate_p = (short_cal_t *)item;
	printf("[%d-%d]-> ",  coordenate_p->start, coordenate_p->end);

        return 0;
}

int cal_location_compare(short_cal_t *a, short_cal_t *b) {
	int res = 1;

	if(a->start == b->start){
		if(a->end == b->end){
			res = 0;
		}else if(a->end < b->end){
			res = -1;
		}else{
			res = 1;
		}
	}else if(a->start < b->start){
		res = -1;
	}

	return res;

}

short_cal_t* cal_location_dup(short_cal_t *a) {
	short_cal_t *exon_location_p = (short_cal_t *)malloc(sizeof(short_cal_t));
	exon_location_p->end = a->end;
	exon_location_p->start = a->start;

	return exon_location_p;
}

void my_cp_list_append(cp_list* list_p, size_t start, size_t end, size_t max_cal_distance){
  unsigned char actualization = 0;
  short_cal_t *item, *item_aux, *new_item_p;
  cp_list_iterator itr;
  
  cp_list_iterator_init(&itr, list_p, COLLECTION_LOCK_NONE);
	
  if (cp_list_item_count(list_p) <= 0) {
    new_item_p = short_cal_new(start, end);
    cp_list_insert(list_p, new_item_p);
  } else {
    item = cp_list_iterator_curr(&itr);
    while ((item != NULL )) {
      if (start < item->start) {
	if (end + max_cal_distance < item->start) {
	  /*********************************************
	   *    Case 1: New item insert before item.   *
           *                                           *
           *        new item     item                  *
           *       |-------| |--------|                *
           ********************************************/
	  new_item_p = short_cal_new(start, end);
	  cp_list_iterator_insert(&itr, new_item_p);
	  cp_list_iterator_prev(&itr);
	} else {
	  /********************************************
           *  Case 2: Actualization item start        *
           *           new item                       *
           *          |-------|   item                *
           *                   |--------|             *                            
           ********************************************/
	  item->start = start;
	  if (end > item->end) {
	    /**************************************************
             *  Case 3: Actualization item start and item end *
             *          new item                              *
             *         |------------|                         *
             *              item                              *    
             *           |--------|                           *                                    
             **************************************************/
	    item->end = end;
	    actualization = 1;
	  }
	}
	break;
      }else {
	if(end <= item->end){
	  /**************************************************                                       
           *  Case 4: The new item don't insert in the list *                             
           *              item                              * 
           *         |-------------|                        * 
           *             new item                           * 
           *            |--------|                          * 
           **************************************************/
	  break;
	}else if (item->end + max_cal_distance >= start){
	  /********************************************                                              
           *  Case 5: Actualization item end          *
           *            item                          *                                              
           *          |-------| new item              *                                            
           *                 |--------|               *                                              
           ********************************************/
	  item->end = end;
	  actualization = 1;
	  break;
	}
      }//End else

      //continue loop...
      cp_list_iterator_next(&itr);
      item = cp_list_iterator_curr(&itr);

    }//End while

    if( item == NULL ){
      /******************************************************* 
       * Case 5: Insert new item at the end of the list      * 
       *                 item    new item                    * 
       *              |-------| |--------|                   *    
       *******************************************************/
      new_item_p = short_cal_new(start, end);
      cp_list_append(list_p, new_item_p);
    }
    
    if(actualization == 1){
      //printf("\tActualization RIGHT items (Next). Current item [%d-%d]\n", item->start, item->end);
      cp_list_iterator_next(&itr);
      item_aux = cp_list_iterator_curr(&itr);
      while(item_aux != NULL){
	//printf("\t\tFusion right items. item->end=%d < item_aux->start=%d?\n", item->end, item_aux->start);
	if (item->end + max_cal_distance< item_aux->start){
	  //printf("\t\tSTOP Actualization\n");
	  break;
	}else{
	  //printf("\t\tCONTINUE Actualization. item->end=%d < item_aux->end=%d?\n", item->end, item_aux->end);
	  if( item->end < item_aux->end ){
	    //printf("\t\tActualization end value %d\n", item_aux->end);
	    item->end = item_aux->end;
	  }
          //printf("\t\tDelete item %d-%d\n", item_aux->start, item_aux->end);
	  cp_list_iterator_remove(&itr);
	  //printf("\t\tDelete OK!\n");
	  item_aux = cp_list_iterator_curr(&itr);
	}                                                                                             
      }
    }
  }//end first else
}

size_t bwt_generate_cal_list_linkedlist(array_list_t *mapping_list,
					cal_optarg_t *cal_optarg,
					array_list_t *cal_list){

  short_cal_t *short_cal_p;
  region_t *region;
  size_t min_cal_size = cal_optarg->min_cal_size;
  size_t max_cal_distance  = cal_optarg->max_cal_distance;
  size_t num_mappings = array_list_size(mapping_list);
  size_t chromosome_id;
  short int strand;
  size_t start, end;  
  cp_list_iterator itr;

  const unsigned char nstrands = 2;
  const unsigned char nchromosomes = 30;
  cp_list ***cals_list = (cp_list ***)malloc(sizeof(cp_list **)*nstrands);

  for (unsigned int i = 0; i < nstrands; i++) {
    cals_list[i] = (cp_list **)malloc(sizeof(cp_list *)*nchromosomes);
    for (unsigned int j = 0; j < nchromosomes; j++) {
      cals_list[i][j] = cp_list_create_list(COLLECTION_MODE_NOSYNC |
					    /*COLLECTION_MODE_COPY |*/
					    COLLECTION_MODE_DEEP |
					    COLLECTION_MODE_MULTIPLE_VALUES,
					    (cp_compare_fn) cal_location_compare,
					    NULL/*(cp_copy_fn) cal_location_dup*/,
					    (cp_destructor_fn) short_cal_free);
    }
  }
    
  // from the results mappings, generates CALs
  for (unsigned int m = 0; m < num_mappings; m++) {
    region = array_list_get(m, mapping_list);
    chromosome_id = region->chromosome_id;
    strand = region->strand;
    start = region->start;
    end = region->end;

    /*
    if (chromosome_id == 3) {
      printf("num_mapping %d\tchromosome = %i\tstrand = %d, start = %d, end = %d\n", m, chromosome_id, strand, start, end);
    }
    */

    my_cp_list_append(cals_list[strand][chromosome_id], start, end, max_cal_distance);
  }

  //Store CALs in Array List for return results
  for (unsigned int j = 0; j < nchromosomes; j++) {
    for (unsigned int i = 0; i < nstrands; i++) {
      cp_list_iterator_init(&itr, cals_list[i][j], COLLECTION_LOCK_NONE);
      short_cal_p = cp_list_iterator_curr(&itr);
      
      while ((short_cal_p != NULL )) {
	if (short_cal_p->end - short_cal_p->start + 1 >= min_cal_size) {
	  array_list_insert(cal_new(j, i, short_cal_p->start, short_cal_p->end), cal_list);
	}
	//short_cal_free(short_cal_p);
	cp_list_iterator_next(&itr);
	short_cal_p = cp_list_iterator_curr(&itr);
      }
      cp_list_destroy(cals_list[i][j]);
    }
  }
  
  for (unsigned int i = 0; i < nstrands; i++) {
    free(cals_list[i]);
  }

  free(cals_list);

  return array_list_size(cal_list);  
}


//-----------------------------------------------------------------------------

size_t bwt_generate_cal_list(array_list_t *mapping_list,
			     cal_optarg_t *cal_optarg,
			     array_list_t *cal_list) {
  cal_t *cal, *cal1;
  region_t *region;
  short int extend, added;
  size_t cal_idx;
  size_t min_cal_size = cal_optarg->min_cal_size;
  size_t max_cal_distance  = cal_optarg->max_cal_distance;
  size_t num_mappings = array_list_size(mapping_list);
  size_t num_cals;
  size_t chromosome_id;
  short int strand;
  size_t start, end;
  
  int print = 0;
  
  //printf("Start generate CALs: %d mappings\n", num_mappings);
  // from the results mappings, generates CALs
  for (unsigned int m = 0; m < num_mappings; m++) {
    region = array_list_get(m, mapping_list);

    chromosome_id = region->chromosome_id;
    strand = region->strand;
    start = region->start;
    end = region->end;

    /*if ( (strand == 1 && chromosome_id == 1 && (start >= 274839 && start <= 274859)) ||
	   (strand == 1 && chromosome_id == 1 && (end >= 274839 && end <= 274859)) ) {
	printf("----> mapping, start = %d, end = %d\n", start, end);
      }*/
    //printf("----> mapping, start = %d, end = %d\n", start, end);
    
    extend = 0;
    added = 0;
    num_cals = array_list_size(cal_list);
    //printf("Process region strand %d chromosome %d [%d-%d]\n", strand, chromosome_id, start, end);
    for (unsigned int c = 0; c < num_cals; c++) {
      cal = array_list_get(c, cal_list);

      /*if ( (strand == 1 && cal->chromosome_id == 1 && (cal->start >= 274839 && cal->start <= 274859)) ||
	   (strand == 1 && cal->chromosome_id == 1 && (cal->end >= 274839 && cal->end <= 274859)) ) {
	//printf("Cal %d - %d\n", cal->start);
	print = 1;
      }*/

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
	      cal_t *aux = array_list_remove_at(c, cal_list);
	      if (aux != NULL) { 
		cal_free(aux); 
	      }
	      break;
	    }
	  }
	}
      }
    } else {
      array_list_insert(cal_new(chromosome_id, strand, start, end),
			cal_list);
    }
  } // end for mappings
  //printf("End\n");
  
  num_cals = array_list_size(cal_list);
  // removing 'small' CALs from list
  for (int i = num_cals - 1; i >= 0 ; i--) {
    cal = array_list_get(i, cal_list);
    //printf("(%d - %d + 1= %d) < %d \n", cal->end, cal->start, cal->end - cal->start + 1, min_cal_size);
    if (cal->end - cal->start + 1 < min_cal_size) {
      cal = array_list_remove_at(i, cal_list);
      cal_free(cal);
    }/*else{
      printf("Not Remove strand %d chromosome %d CAL %d-%d\n", cal->strand, cal->chromosome_id, cal->start, cal->end);
      }*/
  }
  
  return array_list_size(cal_list);
}

//

//-----------------------------------------------------------------------------

size_t bwt_find_cals_from_seq(char *seq, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      cal_optarg_t *cal_optarg, 
			      array_list_t *cal_list) {
  // create list where to store the seed mappings
  array_list_t *mapping_list = array_list_new(100000, 
					      1.25f, 
					      COLLECTION_MODE_SYNCHRONIZED); 
  //printf("Array list complete\n");

  // map seeds
  unsigned int num_mappings = bwt_map_exact_seeds_seq(seq, cal_optarg->seed_size,
				         		cal_optarg->min_seed_size,
							bwt_optarg, index, mapping_list);
  

  //printf("num_mappings = %d, mapping_list_size = %d\n", num_mappings, array_list_size(mapping_list));

  // generate cals from mapped seeds
  unsigned int num_cals = bwt_generate_cal_list(mapping_list,
						cal_optarg, cal_list);
  // free memory for mapping list
  array_list_free(mapping_list, NULL);

  return num_cals;
}

//-----------------------------------------------------------------------------

size_t bwt_find_cals_from_read(fastq_read_t *read, 
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       cal_optarg_t *cal_optarg, 
			       array_list_t *cal_list) {

  return bwt_find_cals_from_seq(read->sequence, 
				    bwt_optarg, index,
				    cal_optarg, cal_list);
}

//-----------------------------------------------------------------------------

size_t bwt_find_cals_from_seqs(char **seqs, 
			       size_t num_reads,
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       cal_optarg_t *cal_optarg, 
			       char *out_status,
			       array_list_t *cal_list) {
  
    unsigned int num_seeds = 0;
    unsigned int num_threads = bwt_optarg->num_threads;
    const unsigned int chunk = MAX(1, num_reads/( num_threads*10));
    //    seed_t **seeds_p = (seed_t **) calloc(num_reads, sizeof(seed_t *));
    size_t num_mappings;
    fastq_read_t fq_read;
    //    seed_t *seeds;
    unsigned int th_id;
    size_t num_cals, total_cals = 0;

    array_list_t **mapping_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
    array_list_t **cals_reads = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);

    for(int i = 0; i < num_reads; i++){
        mapping_list[i] = array_list_new(1000,
                                         1.25f,
                                         COLLECTION_MODE_SYNCHRONIZED);

        cals_reads[i] = array_list_new(100,
                                       1.25f,
                                       COLLECTION_MODE_SYNCHRONIZED);
    }

    omp_set_num_threads(num_threads);

    size_t (*map_seeds)(char *, size_t, size_t, bwt_optarg_t *, bwt_index_t *, array_list_t *);

    if(cal_optarg->num_errors == 0){
      map_seeds = &bwt_map_exact_seeds_seq;
    }else{
      map_seeds = &bwt_map_inexact_seeds_seq;
    }
    
    #pragma omp parallel for private(num_cals) reduction(+:total_cals)
    for(int i = 0; i < num_reads; i++){
        //th_id = omp_get_thread_num();
        // map seeds
        mapping_list[i] = array_list_new(1000,
                                         1.25f,
                                         COLLECTION_MODE_SYNCHRONIZED);

        (*map_seeds)(seqs[i], cal_optarg->seed_size,
                     cal_optarg->min_seed_size,
                     bwt_optarg, index, mapping_list[i]);

        //printf("num_mappings = %d, mapping_list_size = %d\n", num_mappings, array_list_size(mapping_list[th_id]));

        // generate cals from mapped seeds
        num_cals = bwt_generate_cal_list(mapping_list[i],
                                         cal_optarg, cals_reads[i]);


        if(!num_cals){
          out_status[i] = 0;
        }else{
          out_status[i] = 1;
          total_cals += num_cals;
        }
    }
    
    for(int i = 0; i < num_reads; i++){
        array_list_insert((void *)cals_reads[i], cal_list);
        // free memory for mapping list
        array_list_free(mapping_list[i], NULL);
    }

    free(mapping_list);

    return total_cals;

}

//-----------------------------------------------------------------------------

size_t bwt_find_cals_from_batch(fastq_batch_t *batch,
				bwt_optarg_t *bwt_optarg, 
				bwt_index_t *index, 
				fastq_batch_t *unmapped_batch,
				cal_optarg_t *cal_optarg, 
				array_list_t *cal_list) {
  unsigned int num_mappings;
  size_t num_reads = batch->num_reads;
  size_t *unmapped_index = (size_t *)malloc(sizeof(size_t)*num_reads);
  unsigned int num_unmapping = 0;
  unsigned int num_cals = 0;
  size_t read_pos;
  unsigned int length_header;
  unsigned int length_seq;
  size_t id_read;
  size_t total_cals = 0;

  array_list_t **mapping_list = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t **cals_reads = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);

  for(int i = 0; i < num_reads; i++){
      mapping_list[i] = array_list_new(1000,
                                        1.25f,
                                        COLLECTION_MODE_SYNCHRONIZED);

      cals_reads[i] = array_list_new(100,
                                      1.25f,
                                      COLLECTION_MODE_SYNCHRONIZED);
  }

  //  printf("Array list complete\n");
  
  size_t (*map_seeds)(char *, size_t, size_t, bwt_optarg_t *, bwt_index_t *, array_list_t *);

  if(cal_optarg->num_errors == 0){
    map_seeds = &bwt_map_exact_seeds_seq;
  }else{
    map_seeds = &bwt_map_inexact_seeds_seq;
  }

  #pragma omp parallel for private(num_mappings, num_cals) reduction(+:total_cals)
  for (size_t i = 0; i < num_reads; i++) {
    //create list where to store the seed mappings
    mapping_list[i] = array_list_new(1000,
                                     1.25f,
                                     COLLECTION_MODE_SYNCHRONIZED);

    num_mappings = 0;
    num_mappings = map_seeds(&(batch->seq[batch->data_indices[i]]), cal_optarg->seed_size,
                                                cal_optarg->min_seed_size,
                                                bwt_optarg, index, mapping_list[i]);
    if(!num_mappings){
      unmapped_index[num_unmapping] = i;
      num_unmapping++;
    }else{
      // generate cals from mapped seeds
      num_cals = bwt_generate_cal_list(mapping_list[i],
                                       cal_optarg, cals_reads[i]);
      if(!num_cals){
        unmapped_index[num_unmapping] = i;
        num_unmapping++;
      }else{
        total_cals += num_cals;
      }
    }
  }
  
  // map seeds
  read_pos = 0;
  unmapped_batch->header_indices[read_pos] = 0;
  unmapped_batch->data_indices[read_pos] = 0;
  for (int i = 0; i < num_unmapping; i++) {
    id_read = unmapped_index[i];
    length_header = batch->header_indices[id_read] - batch->header_indices[id_read - 1];
    length_seq = batch->data_indices[id_read] - batch->data_indices[id_read - 1];
    read_pos++;
    memcpy(&(unmapped_batch->header[unmapped_batch->header_indices[read_pos - 1]]), &(batch->header[batch->header_indices[id_read - 1]]), length_header);
    memcpy(&(unmapped_batch->seq[unmapped_batch->data_indices[read_pos - 1]]),      &(batch->seq[batch->data_indices[id_read - 1]]),      length_seq);
    memcpy(&(unmapped_batch->quality[unmapped_batch->data_indices[read_pos - 1]]),  &(batch->quality[batch->data_indices[id_read - 1]]),  length_seq);

    unmapped_batch->data_indices[read_pos]   = unmapped_batch->data_indices[read_pos - 1]   + length_seq;
    unmapped_batch->header_indices[read_pos] = unmapped_batch->header_indices[read_pos - 1] + length_header;
  }

  unmapped_batch->num_reads = read_pos;


  // free memory for mapping list
  for(int i = 0; i < num_reads; i++){
    array_list_insert((void *)cals_reads[i], cal_list);
    // free memory for mapping list
    array_list_free(mapping_list[i], NULL);
  }

   free(mapping_list);

  return total_cals;
  
}

//-----------------------------------------------------------------------------
