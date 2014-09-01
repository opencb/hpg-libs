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

//------------------------------------------------------------------------------

size_t bwt_map_exact_seed(uint8_t *seq, size_t seq_len, 
			  size_t start, size_t end,
			  bwt_optarg_t *bwt_optarg,
			  bwt_index_t *index,
			  array_list_t *mapping_list,
			  int id);

//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seq(char *seq, 
			   bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, 
			   array_list_t *mapping_list);

/*size_t bwt_map_inexact_read(fastq_read_t *read, 
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    array_list_t *mapping_list);
*/
size_t bwt_map_inexact_seqs(char **seqs, 
			    size_t num_reads,
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    char *out_status,
			    array_list_t *mapping_list);

//------------------------------------------------------------------------------

size_t bwt_map_inexact_seed(char *seq, size_t seq_len, 
			    size_t start, size_t end,
			    bwt_optarg_t *bwt_optarg,
			    bwt_index_t *index,
			    array_list_t *mapping_list,
			    int seed_id);

//------------------------------------------------------------------------------

char *bwt_error_type(char error_kind);

//------------------------------------------------------------------------------

size_t __bwt_map_inexact_read(fastq_read_t *read, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list);

//------------------------------------------------------------------------------

//alignment_t* add_optional_fields(alignment_t *alignment, size_t n_mappings, int read_length);

//------------------------------------------------------------------------------

void *__bwt_generate_anchor_list(size_t k_start, size_t l_start, int len_calc, bwt_optarg_t *bwt_optarg, 
				 bwt_index_t *index, int type, array_list_t *anchor_list, int type_anchor,
				 int dsp);

//------------------------------------------------------------------------------

void append_seed_region_linked_list(linked_list_t* sr_list,
				    size_t read_start, size_t read_end, 
				    size_t genome_start, size_t genome_end, 
				    int seed_id, int pos_err, int type_err, int type_seeds);

//------------------------------------------------------------------------------

void seed_region_select_linked_list(linked_list_t* sr_list, linked_list_t* sr_duplicate_list, 
				    size_t read_start, size_t read_end,
				    size_t genome_start, size_t genome_end,
				    int seed_id, int pos_err, int type_err, 
				    unsigned char *seeds_ids_array);

//------------------------------------------------------------------------------

void my_cp_list_append_linked_list(linked_list_t* list_p, region_t *region, size_t max_cal_distance, int max_seeds);

//------------------------------------------------------------------------------
// Paratemers for the candidate alignment localizations (CALs)
//------------------------------------------------------------------------------

bwt_err_t *bwt_err_new(int pos, char name) {
  bwt_err_t *bwt_err = (bwt_err_t *)malloc(sizeof(bwt_err_t));
  bwt_err->pos = pos;
  bwt_err->name = name;

  return bwt_err;
  
}


void bwt_err_free(bwt_err_t *p) {
  if (p) { free(p); }
}


cal_optarg_t *cal_optarg_new(const size_t min_cal_size, 
			     const size_t max_cal_distance, 
			     const size_t num_seeds,
			     const size_t min_num_seeds_in_cal,
			     const size_t seed_size,
			     const size_t min_seed_size,
			     const size_t num_errors,
			     const size_t max_intron_size,
			     const size_t min_intron_size) {
			       
  cal_optarg_t *cal_optarg_p = (cal_optarg_t *)malloc(sizeof(cal_optarg_t));

  cal_optarg_p->min_cal_size = min_cal_size;
  cal_optarg_p->max_cal_distance = max_cal_distance;
  cal_optarg_p->num_seeds = num_seeds;
  cal_optarg_p->min_num_seeds_in_cal = min_num_seeds_in_cal;
  cal_optarg_p->min_seed_size = min_seed_size;
  cal_optarg_p->seed_size = seed_size;
  cal_optarg_p->num_errors = num_errors;
  
  cal_optarg_p->max_intron_size = max_intron_size;
  cal_optarg_p->min_intron_size = min_intron_size;

  return cal_optarg_p;
}
    
void cal_optarg_free(cal_optarg_t *optarg) {
  free(optarg);
}

//------------------------------------------------------------------------------

simple_seed_t *simple_seed_new(size_t read_start, size_t read_end) {
  simple_seed_t *ss = (simple_seed_t *)malloc(sizeof(simple_seed_t));
  ss->start = read_start;
  ss->end = read_end;

  return ss;
}

void simple_seed_free(simple_seed_t *simple_seed) {
  free(simple_seed);
}

//------------------------------------------------------------------------------

cal_t *cal_new(const size_t chromosome_id, 
	       const short int strand,
	       const size_t start, 
	       const size_t end,
	       const size_t num_seeds,
	       linked_list_t *sr_list,
	       linked_list_t *sr_duplicate_list) {
		 
  cal_t *cal = (cal_t *)malloc(sizeof(cal_t));  

  cal->chromosome_id = chromosome_id;
  cal->end = end;
  cal->start = start;
  cal->strand = strand;
  cal->num_seeds = num_seeds;
  cal->sr_list = sr_list;
  cal->sr_duplicate_list = sr_duplicate_list;
  cal->read_area = 0;
  cal->info = NULL;
  cal->fill_gaps = 0;
  //cal->num_targets;
  cal->type_seeds = 0;
  cal->candidates_seeds_start = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  cal->candidates_seeds_end = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  //TODO: Good CAL Selector??Mmm...
  //First, seed size sumation  

  //Second, read distance between the first element of the list and the last
  seed_region_t *s_first = linked_list_get_first(cal->sr_list);
  seed_region_t *s_last = linked_list_get_last(cal->sr_list);
  
  //printf("List Size %i:\n", cal->sr_list->size);
  if (s_first && s_last) {
    cal->read_area += (s_last->read_end - s_first->read_start);
    //printf("\t %i - %i = %i\n", s_last->read_end, s_first->read_start, cal->read_area);
  }

  return cal;
}

cal_t *cal_simple_new(const size_t chromosome_id, 
		      const short int strand,
		      const size_t start, 
		      const size_t end) {		 
  cal_t *cal = (cal_t *)malloc(sizeof(cal_t));  

  cal->chromosome_id = chromosome_id;
  cal->start = start;
  cal->end = end;
  cal->strand = strand;
  cal->sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);

  cal->info = NULL;

  return cal;

}

void cal_simple_free(cal_t *cal) {
  if (cal) {
    if (cal->sr_list) linked_list_free(cal->sr_list, (void *) seed_region_free);
    free(cal);
  }
}

void cal_free(cal_t *cal) {
  if (cal) {
    if (cal->sr_list) linked_list_free(cal->sr_list, (void *) seed_region_free);
    if (cal->sr_duplicate_list) linked_list_free(cal->sr_duplicate_list, (void *) seed_region_free);

    if (cal->candidates_seeds_start) array_list_free(cal->candidates_seeds_start, (void *)seed_region_free);
    if (cal->candidates_seeds_end) array_list_free(cal->candidates_seeds_end, (void *)seed_region_free);
    free(cal);
  }

}

void cal_print(cal_t *cal) {
  printf(" CAL (%c)[%lu:%lu-%lu]:\n", cal->strand==0?'+':'-', 
	 cal->chromosome_id, cal->start, cal->end);
  printf("\t SEEDS LIST: ");
  if (cal->sr_list == NULL || cal->sr_list->size == 0) {
    printf(" NULL\n");
  } else {
    for (linked_list_item_t *item = cal->sr_list->first; 
	 item != NULL; item = item->next) {
      seed_region_t *seed = item->item;
      printf(" (%i)[%lu|%i - %i|%lu] ", seed->id, seed->genome_start, seed->read_start, 
	     seed->read_end, seed->genome_end);
    }
    printf("\n");
  }
  /*
  if (cal->sr_duplicate_list != NULL && cal->sr_duplicate_list->size > 0) {
    printf("\t SEEDS DUPLICATE LIST: ");
    for (linked_list_item_t *item = cal->sr_duplicate_list->first; 
	 item != NULL; item = item->next) {
      seed_region_t *seed = item->item;
      printf(" [%lu|%lu - %lu|%lu] ", seed->genome_start, seed->read_start, 
	     seed->read_end, seed->genome_end);
    }
    printf("\n");
  }
  */
}

//------------------------------------------------------------------------------

seed_region_t *seed_region_new(int read_start, int read_end, size_t genome_start, 
			       size_t genome_end, int id, int pos_err, int type_err) {

  seed_region_t *seed_region = (seed_region_t *)malloc(sizeof(seed_region_t));

  seed_region->read_start = read_start;
  seed_region->read_end = read_end;
  seed_region->genome_start = genome_start;
  seed_region->genome_end = genome_end;
  seed_region->id = id;
  seed_region->info = NULL;
  seed_region->fusion_left = 0;
  seed_region->fusion_right = 0;

  seed_region->pos_err = pos_err;
  seed_region->type_err = type_err;
  seed_region->errors_list = NULL;

  return seed_region;
}

void seed_region_free(seed_region_t *seed_region) {
  if (seed_region->info) free(seed_region->info);
  if (seed_region->errors_list) {
    array_list_free(seed_region->errors_list, (void *)bwt_err_free);
    seed_region->errors_list = NULL;
  }

  free(seed_region);
}

void seed_region_simple_free(seed_region_t *seed_region) {
  free(seed_region);
}

//------------------------------------------------------------------------------

short_cal_t *short_cal_new(const size_t start, 
	                   const size_t end,
			   const size_t seq_start,
			   const size_t seq_end,
			   const size_t seq_len,
			   const int max_seeds,
			   const int id,
			   const int pos_err,
			   const int type_err) {

  short_cal_t *short_cal = (short_cal_t *)malloc(sizeof(short_cal_t));  
  
  short_cal->end = end;
  short_cal->start = start;
  short_cal->seq_len = seq_len;
  short_cal->num_seeds = 1;

  short_cal->seeds_ids_array = (unsigned char *)calloc(max_seeds, sizeof(unsigned char));

  if (short_cal->seeds_ids_array == NULL) LOG_FATAL("NO MORE MEMORY\n");
  if (id >= max_seeds) LOG_FATAL_F("STORAGE SEED ID OVERFLOW: %i > %i\n", id, max_seeds);

  short_cal->seeds_ids_array[id] = 1;

  short_cal->sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  seed_region_t *seed_region = seed_region_new(seq_start, seq_end, start, 
					       end, id, pos_err, type_err);

  //  printf("\tInsert New [Seed:=%lu-%lu]\n", seq_start, seq_end);
  linked_list_insert(seed_region, short_cal->sr_list);
  short_cal->sr_duplicate_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);

  return short_cal;
}


void short_cal_free(short_cal_t *short_cal){
  if (short_cal) {
    if (short_cal->seeds_ids_array) free(short_cal->seeds_ids_array);
    if (short_cal->sr_list) linked_list_free(short_cal->sr_list, (void *) seed_region_free);
    if (short_cal->sr_duplicate_list) linked_list_free(short_cal->sr_duplicate_list, (void *) seed_region_free);
    
    free(short_cal);
  }

}

//------------------------------------------------------------------------------

bwt_anchor_t *bwt_anchor_new(int strand, int chromosome, size_t start, size_t end, int type) {
  bwt_anchor_t *bwt_anchor = (bwt_anchor_t *)malloc(sizeof(bwt_anchor_t));

  bwt_anchor->strand = strand;
  bwt_anchor->chromosome = chromosome;
  bwt_anchor->start = start;
  bwt_anchor->end = end;
  bwt_anchor->type = type;

  return bwt_anchor;
}

void bwt_anchor_free(bwt_anchor_t *bwt_anchor) {
  free(bwt_anchor);
}
//------------------------------------------------------------------------------

region_t *region_bwt_new(const size_t chromosome_id, 
			 const short int strand,
			 const size_t start, 
			 const size_t end,
			 const size_t seq_start,
			 const size_t seq_end,
			 const size_t seq_len,
			 const int id) {
		 
  region_t *region = (region_t *) malloc(sizeof(region_t));  

  region->chromosome_id = chromosome_id;
  region->end = end;
  region->start = start;
  region->strand = strand;
  region->seq_start = seq_start;
  region->seq_end = seq_end;
  region->seq_len = seq_len;
  region->id = id;

  return region;
}


void region_bwt_free(region_t *region) {
  free(region);
}

//------------------------------------------------------------------------------

read_cals_t *read_cals_new(fastq_read_t *read) {

    read_cals_t *read_cals = (read_cals_t *) malloc(sizeof(read_cals_t));

    read_cals->cal_list = array_list_new(100000, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
    read_cals->read = read;

    return read_cals;
}

void read_cals_free(read_cals_t *read_cals){
  fastq_read_free(read_cals->read);
  array_list_free(read_cals->cal_list, (void *)cal_free);
  free(read_cals);
}

//------------------------------------------------------------------------------
//Options settings Burrows-Wheeler Transform 
//------------------------------------------------------------------------------

bwt_optarg_t *bwt_optarg_new(const size_t num_errors,
			     const size_t num_threads,
 			     const int filter_read_mappings, 
			     const int filter_seed_mappings) {
  bwt_optarg_t *bwt_optarg = (bwt_optarg_t *) calloc(1, sizeof(bwt_optarg_t));
  
  bwt_optarg->num_errors = num_errors;
  bwt_optarg->num_threads = num_threads;
  bwt_optarg->filter_read_mappings = filter_read_mappings;
  bwt_optarg->filter_seed_mappings = filter_seed_mappings;

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
char *strcpy_capitalize(char *dest, const char *src, size_t n) {
  size_t i;
  
  for (i = 0 ; i < n && src[i] != '\0' ; i++) {
    if (src[i] > 96) {
      dest[i] = src[i] - 32;
    } else {
      dest[i] = src[i];
    }
  }

  for ( ; i < n ; i++) {
    dest[i] = '\0';
  }

  return dest;

}

//-----------------------------------------------------------------------------

void bwt_cigar_cpy(alignment_t *mapping, char *quality) {
  
  unsigned int quality_type;
  size_t quality_len;

  quality_len = strlen(quality);
  quality_type = atoi(mapping->quality);
  //printf("Quality len from batch: %i\n", quality_len);
  free(mapping->quality);
  mapping->quality = (char *)malloc(sizeof(char)*(quality_len + 1));
  //printf("Read:%s\n", mapping->sequence);
/*
  if (quality_type == START_HARD_CLIPPING){
    if (mapping->seq_strand == 0) {
      memcpy(mapping->quality, 
	     quality + 1, 
	     quality_len - 1);
    } else {
      reverse_str(quality + 1,
		  mapping->quality, quality_len - 1);
    }
    mapping->quality[quality_len - 1] = '\0';	    
    //printf("HARD START : %s\n", mapping->quality);
  }else if(quality_type == END_HARD_CLIPPING){
    if (mapping->seq_strand == 0) {
      memcpy(mapping->quality, 
	     quality,
	     quality_len - 1);
    } else {
      reverse_str(quality,
		  mapping->quality, quality_len - 1);
    }
    mapping->quality[quality_len - 1] = '\0';
    //printf("HARD END : %s\n", mapping->quality);
  }else{*/
    //printf("ELSE....\n");
  if (mapping->seq_strand == 0) {
       memcpy(mapping->quality, quality, quality_len);
  } else {
       reverse_str(quality,
		   mapping->quality, quality_len);
  }
  //mapping->quality[quality_len] = '\0';
  //printf("(%i)NORMAL : %s\n", mapping->seq_strand, mapping->quality);
//}
  //array_list_insert( mapping, mapping_list);
}
//-----------------------------------------------------------------------------

void bwt_cigar_cpy_batch(alignment_t *mapping, size_t read_i, fastq_batch_t *batch) {

  unsigned int quality_type;
  size_t quality_len;
  quality_len = batch->data_indices[read_i + 1] - batch->data_indices[read_i] - 1;
  quality_type = atoi(mapping->quality);
  //printf("Quality len from batch: %i\n", quality_len);                                                                                                                                                                                                                                                                                    
  free(mapping->quality);
  mapping->quality = (char *)malloc(sizeof(char)*(quality_len + 1));
  //printf("Read:%s\n", mapping->sequence);                                                                                                                                                                                                                                                                                                 

  if (quality_type == START_HARD_CLIPPING){
    if (mapping->seq_strand == 0) {
      memcpy(mapping->quality ,
             &(batch->quality[batch->data_indices[read_i]]) + 1,
             quality_len - 1);
    } else {
      reverse_str(&(batch->quality[batch->data_indices[read_i]]) + 1,
                  mapping->quality, quality_len - 1);
    }
    mapping->quality[quality_len - 1] = '\0';
    //printf("HARD START : %s\n", mapping->quality);                                                                                                                                                                                                                                                                                        
  }else if(quality_type == END_HARD_CLIPPING){
    if (mapping->seq_strand == 0) {
      memcpy(mapping->quality,
             &(batch->quality[batch->data_indices[read_i]]),
             quality_len - 1);
    } else {
      reverse_str(&(batch->quality[batch->data_indices[read_i]]),
                  mapping->quality, quality_len - 1);
    }
    mapping->quality[quality_len - 1] = '\0';
    //printf("HARD END : %s\n", mapping->quality);                                                                                                                                                                                                                                                                                          
  }else{
    //printf("ELSE....\n");                                                                                                                                                                                                                                                                                                                 
    if (mapping->seq_strand == 0) {
      memcpy(mapping->quality, &(batch->quality[batch->data_indices[read_i]]), quality_len);
    } else {
      reverse_str(&(batch->quality[batch->data_indices[read_i]]),
                  mapping->quality, quality_len);
    }
    //mapping->quality[quality_len] = '\0';                                                                                                                                                                                                                                                                                                 
    //printf("(%i)NORMAL : %s\n", mapping->seq_strand, mapping->quality);                                                                                                                                                                                                                                                                   
  }
  //array_list_insert( mapping, mapping_list);                                                                                                                                                                                                                                                                                              
}

//-----------------------------------------------------------------------------

unsigned int alignmentcmp(alignment_t *alignment_1, alignment_t *alignment_2) {
  
  size_t cigar1_len = strlen(alignment_1->cigar);
  size_t cigar2_len = strlen(alignment_2->cigar);

  if (alignment_1->map_quality > alignment_2->map_quality) { return 1; }
  else if (alignment_1->map_quality < alignment_2->map_quality) { return 2; }
  else {
    if (alignment_1->num_cigar_operations == 1 ||
	alignment_2->num_cigar_operations == 1) {
      if (alignment_1->cigar[cigar1_len - 1] == '=' && 
	  alignment_2->cigar[cigar1_len - 1] == '=') { return 0; }
      
      if (alignment_1->cigar[cigar1_len - 1] == '=') { return 1; }
      else if (alignment_2->cigar[cigar2_len - 1] == '=') { return 2; }
      
      if (alignment_1->cigar[cigar1_len - 1] == 'M' && 
	  alignment_2->cigar[cigar1_len - 1] == 'M') { return 0; }
      
      if (alignment_1->cigar[cigar1_len - 1] == 'M') { return 1; }
      else if (alignment_2->cigar[cigar2_len - 1] == 'M') { return 2; }
      
    }else {
      if (alignment_1->num_cigar_operations < 
	  alignment_2->num_cigar_operations) {
	return 1;
      } else if (alignment_1->num_cigar_operations > 
		 alignment_2->num_cigar_operations) {
	return 2;
      } else { return 0; }
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------

char* reverse_str(char *src, char *dsp, size_t length) {
  size_t j = length - 1;
  size_t i = 0;

  memcpy(dsp, src, length);
  
  for (i = 0; i < length; i++) {
    dsp[i] = src[j--];
  }
  dsp[i] = '\0';
 
  return dsp;
}

//-----------------------------------------------------------------------------

void seq_reverse_complementary(char *seq, unsigned int len){
  unsigned int j = 0;
  
  char *seq_tmp = (char *)malloc(len*sizeof(char));
  memcpy(seq_tmp, seq, len);
  for (int i = len - 1; i >= 0; i--){
    if (seq_tmp[i] == 'A' || seq_tmp[i] == 'a') {
      seq[j] = 'T';
    }
    else if (seq_tmp[i] == 'C' || seq_tmp[i] == 'c'){
       seq[j] = 'G';
    }
    else if (seq_tmp[i] == 'G' || seq_tmp[i] == 'g'){
       seq[j] = 'C';
    }
    else if (seq_tmp[i] == 'T' || seq_tmp[i] == 't'){
       seq[j] = 'A';
    } else {
       seq[j] = 'N';
    }
    j++;  
  }

  free(seq_tmp);

}

//-----------------------------------------------------------------------------

bwt_index_t *bwt_index_new(const char *dirname, bool inverse_sa) {
  bwt_index_t *index  = (bwt_index_t*) malloc(sizeof(bwt_index_t));

  char *nucleotides = (char *) malloc(128 * sizeof(char));
  read_config(nucleotides, &(index->bwt_config.duplicate_strand), dirname);
  bwt_init_replace_table(&(index->bwt_config), nucleotides);

  index->karyotype    = exome_new();
  index->backward     = (bwt_index*) malloc(sizeof(bwt_index));
  index->forward      = (bwt_index*) malloc(sizeof(bwt_index));

  if (!index->bwt_config.duplicate_strand) {
    index->backward_rev = (bwt_index*) malloc(sizeof(bwt_index));
    index->forward_rev  = (bwt_index*) malloc(sizeof(bwt_index));
  }

  index->bwt_config.inverse_sa = inverse_sa;
  index->dirname = strdup(dirname);

  if (index->bwt_config.duplicate_strand) {
    load_bwt_index(NULL, index->backward, dirname, 1, inverse_sa, index->bwt_config);
    load_bwt_index(NULL, index->forward, dirname, 0, inverse_sa, index->bwt_config);
  } else {
    load_bwt_index(index->backward_rev, index->backward, dirname, 1, inverse_sa, index->bwt_config);
    load_bwt_index(index->forward_rev, index->forward, dirname, 0, inverse_sa, index->bwt_config);
  }

  load_exome_file(index->karyotype, dirname);

  return index;

}

//-----------------------------------------------------------------------------

void bwt_index_free(bwt_index_t *index) {

  if (index == NULL) return;

  free(index->dirname);

  if (index->bwt_config.duplicate_strand) {
    free_bwt_index(NULL, index->backward, index->bwt_config.inverse_sa);
    free_bwt_index(NULL, index->forward, index->bwt_config.inverse_sa);
  } else {
    free_bwt_index(index->backward_rev, index->backward, index->bwt_config.inverse_sa);
    free_bwt_index(index->forward_rev, index->forward, index->bwt_config.inverse_sa);
  }

  bwt_free_replace_table(&index->bwt_config);

  free(index->backward);
  free(index->forward);
  free(index->backward_rev);
  free(index->forward_rev);
  exome_free(index->karyotype);
  
  free(index);
}

//-----------------------------------------------------------------------------

void bwt_generate_index_files(char *ref_file, char *output_dir, 
			      unsigned int s_ratio, bool duplicate_strand, 
			      char *bases) {
  ref_vector X, Xi, B, Bi;
  vector C, C1;
  comp_matrix O, Oi;
  comp_vector S, R, Si, Ri;

  exome *ex = exome_new();
  bwt_config_t bwt_config;
  
  save_config(bases, duplicate_strand, output_dir);
  bwt_config.duplicate_strand = duplicate_strand;
  bwt_init_replace_table(&bwt_config, bases);

  // Calculating BWT
  encode_reference(&X, ex, ref_file, &bwt_config);
  save_ref_vector(&X, output_dir, "X");
  save_exome_file(ex, duplicate_strand, output_dir);
  print_vector(X.vector, X.n);

  calculate_and_save_B(&X, output_dir, "B");

  read_ref_vector(&B, output_dir, "B");
  print_vector(B.vector, B.n);

  calculate_C(&C, &C1, &B, bwt_config.nA);
  print_vector(C.vector, C.n);
  print_vector(C1.vector, C1.n);
  save_vector(&C, output_dir, "C");
  save_vector(&C1,output_dir, "C1");

  calculate_O(&O, &B, bwt_config.nA);
  print_comp_matrix(O);
  save_comp_matrix(&O, output_dir, "O");

  calculate_S_and_R(&S, &R, &B, &C, &O, s_ratio);
  print_vector(S.vector, S.n);
  print_vector(R.vector, S.n);
  save_comp_vector(&S, output_dir, "S");
  save_comp_vector(&R, output_dir, "R");

  free(B.vector);
  free_comp_matrix(NULL, &O);
  free(S.vector);
  free(R.vector);

  read_ref_vector(&Xi, output_dir, "X");
  revstring(Xi.vector, Xi.n);
  save_ref_vector(&Xi, output_dir, "Xi");
  print_vector(Xi.vector, Xi.n);

  calculate_and_save_B(&Xi, output_dir, "Bi");

  read_ref_vector(&Bi, output_dir, "Bi");
  print_vector(Bi.vector, Bi.n);

  calculate_O(&Oi, &Bi, bwt_config.nA);
  print_comp_matrix(Oi);
  save_comp_matrix(&Oi, output_dir, "Oi");

  calculate_S_and_R(&Si, &Ri, &Bi, &C, &Oi, s_ratio);
  print_vector(Si.vector, Si.n);
  print_vector(Ri.vector, Si.n);
  save_comp_vector(&Si, output_dir, "Si");
  save_comp_vector(&Ri, output_dir, "Ri");

  free(Bi.vector);
  free_comp_matrix(NULL, &Oi);
  free(Si.vector);
  free(Ri.vector);

  free(C.vector);
  free(C1.vector);
  exome_free(ex);
}

//-----------------------------------------------------------------------------
// general functions
//-----------------------------------------------------------------------------
/*
size_t bwt_map_seq(char *seq, 
		   bwt_optarg_t *bwt_optarg, 
		   bwt_index_t *index, 
		   array_list_t *mapping_list) {
  
  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_seq(seq, bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_seq(seq, bwt_optarg, index, mapping_list);  
}
*/
//-----------------------------------------------------------------------------
/*
size_t bwt_map_read(fastq_read_t *read, 
		    bwt_optarg_t *bwt_optarg, 
		    bwt_index_t *index, 
		    array_list_t *mapping_list) {
  
  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_read(read, bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_read(read, bwt_optarg, index, mapping_list);  
}
*/
//-----------------------------------------------------------------------------

/*size_t bwt_map_seqs(char **seqs, 
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
			      }*/

//-----------------------------------------------------------------------------
/*
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
}*/

//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------

size_t bwt_map_exact_seq(char *seq, 
			 bwt_optarg_t *bwt_optarg, 
			 bwt_index_t *index, 
			 array_list_t *mapping_list) {
  printf("\tIn function ...\n");
  unsigned int len = strlen(seq);
  int start = 0;
  int end = len - 1;
  uint8_t *code_seq = (uint8_t *) malloc(len * sizeof(uint8_t));
  result result;
  size_t l_aux, k_aux;
  unsigned int num_mappings = 0;
  char plusminus[2] = "-+";
  size_t idx, key, error, pos;
  alignment_t *alignment;
  char *cigar_p, *seq_dup, *seq_strand;
  unsigned int start_mapping;
  //char *chromosome = (char *)malloc(sizeof(char)*10);
  seq_strand = strdup(seq);

  encode_bases(code_seq, seq, len, index->bwt_config.table);
  //replaceBases(seq, code_seq, len);
  struct timeval t_start, t_end;
  char *optional_fields /*= (char *)malloc(sizeof(char)*256)*/;
  //printf("---> EXACT Search Read (%d): %s\n", len, seq);

  
  for (short int type = 1; type >= 0; type--) {
    result.k = 0;
    result.l = size_SA(index->backward) - 1;
    result.start = start;
    result.end = end;
    
    if (type == 1) {
      result.pos = end;
      //printf("Before call BWT\n");
      start_timer(t_start);
      BWExactSearchBackward(code_seq, index->backward, &result, NULL);
      stop_timer(t_start, t_end, time_bwt);
      //printf("After call BWT\n");
    } else {
      result.pos = start;
      start_timer(t_start);
      BWExactSearchForward(code_seq, index->forward_rev, &result, NULL);
      stop_timer(t_start, t_end, time_bwt);
      if(l_aux - k_aux + 1 < bwt_optarg->filter_read_mappings) {
	seq_reverse_complementary(seq_strand, len);
      }
      
    }
      
    start_timer(t_start);
    k_aux = result.k;
    l_aux = result.l;
    //printf("\tk=%d - l=%d\n", k_aux, l_aux);      
    if (l_aux - k_aux + 1 >= bwt_optarg->filter_read_mappings) {
      l_aux = k_aux + 1;
    }
    
    
    for (size_t j = k_aux; j <= l_aux; j++) {
      
      key = get_SA(j, index->backward);
      
      idx = binsearch(index->karyotype->offset, index->karyotype->size, key);
      //chromosome = index->karyotype->chromosome + (idx-1) * IDMAX;
      
      if(key + len <= index->karyotype->offset[idx]) {
	start_mapping = index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]);
	
	cigar_p = (char *)malloc(sizeof(char)*len);
	sprintf(cigar_p, "%d=%c", len, '\0');
	
	// save all into one alignment structure and insert to the list
	alignment = alignment_new();
	alignment_init_single_end(NULL, strdup(seq_strand), NULL, !type, 
				  idx - 1,				      
				  start_mapping, 
				  cigar_p, 1, 254, 1, (num_mappings > 0), 0, optional_fields, alignment);
	  
	if(!array_list_insert((void*) alignment, mapping_list)){
	  printf("Error to insert item into array list\n");
	}
	
	num_mappings++;
      }
    }
    
    stop_timer(t_start, t_end, time_search);
  }
  
  /*
  printf("Time BWT Search %.3fs\n", time_bwt/1000000);
  printf("Time S Search %.3fs\n", time_search/1000000);
  */
  free(seq_strand);
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
      if (mapping->seq_strand == 0) {
	memcpy(mapping->quality, read->quality, quality_len);
      }else {
	reverse_str(read->quality,
		    mapping->quality, quality_len);
      }
      mapping->quality[quality_len - 1] = '\0';
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
    individual_mapping_list_p[i] = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
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
      //printf("@DNA\n%s\n+\n%s\n", &(batch->seq[batch->data_indices[i]]), &(batch->seq[batch->data_indices[i]]));
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

size_t bwt_map_exact_seed(uint8_t *seq, size_t seq_len,
			  size_t seq_start, size_t seq_end,
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index,
			  array_list_t *mapping_list,
			  int id/*, int *limit_exceeded*/) {
  
  //printf("Process New Seeds\n");
  //printf("seq_start = %lu, seq_end = %lu\n", seq_start, seq_end);
  region_t *region;
  size_t start = 0;
  size_t end = seq_end - seq_start;
  result result;
  uint8_t *code_seq = &seq[seq_start];
  size_t len = seq_end - seq_start;
  size_t num_mappings = 0;
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
  size_t aux_seq_start, aux_seq_end;
  array_list_t *mappings = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, 
						     COLLECTION_MODE_ASYNCHRONIZED);
  int discard_seed = 0;
  int actual_mappings = 0;
  struct timeval t_start, t_end;
  //*limit_exceeded = 0;

  for (short int type = 1; type >= 0; type--) {
    result.k = 0;
    result.l = size_SA(index->backward) - 1;
    result.start = start;
    result.end = end;
    //printf("Search...\n");
    if (type == 1) {
      // strand +
      aux_seq_start = seq_start;
      aux_seq_end = seq_end;
      result.pos = end;	
      BWExactSearchBackward(code_seq, index->backward, &result, NULL);
    } else {
      // strand -
      aux_seq_start = seq_len - seq_end - 1;
      aux_seq_end = seq_len - seq_start - 1;
      result.pos = start;
      BWExactSearchForward(code_seq, index->backward_rev, &result, NULL);
    }
    //printf("Search end\n");
    //start_timer(t_start);

    //printf("results...\n");
    k_aux = result.k;
    l_aux = result.l;
    actual_mappings += (result.l - result.k + 1);

    if (actual_mappings > bwt_optarg->filter_seed_mappings) {
      //printf("Limit exceded LIMIT = %i, MAP = %i\n", bwt_optarg->filter_seed_mappings, actual_mappings);
      continue;
      //k_aux = result.k;
      //l_aux = result.k + 10; 
      //break;
    } else {
      k_aux = result.k;
      l_aux = result.l;
    }
      
    for (size_t j = k_aux; j <= l_aux; j++) {
      
      key = get_SA(j, index->backward);
      idx = binsearch(index->karyotype->offset, index->karyotype->size, key);     
      
      if (key + len <= index->karyotype->offset[idx]) {
	start_mapping = index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]) + 1;
	// save all into one alignment structure and insert to the list
	//printf("\t\t Region%i[%i:%lu-%lu]\n", !type, idx, start_mapping, start_mapping + len);
	region = region_bwt_new(idx, !type, start_mapping, start_mapping + len, aux_seq_start, aux_seq_end, seq_len, id);

	if (!array_list_insert((void*) region, mapping_list)){
	  printf("Error to insert item into array list\n");
	}
	
	num_mappings++;
      }
    }
  }

  array_list_free(mappings, NULL);

  return num_mappings;  
}

size_t bwt_map_exact_seed_by_region(char *seq, size_t seq_len,
				    size_t seq_start, size_t seq_end,
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index,
				    array_list_t *mapping_list,
				    int id, 
				    int strand_target,
				    int chromosome_target,
				    size_t start_target, 
				    size_t end_target) {  
  region_t *region;
  size_t start = 0;
  size_t end = seq_end - seq_start;
  result result;
  char *code_seq = &seq[seq_start];

  size_t len = seq_end - seq_start;
  size_t num_mappings = 0;
  char plusminus[2] = "-+";
  size_t idx, key, direction, error, pos;
  results_list *r_list;
  size_t l_aux, k_aux;
  alignment_t *alignment;

  size_t start_mapping;
  size_t aux_seq_start, aux_seq_end;
  array_list_t *mappings = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, 
						     COLLECTION_MODE_ASYNCHRONIZED);
  int discard_seed = 0;
  int actual_mappings = 0;
  struct timeval t_start, t_end;

  for (short int type = 1; type >= 0; type--) {
    if (!type != strand_target) { continue; }

    result.k = 0;
    result.l = size_SA(index->backward) - 1;
    result.start = start;
    result.end = end;
    if (type == 1) {
      // strand +
      aux_seq_start = seq_start;
      aux_seq_end = seq_end;
      result.pos = end;	
      start_timer(t_start);
      BWExactSearchBackward(code_seq, index->backward, &result, NULL);
      //BWExactSearchBackward(code_seq, start, end, &index->h_C, &index->h_C1, &index->h_O, result_p);
      stop_timer(t_start, t_end, time_bwt_seed);
    } else {
      // strand -
      aux_seq_start = seq_len - seq_end - 1;
      aux_seq_end = seq_len - seq_start - 1;

      result.pos = start;      
      start_timer(t_start);
      BWExactSearchForward(code_seq, index->backward_rev, &result, NULL);
      //BWExactSearchForward(code_seq, start, end, &index->h_rC, &index->h_rC1, &index->h_rO, result_p);
      stop_timer(t_start, t_end, time_bwt_seed);
    }

    //start_timer(t_start);

    k_aux = result.k;
    l_aux = result.l;
    actual_mappings += (result.l - result.k + 1);

    if (actual_mappings > 500) {
      k_aux = result.k;
      l_aux = result.k + 10;
    } else {
      k_aux = result.k;
      l_aux = result.l;
    }
      
    for (size_t j = k_aux; j <= l_aux; j++) {            
      key = get_SA(j, index->backward);
      idx = binsearch(index->karyotype->offset, index->karyotype->size, key);     
      if (idx != chromosome_target) { continue; }

      if (key + len <= index->karyotype->offset[idx]) {
	start_mapping = index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]) + 1;
	if (start_mapping < start_target || start_mapping > end_target) { continue; }
	// save all into one alignment structure and insert to the list
	printf(" \t::: %lu:%lu :::\n", idx, start_mapping);
	region = region_bwt_new(idx, !type, start_mapping, start_mapping + len, aux_seq_start, aux_seq_end, seq_len, id);
	
	if (!array_list_insert((void*) region, mapping_list)){
	  printf("Error to insert item into array list\n");
	}
	
	num_mappings++;
      }
    }
    //stop_timer(t_start, t_end, time_search_seed);
  }
  //  free(result_p);
  array_list_free(mappings, NULL);

  return num_mappings;  
}


//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seed(char *code_seq, size_t seq_len, 
			    size_t seq_start, size_t seq_end,
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    array_list_t *mapping_list,
			    int seed_id) {
  region_t *region; 
  size_t len = seq_end - seq_start + 1;
  const int max_mappings = 50;
  //  size_t len = strlen(seq);
  //  size_t start = 0;
  //  size_t end = len - 1;
  size_t start_mapping;
  //  char *code_seq = (char *) calloc(len, sizeof(char));
  //  replaceBases(seq, code_seq, len);

  // compare the vectors k and l to get mappings in the genome
  int aux_seq_start, aux_seq_end;
  unsigned int num_mappings = 0;
  char plusminus[2] = "-+";
  size_t idx, key, direction, error;
  //results_list *r_list;
  results_list r_list;
  result *r;
  alignment_t *alignment;
  size_t k_start, l_start;
  size_t len_calc;
  array_list_t *tmp_mapping_list = array_list_new(bwt_optarg->filter_read_mappings, 1.25f,
						  COLLECTION_MODE_ASYNCHRONIZED);
  new_results_list(&r_list, 20000);
  result res;

  // calculate vectors k and l
  intmax_t *k0 = (intmax_t *) malloc(len * sizeof(intmax_t));
  intmax_t *l0 = (intmax_t *) malloc(len * sizeof(intmax_t));
  intmax_t *k1 = (intmax_t *) malloc(len * sizeof(intmax_t));
  intmax_t *l1 = (intmax_t *) malloc(len * sizeof(intmax_t));
  intmax_t *ki0 = (intmax_t *) malloc(len * sizeof(intmax_t));
  intmax_t *li0 = (intmax_t *) malloc(len * sizeof(intmax_t));
  intmax_t *ki1 = (intmax_t *) malloc(len * sizeof(intmax_t));
  intmax_t *li1 = (intmax_t *) malloc(len * sizeof(intmax_t));
  
  int back0_nt, forw0_nt, back1_nt, forw1_nt;
  
  int pos, final_pos, final_err;

  BWExactSearchVectorBackward(code_seq, seq_start, seq_end, 0, size_SA(index->backward)-1,
			      k1, l1, index->backward);
  
  BWExactSearchVectorForward(code_seq, seq_start, seq_end, 0, size_SA(index->forward)-1,
			     ki1, li1, index->forward);
  
  BWExactSearchVectorForward(code_seq, seq_start, seq_end, 0, size_SA(index->backward_rev)-1,
			     k0, l0, index->backward_rev);
  
  BWExactSearchVectorBackward(code_seq, seq_start, seq_end, 0, size_SA(index->forward_rev)-1,
			      ki0, li0, index->forward_rev);
  
  for (int type = 1; type >= 0; type--) {
    //    r_list = (results_list *) calloc(1, sizeof(results_list));
    //new_results_list(r_list, 2000);
    r_list.num_results = 0;
    //r_list->n = 0;
    r_list.read_index = 0;

    if (type == 1) {
      BWSearch1VectorHelper(code_seq, seq_start, seq_end, k1, l1, ki1, li1, 
			    index->backward, index->forward, 
			    &r_list, index->bwt_config.nA);
      //ALERT: This not works!!!
      //bound_result(&res, seq_start, seq_end);
      //change_result(&res, 0, size_SA(index->backward)-1, 0);
      //BWSearch1CPU(code_seq, index->backward, index->forward, &res, &r_list, NULL, index->bwt_config.nA);
      aux_seq_start = seq_start;
      aux_seq_end = seq_end;
    } else {
      BWSearch1VectorHelper(code_seq, seq_start, seq_end, ki0, li0, k0, l0, 
			    index->forward_rev, index->backward_rev, 
			    &r_list, index->bwt_config.nA); 
      //ALERT: This not works!!!
      //bound_result(&res, seq_start, seq_end);
      //change_result(&res, 0, size_SA(index->backward)-1, 0);
      //BWSearch1CPU(code_seq, index->forward_rev, 
      //		   index->backward_rev, &res, &r_list, NULL, index->bwt_config.nA);

      aux_seq_start = seq_len - seq_end - 1;
      aux_seq_end = seq_len - seq_start - 1;
    }

    //printf("Num results: %i\n", r_list.num_results);
    for (size_t ii = 0; ii < r_list.num_results; ii++) {
      //for (size_t ii = 0; ii < r_list->n; ii++) {
      r = &r_list.list[ii];

      if(!r->num_mismatches)
	error = 0;
      else
	error = r->err_kind[0];

      len_calc = len;
      if (error == DELETION) {
	len_calc--;
      } else if (error == INSERTION) {
	len_calc++;
      }

      //Bug, correct the cigar I & D
      pos = r->err_pos[0];      
      if (error == DELETION || pos == 0 || pos == len -1) {
	if (type) {
	  final_pos = pos;
	} else {
	  final_pos = len - 1 - pos;
	}
      } else {
	if (error == INSERTION) {
	  if (type) {
	    if (r->dir) { final_pos = pos; }
	    else { final_pos = pos + 1; }
	  } else {
	    if (r->dir) { final_pos = len - pos; }
	    else { final_pos = len - pos - 1; }	    
	  }
	}
      }

      if (error == MISMATCH) {
	final_err = -1;
      } else {
	if (error == DELETION) {
	  final_err = SEED_INSERTION;
	} else {
	  final_err = SEED_DELETION;
	}
      }
      final_err += seq_start;
      //End

      if ((r->l - r->k + 1) < max_mappings) {//bwt_optarg->filter_read_mappings) {
	k_start = r->k;
	l_start = r->l;
      }else  {
	//printf("LIMIT EXCEEDED %lu\n", (r->l - r->k + 1));
	k_start = r->k;
	l_start = r->k + 1;
      }

      if (type) {
	direction = r->dir;
      } else {
	direction = !r->dir;
      }
      
      //k_start = r->k;
      //l_start = r->l;      

      //printf("%lu - %lu\n", k_start, l_start);
      for (size_t j = k_start; j <= l_start; j++) {

	key = (direction)
	  ? size_SA(index->forward) - get_SA(j, index->forward) - len_calc - 1
	  : get_SA(j, index->backward);

	idx = binsearch(index->karyotype->offset, index->karyotype->size, key);
	if(key + len_calc <= index->karyotype->offset[idx]) {
	    //printf("\tvalue idx=%d\n", idx);
	  /*printf("\t%s\t%c\t%s %u %s error: %s, pos: %i, base: %i ",
		 "nothing", plusminus[type],
		 index->karyotype->chromosome + (idx-1) * IDMAX,
		 index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]),
		 seq, bwt_error_type(r->err_kind[0]), r->position[0], r->base[0]);
	  */
	  start_mapping = index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]) + 1;
	  // save all into one alignment structure and insert to the list
	  region = region_bwt_new(idx, !type, start_mapping, start_mapping + len - 1, aux_seq_start, aux_seq_end, seq_len, seed_id);
	  
	  region->pos_err = final_pos;
	  region->type_err = final_err;

	  //printf("Report Region in position [%i:%lu-%lu]\n", idx, start_mapping, start_mapping + len - 1);
	  //printf("BWT-ERR: %i, %i\n", region->pos_err, region->type_err);
	    
	  if(!array_list_insert((void*) region, tmp_mapping_list)){
	    printf("Error to insert item into array list\n");
	  }
	  num_mappings++;
	}
      }
    }

    //
    //
    //Search for equal BWT mappings and set the mappings that will be delete
    int n_mappings = array_list_size(tmp_mapping_list);
    region_t *reg_1, *reg_2;
    unsigned int *delete_mark = (unsigned int *)calloc(n_mappings, sizeof(unsigned int));
    const int max_distance = 10;
    //printf("------------------Num mappings %i---------------\n", n_mappings);     
    for (int a1 = n_mappings - 1; a1 >= 1; a1--) {
      if (!delete_mark[a1]) {
	reg_1 = array_list_get(a1, tmp_mapping_list);
	//printf("Alig1(%i): %i chromosome, %i seq_strand, %s cigar:\n", a1, alig_1->chromosome, alig_1->seq_strand, alig_1->cigar);
	for (int a2 = a1 - 1; a2 >= 0; a2--) {
	  reg_2 = array_list_get(a2, tmp_mapping_list);
	  size_t dist = abs(reg_1->start - reg_2->start);
	  //printf("\t Alig2(%i): %i chromosome, %i seq_strand, %i dist, %i delete mark, %s cigar\n", a2, alig_2->chromosome, alig_2->seq_strand, dist, delete_mark[a2],alig_2->cigar );
	  if (reg_1->chromosome_id == reg_2->chromosome_id &&
	      dist < max_distance && 
	      !delete_mark[a2]) {
	    //Same chromosome && same position
	    delete_mark[a2] = 1;
	  }
	}
      }
    }
    //Delete all set mappings
    for (int m = n_mappings - 1; m >= 0; m--) {
      reg_1 = array_list_remove_at(m, tmp_mapping_list);
      if (delete_mark[m]) {
	region_bwt_free(reg_1);
      } else {
	array_list_insert((void*) reg_1, mapping_list);
	num_mappings++;
      }
    }
    free(delete_mark);
  } // end for type 

  //free(code_seq);
  free(r_list.list);
  free(k0);
  free(l0);
  free(k1);
  free(l1);
  free(ki0);
  free(li0);
  free(ki1);
  free(li1);

  array_list_free(tmp_mapping_list, NULL);

  return num_mappings;

}

size_t bwt_map_inexact_seed_by_region(char *code_seq, size_t seq_len, 
				      size_t seq_start, size_t seq_end,
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list,
				      int seed_id, int chromosome_id,
				      size_t start_pos, 
				      size_t end_pos, int strand_target) {  
  region_t *region; 

  size_t len = seq_end - seq_start + 1;

  const int max_mappings = 50;
  //  size_t len = strlen(seq);
  //  size_t start = 0;
  //  size_t end = len - 1;
  size_t start_mapping;
  //  char *code_seq = (char *) calloc(len, sizeof(char));
  //  replaceBases(seq, code_seq, len);

  // compare the vectors k and l to get mappings in the genome
  int aux_seq_start, aux_seq_end;
  unsigned int num_mappings = 0;
  char plusminus[2] = "-+";
  size_t idx, key, direction, error, pos;
  //results_list *r_list;
  results_list r_list;
  result *r;
  alignment_t *alignment;
  size_t k_start, l_start;
  size_t len_calc;
  array_list_t *tmp_mapping_list = array_list_new(bwt_optarg->filter_read_mappings, 1.25f,
						  COLLECTION_MODE_ASYNCHRONIZED);
  new_results_list(&r_list, 20000);
  result res;

  for (int type = 1; type >= 0; type--) {
    if (!type != strand_target) { continue; }
    //    r_list = (results_list *) calloc(1, sizeof(results_list));
    //new_results_list(r_list, 2000);
    r_list.num_results = 0;
    //r_list->n = 0;
    r_list.read_index = 0;

     if (type == 1) {
      bound_result(&res, seq_start, seq_end);
      change_result(&res, 0, size_SA(index->backward)-1, 0);
      BWSearch1CPU(code_seq, index->backward, index->forward, &res, &r_list, NULL, index->bwt_config.nA);

      aux_seq_start = seq_start;
      aux_seq_end = seq_end;
    } else {
      bound_result(&res, seq_start, seq_end);
      change_result(&res, 0, size_SA(index->backward)-1, 0);
      BWSearch1CPU(code_seq, index->forward_rev, 
		   index->backward_rev, &res, &r_list, NULL, index->bwt_config.nA);

      aux_seq_start = seq_len - seq_end - 1;
      aux_seq_end = seq_len - seq_start - 1;
    }

    for (size_t ii = 0; ii < r_list.num_results; ii++) {
      //for (size_t ii = 0; ii < r_list->n; ii++) {
      r = &r_list.list[ii];

      if(!r->num_mismatches)
	error = 0;
      else
	error = r->err_kind[0];

      len_calc = len;
      if (error == DELETION) {
	len_calc--;
      } else if (error == INSERTION) {
	len_calc++;
      }

      if ((r->l - r->k + 1) < max_mappings) {//bwt_optarg->filter_read_mappings) {
	k_start = r->k;
	l_start = r->l;
      }else  {
	k_start = r->k;
	l_start = r->k + 1;
      }

      if (type) {
	direction = r->dir;
      } else {
	direction = !r->dir;
      }
      
      //k_start = r->k;
      //l_start = r->l;      

      for (size_t j = k_start; j <= l_start; j++) {

	key = (direction)
	  ? size_SA(index->forward) - get_SA(j, index->forward) - len_calc - 1
	  : get_SA(j, index->backward);

	idx = binsearch(index->karyotype->offset, index->karyotype->size, key);
	if (idx != chromosome_id) { continue; }

	if(key + len_calc <= index->karyotype->offset[idx]) {
	    //printf("\tvalue idx=%d\n", idx);
	  /*printf("\t%s\t%c\t%s %u %s error: %s, pos: %i, base: %i ",
		 "nothing", plusminus[type],
		 index->karyotype->chromosome + (idx-1) * IDMAX,
		 index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]),
		 seq, bwt_error_type(r->err_kind[0]), r->position[0], r->base[0]);
	  */
	  start_mapping = index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]);

	  //printf("(%i)Report Region in position [%i:%lu|%i-%i|%lu]\n", !type, idx, start_mapping, 
	  //	 aux_seq_start, aux_seq_end, start_mapping + end);
	  if (start_mapping < start_pos || start_mapping > end_pos ) { continue; }
	  // save all into one alignment structure and insert to the list
	  //printf("%i:%lu-%lu , %i-%i\n", idx, start_mapping, start_mapping + len, aux_seq_start, aux_seq_end);	  
	  region = region_bwt_new(idx, !type, start_mapping, start_mapping + len-1, aux_seq_start, aux_seq_end, seq_len, seed_id);
	    
	  if(!array_list_insert((void*) region, tmp_mapping_list)){
	    printf("Error to insert item into array list\n");
	  }
	  num_mappings++;
	}
      }
    }

    //
    //
    //Search for equal BWT mappings and set the mappings that will be delete
    int n_mappings = array_list_size(tmp_mapping_list);
    region_t *reg_1, *reg_2;
    unsigned int *delete_mark = (unsigned int *)calloc(n_mappings, sizeof(unsigned int));
    const int max_distance = 10;
    //printf("------------------Num mappings %i---------------\n", n_mappings);     
    for (int a1 = n_mappings - 1; a1 >= 1; a1--) {
      if (!delete_mark[a1]) {
	reg_1 = array_list_get(a1, tmp_mapping_list);
	//printf("Alig1(%i): %i chromosome, %i seq_strand, %s cigar:\n", a1, alig_1->chromosome, alig_1->seq_strand, alig_1->cigar);
	for (int a2 = a1 - 1; a2 >= 0; a2--) {
	  reg_2 = array_list_get(a2, tmp_mapping_list);
	  size_t dist = abs(reg_1->start - reg_2->start);
	  //printf("\t Alig2(%i): %i chromosome, %i seq_strand, %i dist, %i delete mark, %s cigar\n", a2, alig_2->chromosome, alig_2->seq_strand, dist, delete_mark[a2],alig_2->cigar );
	  if (reg_1->chromosome_id == reg_2->chromosome_id &&
	      dist < max_distance && 
	      !delete_mark[a2]) {
	    //Same chromosome && same position
	    delete_mark[a2] = 1;
	  }
	}
      }
    }
    //Delete all set mappings
    for (int m = n_mappings - 1; m >= 0; m--) {
      reg_1 = array_list_remove_at(m, tmp_mapping_list);
      if (delete_mark[m]) {
	region_bwt_free(reg_1);
      } else {
	array_list_insert((void*) reg_1, mapping_list);
	num_mappings++;
      }
    }
    free(delete_mark);
  } // end for type 

  //free(code_seq);
  free(r_list.list);
  array_list_free(tmp_mapping_list, NULL);

  return num_mappings;

}

//-----------------------------------------------------------------------------

void *__bwt_generate_anchor_list(size_t k_start, size_t l_start, int len_calc, bwt_optarg_t *bwt_optarg, 
				 bwt_index_t *index, int type, array_list_t *anchor_list, int type_anchor,
				 int dsp) {

     size_t idx, key, direction;
     size_t start_mapping;
     const int MAX_BWT_ANCHORS = 50;
     bwt_anchor_t *bwt_anchor;
     
     //len_calc = 100;

     //printf("dsp:%i, nt:%i, len_calc:%i\n", dsp, len_calc, len_calc + dsp);

     if (type) { //(+)
       if (type_anchor == BACKWARD_ANCHOR) {
	 //printf("(+) BACKWARD %i\n", len_calc);
	 direction = BACKWARD_ANCHOR;
       } else {
	 //printf("(+) FORWARD %i\n", len_calc);
	 direction = FORWARD_ANCHOR; 
       }
     } else { //(-)
       if (type_anchor == BACKWARD_ANCHOR) {
	 //printf("(-) BACKWARD %i\n", len_calc);
	 direction = BACKWARD_ANCHOR; 
	 //len_calc += dsp;
       } else {
	 //printf("(-) FORWARD %i\n", len_calc);
	 direction = FORWARD_ANCHOR;
	 //len_calc += dsp;
       }
     }
  
     //printf("k_start = %lu, l_start = %lu < %lu\n", k_start, l_start, MAX_BWT_ANCHORS);
     if (l_start - k_start < MAX_BWT_ANCHORS) {
       for (size_t j = k_start; j <= l_start; j++) {

	 key = (direction)
	   ? size_SA(index->forward) - get_SA(j, index->forward) - len_calc - 1
	   : get_SA(j, index->backward);

	 idx = binsearch(index->karyotype->offset, index->karyotype->size, key);
	 if(key + len_calc <= index->karyotype->offset[idx]) {
	   start_mapping = index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]);	
	   bwt_anchor = bwt_anchor_new(!type, idx - 1, start_mapping + 1, start_mapping + len_calc, type_anchor);
	   //printf("anchor: %i:%lu-%lu\n", idx, bwt_anchor->start, bwt_anchor->end);
	   array_list_insert(bwt_anchor, anchor_list);
	   //printf("\tTest k-l: %i:%lu\n", idx, start_mapping);
	 }
       }
     }

}

//-----------------------------------------------------------------------------

size_t bwt_map_inexact_read(fastq_read_t *read, 
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    array_list_t *mapping_list) {

  return __bwt_map_inexact_read(read, 
			       bwt_optarg, 
			       index, 
			       mapping_list);  
}

size_t __bwt_map_inexact_read(fastq_read_t *read, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list) {

  //return array_list_size(mapping_list);
   
     char *seq = read->sequence;
     alignment_t *alignment;
     size_t len = strlen(seq);

     if (len < 5) {
	  char aux[len + 2];
	  sprintf(aux, "%luX", len);

	  char *quality_clipping = (char *)malloc(sizeof(char)*50);
	  sprintf(quality_clipping, "%i", NONE_HARD_CLIPPING);

	  alignment = alignment_new();
	  alignment_init_single_end(NULL,
				    strdup(seq),
				    quality_clipping,
				    0,
				    -1,
				    -1,
				    strdup(aux), 1, 0, 0, 0, 0, NULL, alignment);
	  array_list_insert((void*) alignment, mapping_list);
	  return 1;
     }
 
     char *seq_dup, *seq_strand;
     size_t start = 0;
     size_t end = len - 1;
     size_t len_calc = len;

     uint8_t *code_seq = (uint8_t *) malloc(len * sizeof(uint8_t));
     encode_bases(code_seq, seq, len, index->bwt_config.table);

     // calculate vectors k and l
     intmax_t *k0 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *l0 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *k1 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *l1 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *ki0 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *li0 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *ki1 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *li1 = (intmax_t *) malloc(len * sizeof(intmax_t));

     int back0_nt, forw0_nt, back1_nt, forw1_nt;

     back1_nt = BWExactSearchVectorBackward(code_seq, start, end, 0, size_SA(index->backward)-1,
					    k1, l1, index->backward);
     
     forw1_nt = BWExactSearchVectorForward(code_seq, start, end, 0, size_SA(index->forward)-1,
					   ki1, li1, index->forward);
  
     back0_nt = BWExactSearchVectorForward(code_seq, start, end, 0, size_SA(index->backward_rev)-1,
					   k0, l0, index->backward_rev);
  
     forw0_nt = BWExactSearchVectorBackward(code_seq, start, end, 0, size_SA(index->forward_rev)-1,
					    ki0, li0, index->forward_rev);

     //compare the vectors k and l to get mappings in the genome
     size_t num_mappings = 0;
     char plusminus[2] = "-+";
     size_t idx, key, direction;
     char error;
     int pos;
     //results_list *r_list;
     results_list r_list;
     result *r;
     char error_debug = 0;
     char *cigar_dup;
     char cigar[1024];
     size_t cigar_len, num_cigar_ops;
     array_list_t *mapping_list_filter;
     alignment_t *best_alignment, *aux_alignment;
     size_t best_pos, array_size;
     //size_t i, j, z;
     size_t *allocate_pos_alignments;
     size_t k_start, l_start;

     size_t start_mapping;
     size_t tot_alignments = 0;
     const int MAX_BWT_ALIGNMENTS = 10;
     int filter_exceeded = 0;
     //seq_dup = (char *)malloc(sizeof(char)*(len + 1));
     seq_strand = seq;
     error = MISMATCH;


     char *quality_clipping = (char *) malloc(sizeof(char) * 50);
     seq_dup = (char *) malloc(sizeof(char) * (len + 1));

     new_results_list(&r_list, bwt_optarg->filter_read_mappings);

     array_list_t *tmp_mapping_list = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, 
						     COLLECTION_MODE_ASYNCHRONIZED);

     for (int type = 1; type >= 0; type--) {
       //printf("Process Strand(%i)\n", !type);
	  r_list.num_results = 0;
	  r_list.read_index = 0;

	  //    printf("*** bwt.c: calling BWSearch1 with type = %d...\n", type);
	  if (type == 1) {
	       BWSearch1VectorHelper(code_seq, start, end, k1, l1, ki1, li1, 
				     index->backward, index->forward, 
				     &r_list, index->bwt_config.nA);
	  } else {
	       BWSearch1VectorHelper(code_seq, start, end, ki0, li0, k0, l0, 
			 index->forward_rev, index->backward_rev, 
				     &r_list, index->bwt_config.nA);      
	       if (r_list.num_results) {
		 seq_strand = read->revcomp;
		 //seq_reverse_complementary(seq_strand, len);
	       }
	  }

	  //printf("*** bwt.c: calling BWSearch1 with type = %d (num_results = %d). Done !!\n", type, r_list.num_results);

	  for (size_t ii = 0; ii < r_list.num_results; ii++) {
	    r = &r_list.list[ii];

	       //printf("Errors number %d\n", r->num_mismatches);
	       //if(r == NULL){printf("ERROR: Result list position null");exit(-1);}
	       if(!r->num_mismatches)
		    error = 0;
	       else
		    error = r->err_kind[0];

	       pos = r->err_pos[0];
	       //pos = r->position[0];
	  

	       if (type) {
		    direction = r->dir;
	       } else {
		    direction = !r->dir;
	       }

	       len_calc = len;
	       if (error == DELETION) {
		    len_calc--;
	       } else if (error == INSERTION) {
		    len_calc++;
	       }

	       if (error != 0) { continue; }

	       // generating cigar
	       sprintf(quality_clipping, "%i", NONE_HARD_CLIPPING);
	       if (error == 0) {
		 sprintf(cigar, "%luM%c", len, '\0');
		 num_cigar_ops = 1;
		 memcpy(seq_dup, seq_strand, len);
		 seq_dup[len] = '\0';
		 
	       } else if (error == MISMATCH) {
		 if (pos == 0) {
		   //Positive strand
		   if(type) { 
		     sprintf(cigar, "1S%luM%c", len-1, '\0'); 
		     start_mapping++;
		   }
		   else { 
		     sprintf(cigar, "%luM1S%c", len-1, '\0'); 
		   }
		   num_cigar_ops = 2;
		 } else if (pos == len - 1) {
		   //Positive strand
		   if(type) { 
		     sprintf(cigar, "%luM1S%c", len - 1, '\0'); 
		   }
		   else{ 
		     sprintf(cigar, "1S%luM%c", len-1, '\0'); 
		     start_mapping++;
		   }
		   num_cigar_ops = 2;
		 } else {
		   sprintf(cigar, "%luM%c", len, '\0');
		   num_cigar_ops = 1;
		 }
		 memcpy(seq_dup, seq_strand, len);
		 seq_dup[len] = '\0';
		 //printf("MISMATCH\n");
		 
	       } else if (error == INSERTION) {
		 //printf("INSERTION\n");
		 if (pos == 0) {
		   if(type) {
		     sprintf(cigar, "1M1D%luM%c", len - 1, '\0'); 
		   }
		   else{ 
		     sprintf(cigar, "%luM1D1M%c", len - 1, '\0'); 
		   }	      
		 } else if (pos == len - 1) {
		   if(type) { 
		     sprintf(cigar, "%luM1D1M%c", len - 1, '\0'); 
		   }
		   else{ 
		     sprintf(cigar, "1M1D%luM%c", len - 1, '\0'); 
		   }
		 } else {		   
		   if(type) {
		     if(r->dir)
		       sprintf(cigar, "%iM1D%luM%c", pos, len - pos, '\0');
		     else
		       sprintf(cigar, "%iM1D%luM%c", pos + 1, len - pos - 1, '\0');
		   } else { 
		     if(r->dir)
		       sprintf(cigar, "%luM1D%dM%c", len - pos, pos, '\0');
		     else
		       sprintf(cigar, "%luM1D%dM%c", len - pos - 1, pos + 1, '\0');
		   }
		 }
		 num_cigar_ops = 3;
		 memcpy(seq_dup, seq_strand, len);
		 seq_dup[len] = '\0';
	       } else if (error == DELETION) {	     
		 //printf("DELETION\n");
		 if (pos == 0) {
		   if(type) { sprintf(cigar, "1I%luM%c", len -1, '\0'); }
		   else{ 
		     sprintf(cigar, "%luM1I%c", len -1, '\0'); 
		     //		   start_mapping++;
		   }
		   
		   num_cigar_ops = 2;		
		 } else if (pos == len - 1) {
		   if(type) { 
		     sprintf(cigar, "%luM1I%c", len -1, '\0'); 
		     //		   start_mapping++;
		   }
		   else{ 
		     sprintf(cigar, "1I%luM%c", len -1, '\0'); 
		   }
		   num_cigar_ops = 2;
		 } else {
		   if(type) { 
		     sprintf(cigar, "%dM1I%luM%c", pos, len - pos - 1, '\0'); 
		   } else { 
		     sprintf(cigar, "%luM1I%dM%c", len - pos - 1, pos, '\0'); 
		   }
		   num_cigar_ops = 3;
		 }
		 memcpy(seq_dup, seq_strand , len );
		 seq_dup[len] = '\0';
		 
	       }else{
		 printf("NUM MAPPINGS %lu -> POS %d -> ERROR %d -> (%lu):%s", num_mappings, pos, error, len, seq);
		 continue;
		 //exit(-1);
		 //error_debug = 1;
	       }
	       k_start = r->k;
	       l_start = r->l;
	       //}

	       //printf("%lu-%lu :: num_results %lu\n", k_start, l_start, r_list.num_results);
	       tot_alignments += (l_start - k_start);
	       if (tot_alignments >  bwt_optarg->filter_read_mappings) {
		 filter_exceeded = 1;
		 LOG_DEBUG_F("Filter exceeded: num. read mappings: %i (total: %i > filter %i)\n", 
			     l_start - k_start, tot_alignments, bwt_optarg->filter_read_mappings);
		 break;
	       }

	       for (size_t j = k_start; j <= l_start; j++) {
		 key = (direction)
		   ? size_SA(index->forward) - get_SA(j, index->forward) - len_calc - 1
		   : get_SA(j, index->backward);

		 idx = binsearch(index->karyotype->offset, index->karyotype->size, key);
		 if(key + len_calc <= index->karyotype->offset[idx]) {
		   //start_mapping = index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]) + 1;
		   start_mapping = index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]); //¿+/-1?
		   // save all into one alignment structure and insert to the list

		   //printf("*****Alignments %i:%lu\n", idx, start_mapping);
		   alignment = alignment_new();
		   alignment_init_single_end(strdup(read->id), strdup(seq_dup), strdup(quality_clipping), !type, 
					     idx - 1, //index->karyotype->chromosome + (idx-1) * IDMAX,
					     start_mapping, //index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]), 
					     strdup(cigar), num_cigar_ops, error == 0 ? 0 : 1, 1, (num_mappings > 0), 0, NULL, alignment);

		   array_list_insert((void*) alignment, tmp_mapping_list);
		 }
	       }//end for k and l
	  }//end for 
	  
	  if (filter_exceeded) {
	    array_list_clear(tmp_mapping_list, (void *)alignment_free);
	    array_list_set_flag(2, mapping_list);
	    break;
	  }

	  //
	  //
	  //Search for equal BWT mappings and set the mappings that will be delete
	  int n_mappings = array_list_size(tmp_mapping_list);
	  alignment_t *alig_1, *alig_2;
	  unsigned int *delete_mark = (unsigned int *)calloc(n_mappings, sizeof(unsigned int));
	  const int max_distance = 10;
	  //printf("------------------Num mappings %i---------------\n", n_mappings);
     
	  for (int a1 = n_mappings - 1; a1 >= 1; a1--) {
	    if (!delete_mark[a1]) {
	      alig_1 = array_list_get(a1, tmp_mapping_list);
	      //printf("Alig1(%i): %i chromosome, %i seq_strand, %s cigar:\n", a1, alig_1->chromosome, alig_1->seq_strand, alig_1->cigar);
	      for (int a2 = a1 - 1; a2 >= 0; a2--) {
		alig_2 = array_list_get(a2, tmp_mapping_list);
		size_t dist = abs(alig_1->position - alig_2->position);
		//printf("\t Alig2(%i): %i chromosome, %i seq_strand, %i dist, %i delete mark, %s cigar\n", a2, alig_2->chromosome, alig_2->seq_strand, dist, delete_mark[a2],alig_2->cigar );
		if (alig_1->chromosome == alig_2->chromosome &&
		    dist < max_distance && 
		    !delete_mark[a2]) {
		  //Same chromosome && same position
		  if (alig_1->num_cigar_operations < alig_2->num_cigar_operations) {
		    //printf("\tSet read %i\n", a2);
		    delete_mark[a2] = 1;
		  } else {
		    //printf("\tSet read %i\n", a2);
		    delete_mark[a1] = 1;
		  }
		}
	      }
	    }
	  }
	  //Delete all set mappings
	  int primary_delete = 0;
	  for (int m = n_mappings - 1; m >= 0; m--) {
	    alig_1 = array_list_remove_at(m, tmp_mapping_list);
	    if (delete_mark[m]) {
	      if (!is_secondary_alignment(alig_1)) { primary_delete = 1; }
	      alignment_free(alig_1);
	    } else {
	      array_list_insert((void*) alig_1, mapping_list);
	      set_secondary_alignment(num_mappings > 0, alig_1);
	      num_mappings++;
	    }
	  }
	  if (primary_delete) { 
	    alig_1 = array_list_get(0, mapping_list); 
	    set_secondary_alignment(0, alig_1);
	  }
	  free(delete_mark);	  	  
     } // end for type 

     if (array_list_size(mapping_list) == 0) {
       if (array_list_get_flag(mapping_list) == 1) {

	 //====================================================================================
	 //printf("BWT(+): FORWARD(k-l)%i:%lu-%lu BACKWARD(k-l)%i:%lu-%lu\n", forw1_nt, ki1[forw1_nt], li1[forw1_nt], back1_nt, k1[back1_nt], l1[back1_nt]);
	 //printf("BWT(-): FORWARD(k-l)%i:%lu-%lu BACKWARD(k-l)%i:%lu-%lu\n", forw0_nt, ki0[forw0_nt], li0[forw0_nt], back0_nt, k0[back0_nt], l0[back0_nt]);
	 array_list_t *forward_anchor_list, *backward_anchor_list;       
	 int type = 1;//(+)

	 //printf("BACKWARD (+)\n");
	 __bwt_generate_anchor_list(k1[end - back1_nt + 1], 
				    l1[end - back1_nt + 1], back1_nt, bwt_optarg, 
				    index, type, mapping_list, BACKWARD_ANCHOR, 0);
	 //printf("FORWARD (+)\n");
	 __bwt_generate_anchor_list(ki1[start + forw1_nt - 1], 
				    li1[start + forw1_nt - 1], forw1_nt, bwt_optarg, 
				    index, type,  mapping_list, FORWARD_ANCHOR, 0);
	 type = 0;//(-)
	 //printf("BACKWARD (-)\n");
	 /* __bwt_generate_anchor_list(last_k0, last_l0, back0_nt, bwt_optarg, 
				    index, type,  mapping_list, BACKWARD_ANCHOR, back0_nt - forw0_nt);
	 //printf("FORWARD (-)\n");
	 __bwt_generate_anchor_list(last_ki0, last_li0, forw0_nt, bwt_optarg, 
	 index, type,  mapping_list, FORWARD_ANCHOR, back0_nt - forw0_nt);*/

	 __bwt_generate_anchor_list(k0[start + back0_nt - 1], 
				    l0[start + back0_nt - 1], back0_nt, bwt_optarg, 
				    index, type,  mapping_list, BACKWARD_ANCHOR, forw0_nt - back0_nt);
	 //printf("FORWARD (-)\n");
	 __bwt_generate_anchor_list(ki0[end - forw0_nt + 1], 
				    li0[end - forw0_nt + 1], forw0_nt, bwt_optarg, 
				    index, type,  mapping_list, FORWARD_ANCHOR, forw0_nt - back0_nt);

       }
     } else {
       int header_len;
       size_t num_mapping = array_list_size(mapping_list);
       //printf("BWT:Report Total mappings %i\n", num_mappings);
       for (int i = 0; i < num_mapping; i++) {
	 alignment = array_list_get(i, mapping_list);
	 //header_len = strlen(read->id);
	 //alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	 //get_to_first_blank(read->id, header_len, alignment->query_name);
	 bwt_cigar_cpy(alignment, read->quality);
	 //alignment->quality = strdup(&(batch->quality[batch->data_indices[i]]));                   
	 //************************* OPTIONAL FIELDS ***************************//
	 alignment = add_optional_fields(alignment, num_mappings, len);
       }
       array_list_set_flag(0, mapping_list);
       //printf("########## EXACT! #########\n");
     }

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

     free(seq_dup);
     free(quality_clipping);
     array_list_free(tmp_mapping_list, NULL);

     return array_list_size(mapping_list);

}

//-----------------------------------------------------------------------------

size_t bwt_map_inexact_read_2(fastq_read_t *read, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list) {
  
     char *seq = read->sequence;
     alignment_t *alignment;
     size_t len = strlen(seq);

     if (len < 5) {
	  char aux[len + 2];
	  sprintf(aux, "%luX", len);

	  char *quality_clipping = (char *)malloc(sizeof(char)*50);
	  sprintf(quality_clipping, "%i", NONE_HARD_CLIPPING);

	  alignment = alignment_new();
	  alignment_init_single_end(NULL,
				    strdup(seq),
				    quality_clipping,
				    0,
				    -1,
				    -1,
				    strdup(aux), 1, 0, 0, 0, 0, NULL, alignment);
	  array_list_insert((void*) alignment, mapping_list);
	  return 1;
     }
 
     char *seq_dup, *seq_strand;
     size_t start = 0;
     size_t end = len - 1;
     size_t len_calc = len;

     uint8_t *code_seq = (uint8_t *) malloc(len * sizeof(uint8_t));
     encode_bases(code_seq, seq, len, index->bwt_config.table);

     // calculate vectors k and l
     /*intmax_t *k0 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *l0 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *k1 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *l1 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *ki0 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *li0 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *ki1 = (intmax_t *) malloc(len * sizeof(intmax_t));
     intmax_t *li1 = (intmax_t *) malloc(len * sizeof(intmax_t));
     */
     int back0_nt, forw0_nt, back1_nt, forw1_nt;
     /*
     back1_nt = BWExactSearchVectorBackward(code_seq, start, end, 0, size_SA(index->backward)-1,
					    k1, l1, index->backward);
     
     forw1_nt = BWExactSearchVectorForward(code_seq, start, end, 0, size_SA(index->forward)-1,
					   ki1, li1, index->forward);
  
     back0_nt = BWExactSearchVectorForward(code_seq, start, end, 0, size_SA(index->backward_rev)-1,
					   k0, l0, index->backward_rev);
  
     forw0_nt = BWExactSearchVectorBackward(code_seq, start, end, 0, size_SA(index->forward_rev)-1,
					    ki0, li0, index->forward_rev);
     */
     //compare the vectors k and l to get mappings in the genome
     size_t num_mappings = 0;
     char plusminus[2] = "-+";
     size_t idx, key, direction;
     char error;
     int pos;
     //results_list *r_list;
     result *r;
     char error_debug = 0;
     char *cigar_dup;
     char cigar[1024];
     size_t cigar_len, num_cigar_ops;
     array_list_t *mapping_list_filter;
     alignment_t *best_alignment, *aux_alignment;
     size_t best_pos, array_size;
     //size_t i, j, z;
     size_t *allocate_pos_alignments;
     size_t k_start, l_start;

     size_t start_mapping;
     size_t tot_alignments = 0;
     const int MAX_BWT_ALIGNMENTS = 10;
     int filter_exceeded = 0;
     results_list rl_prev, rl_next, rl_prev_i, rl_next_i, rl_final, rl_anchors, rl_anchors_r;
     
     //seq_dup = (char *)malloc(sizeof(char)*(len + 1));
     seq_strand = strdup(seq);
     error = MISMATCH;


     char *quality_clipping = (char *) malloc(sizeof(char) * 50);
     seq_dup = (char *) malloc(sizeof(char) * (len + 1));

     new_results_list(&rl_prev, bwt_optarg->filter_read_mappings);
     new_results_list(&rl_next, bwt_optarg->filter_read_mappings);
     new_results_list(&rl_prev_i, bwt_optarg->filter_read_mappings);
     new_results_list(&rl_next_i, bwt_optarg->filter_read_mappings);
     new_results_list(&rl_final, bwt_optarg->filter_read_mappings);

     array_list_t *tmp_mapping_list = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, 
						     COLLECTION_MODE_ASYNCHRONIZED);



     int num_errors = 1;
     int frag_size = 30;
     int aux_frag_size = len / (num_errors + 1);

     if (aux_frag_size < frag_size)
       aux_frag_size = frag_size;

     //printf("aux_frag_size = %i, len = %i\n", aux_frag_size, len);
     new_results_list(&rl_anchors,   (num_errors + 1));
     new_results_list(&rl_anchors_r, (num_errors + 1));

     rl_anchors.num_results = 0;
     rl_anchors.read_index = 0;	        
     rl_anchors_r.num_results = 0;
     rl_anchors_r.read_index = 0;	        

     for (int type = 1; type >= 0; type--) {

       rl_final.num_results = 0;
       rl_final.read_index = 0;	        
       //printf("Process Strand(%i)\n", !type);
       //    printf("*** bwt.c: calling BWSearch1 with type = %d...\n", type);
	  if (type == 1) {
	    BWSearchCPU(code_seq,
			len, 
			index->backward,
			index->forward, 
			&rl_prev,
			&rl_next,
			&rl_prev_i,
			&rl_next_i,
			&rl_final,
			&rl_anchors,
			aux_frag_size,
			1,
			index->bwt_config.nA
			);
	    // BWSearch1VectorHelper(code_seq, start, end, k1, l1, ki1, li1, 
	    //			     index->backward, index->forward, 
	    //			     &r_list, index->bwt_config.nA);
	  } else {

	    	    BWSearchCPU(code_seq,
			len, 
			index->forward_rev,
			index->backward_rev, 
			&rl_prev,
			&rl_next,
			&rl_prev_i,
			&rl_next_i,
			&rl_final,
			&rl_anchors_r,
			aux_frag_size,
			0,
			index->bwt_config.nA
			);
	    
	    //BWSearch1VectorHelper(code_seq, start, end, ki0, li0, k0, l0, 
	    //		 index->forward_rev, index->backward_rev, 
	    //			     &r_list, index->bwt_config.nA);      

	       if (rl_final.num_results) {
		    seq_reverse_complementary(seq_strand, len);
	       }
	  }

	  //printf("*** bwt.c: calling BWSearch1 with type = %d (num_results = %d). Done !!\n", type, r_list.num_results);

	  for (size_t ii = 0; ii < rl_final.num_results; ii++) {
	       r = &rl_final.list[ii];

	       //printf("Errors number %d\n", r->num_mismatches);
	       //if(r == NULL){printf("ERROR: Result list position null");exit(-1);}
	       if(!r->num_mismatches)
		    error = 0;
	       else
		    error = r->err_kind[0];

	       pos = r->err_pos[0];
	       //pos = r->position[0];
	  

	       if (type) {
		    direction = r->dir;
	       } else {
		    direction = !r->dir;
	       }

	       len_calc = len;
	       if (error == DELETION) {
		    len_calc--;
	       } else if (error == INSERTION) {
		    len_calc++;
	       }

	       // generating cigar
	       sprintf(quality_clipping, "%i", NONE_HARD_CLIPPING);
	       if (error == 0) {
		 sprintf(cigar, "%lu=%c", len, '\0');
		 num_cigar_ops = 1;
		 memcpy(seq_dup, seq_strand, len);
		 seq_dup[len] = '\0';
		 
	       } else if (error == MISMATCH) {
		 if (pos == 0) {
		   //Positive strand
		   if(type) { 
		     sprintf(cigar, "1S%luM%c", len-1, '\0'); 
		     start_mapping++;
		   }
		   else { 
		     sprintf(cigar, "%luM1S%c", len-1, '\0'); 
		   }
		   num_cigar_ops = 2;
		 } else if (pos == len - 1) {
		   //Positive strand
		   if(type) { 
		     sprintf(cigar, "%luM1S%c", len - 1, '\0'); 
		   }
		   else{ 
		     sprintf(cigar, "1S%luM%c", len-1, '\0'); 
		     start_mapping++;
		   }
		   num_cigar_ops = 2;
		 } else {
		   sprintf(cigar, "%luM%c", len, '\0');
		   num_cigar_ops = 1;
		 }
		 memcpy(seq_dup, seq_strand, len);
		 seq_dup[len] = '\0';
		 //printf("MISMATCH\n");
		 
	       } else if (error == INSERTION) {
		 //printf("INSERTION\n");
		 if (pos == 0) {
		   if(type) {
		     sprintf(cigar, "1M1D%luM%c", len - 1, '\0'); 
		   }
		   else{ 
		     sprintf(cigar, "%luM1D1M%c", len - 1, '\0'); 
		   }	      
		 } else if (pos == len - 1) {
		   if(type) { 
		     sprintf(cigar, "%luM1D1M%c", len - 1, '\0'); 
		   }
		   else{ 
		     sprintf(cigar, "1M1D%luM%c", len - 1, '\0'); 
		   }
		 } else {
		   
		   if(type) {
		     if(r->dir)
		       sprintf(cigar, "%iM1D%luM%c", pos, len - pos, '\0');
		     else
		       sprintf(cigar, "%iM1D%luM%c", pos + 1, len - pos - 1, '\0');
		   } else { 
		     if(r->dir)
		       sprintf(cigar, "%luM1D%dM%c", len - pos, pos, '\0');
		     else
		       sprintf(cigar, "%luM1D%dM%c", len - pos - 1, pos + 1, '\0');
		   }
		 }
		 num_cigar_ops = 3;
		 memcpy(seq_dup, seq_strand, len);
		 seq_dup[len] = '\0';
	       } else if (error == DELETION) {	     
		 //printf("DELETION\n");
		 if (pos == 0) {
		   if(type) { sprintf(cigar, "1I%luM%c", len -1, '\0'); }
		   else{ 
		     sprintf(cigar, "%luM1I%c", len -1, '\0'); 
		     //		   start_mapping++;
		   }
		   
		   num_cigar_ops = 2;		
		 } else if (pos == len - 1) {
		   if(type) { 
		     sprintf(cigar, "%luM1I%c", len -1, '\0'); 
		     //		   start_mapping++;
		   }
		   else{ 
		     sprintf(cigar, "1I%luM%c", len -1, '\0'); 
		   }
		   num_cigar_ops = 2;
		 } else {
		   if(type) { 
		     sprintf(cigar, "%dM1I%luM%c", pos, len - pos - 1, '\0'); 
		   }
		   else{ 
		     sprintf(cigar, "%luM1I%dM%c", len - pos - 1, pos, '\0'); 
		   }
		   num_cigar_ops = 3;
		 }
		 memcpy(seq_dup, seq_strand , len );
		 seq_dup[len] = '\0';
		 
	       }else{
		 printf("NUM MAPPINGS %lu -> POS %d -> ERROR %d -> (%lu):%s", num_mappings, pos, error, len, seq);
		 continue;
		 //exit(-1);
		 //error_debug = 1;
	       }
	       k_start = r->k;
	       l_start = r->l;
	       //}

	       //printf("%lu-%lu :: num_results %lu\n", k_start, l_start, r_list.num_results);
	       tot_alignments += (l_start - k_start);
	       if (tot_alignments >  bwt_optarg->filter_read_mappings) {
		 filter_exceeded = 1;
		 LOG_DEBUG_F("Filter exceeded: num. read mappings: %i (total: %i > filter %i)\n", 
			     l_start - k_start, tot_alignments, bwt_optarg->filter_read_mappings);
		 break;
	       }

	       for (size_t j = k_start; j <= l_start; j++) {
		 key = (direction)
		   ? size_SA(index->forward) - get_SA(j, index->forward) - len_calc - 1
		   : get_SA(j, index->backward);

		 idx = binsearch(index->karyotype->offset, index->karyotype->size, key);
		 if(key + len_calc <= index->karyotype->offset[idx]) {
		   start_mapping = index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]) + 1;
		   // save all into one alignment structure and insert to the list
		   //printf("*****Alignments %i:%lu\n", idx, start_mapping);
		   alignment = alignment_new();
		   alignment_init_single_end(NULL, strdup(seq_dup), strdup(quality_clipping), !type, 
					     idx - 1, //index->karyotype->chromosome + (idx-1) * IDMAX,
					     start_mapping, //index->karyotype->start[idx-1] + (key - index->karyotype->offset[idx-1]), 
					     strdup(cigar), num_cigar_ops, 254, 1, (num_mappings > 0), 0, NULL, alignment);

		   array_list_insert((void*) alignment, tmp_mapping_list);
		 }
	       }//end for k and l
	  }//end for 
	  
	  if (filter_exceeded) {
	    array_list_clear(tmp_mapping_list, (void *)alignment_free);
	    array_list_set_flag(2, mapping_list);
	    break;
	  }

	  //
	  //
	  //Search for equal BWT mappings and set the mappings that will be delete
	  int n_mappings = array_list_size(tmp_mapping_list);
	  alignment_t *alig_1, *alig_2;
	  unsigned int *delete_mark = (unsigned int *)calloc(n_mappings, sizeof(unsigned int));
	  const int max_distance = 10;
	  //printf("------------------Num mappings %i---------------\n", n_mappings);
     
	  for (int a1 = n_mappings - 1; a1 >= 1; a1--) {
	    if (!delete_mark[a1]) {
	      alig_1 = array_list_get(a1, tmp_mapping_list);
	      //printf("Alig1(%i): %i chromosome, %i seq_strand, %s cigar:\n", a1, alig_1->chromosome, alig_1->seq_strand, alig_1->cigar);
	      for (int a2 = a1 - 1; a2 >= 0; a2--) {
		alig_2 = array_list_get(a2, tmp_mapping_list);
		size_t dist = abs(alig_1->position - alig_2->position);
		//printf("\t Alig2(%i): %i chromosome, %i seq_strand, %i dist, %i delete mark, %s cigar\n", a2, alig_2->chromosome, alig_2->seq_strand, dist, delete_mark[a2],alig_2->cigar );
		if (alig_1->chromosome == alig_2->chromosome &&
		    dist < max_distance && 
		    !delete_mark[a2]) {
		  //Same chromosome && same position
		  if (alig_1->num_cigar_operations < alig_2->num_cigar_operations) {
		    //printf("\tSet read %i\n", a2);
		    delete_mark[a2] = 1;
		  } else {
		    //printf("\tSet read %i\n", a2);
		    delete_mark[a1] = 1;
		  }
		}
	      }
	    }
	  }
	  //Delete all set mappings
	  int primary_delete = 0;
	  for (int m = n_mappings - 1; m >= 0; m--) {
	    alig_1 = array_list_remove_at(m, tmp_mapping_list);
	    if (delete_mark[m]) {
	      if (!is_secondary_alignment(alig_1)) { primary_delete = 1; }
	      alignment_free(alig_1);
	    } else {
	      array_list_insert((void*) alig_1, mapping_list);
	      set_secondary_alignment(num_mappings > 0, alig_1);
	      num_mappings++;
	    }
	  }
	  if (primary_delete) { 
	    alig_1 = array_list_get(0, mapping_list); 
	    set_secondary_alignment(0, alig_1);
	  }
	  free(delete_mark);	  	  
     } // end for type 

     if (array_list_size(mapping_list) == 0) {
       if (array_list_get_flag(mapping_list) == 1) {

	 //====================================================================================
	 //printf("BWT(+): FORWARD(k-l)%i:%lu-%lu BACKWARD(k-l)%i:%lu-%lu\n", forw1_nt, ki1[forw1_nt], li1[forw1_nt], back1_nt, k1[back1_nt], l1[back1_nt]);
	 //printf("BWT(-): FORWARD(k-l)%i:%lu-%lu BACKWARD(k-l)%i:%lu-%lu\n", forw0_nt, ki0[forw0_nt], li0[forw0_nt], back0_nt, k0[back0_nt], l0[back0_nt]);

	 size_t k_start, l_start, k_start_r, l_start_r;

	 //printf("Anchors Results (%i): \n", rl_anchors.num_results);
	 for (uintmax_t r1 = 0; r1 <  rl_anchors.num_results; r1++) {

	   r = &rl_anchors.list[r1];
	   //printf("Anchor %i: start = %i, end = %i, pos = %i, dir = %i, len = %i\n", r1, (int)r->start, (int)r->end, (int)r->pos, (int)r->dir, (int)(r->end - r->start + 1));

	   if (r->end == len-1) {
	     __bwt_generate_anchor_list(r->k, 
					r->l, r->end - r->start + 1, bwt_optarg, 
					index, 1/*+*/, mapping_list, BACKWARD_ANCHOR, 0);	     
	   } else if (r->start == 0) {
	     __bwt_generate_anchor_list(r->k,
					r->l, r->end - r->start + 1, bwt_optarg, 
					index, 1/*+*/,  mapping_list, FORWARD_ANCHOR, 0);
	   }
	 }

	 for (uintmax_t r1 = 0; r1 <  rl_anchors_r.num_results; r1++) {
	   r = &rl_anchors_r.list[r1];
	   //printf("Anchor %i: start = %i, end = %i, pos = %i, dir = %i, len = %i\n", r1, (int)r->start, (int)r->end, (int)r->pos, (int)r->dir, (int)(r->end - r->start + 1));

	   if (r->end == len-1) {
	     __bwt_generate_anchor_list(r->k, 
					r->l, r->end - r->start + 1, bwt_optarg, 
					index, 0/*-*/,  mapping_list, FORWARD_ANCHOR, 0);

	   } else if (r->start == 0) {
	     __bwt_generate_anchor_list(r->k,
					r->l, r->end - r->start + 1, bwt_optarg, 
					index, 0/*-*/,  mapping_list, BACKWARD_ANCHOR, 0);
	   }
	 }
       }

     } else {
       int header_len;
       size_t num_mapping = array_list_size(mapping_list);
       //printf("BWT:Report Total mappings %i\n", num_mappings);
       for (int i = 0; i < num_mapping; i++) {
	 alignment = array_list_get(i, mapping_list);
	 header_len = strlen(read->id);
	 alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	 get_to_first_blank(read->id, header_len, alignment->query_name);
	 bwt_cigar_cpy(alignment, read->quality);
	 //alignment->quality = strdup(&(batch->quality[batch->data_indices[i]]));                                                                                     
	 //************************* OPTIONAL FIELDS ***************************//
	 alignment = add_optional_fields(alignment, num_mappings, len);
       }
       array_list_set_flag(0, mapping_list);
       //printf("########## EXACT! #########\n");
     }

     free(rl_prev.list);
     free(rl_next.list);
     free(rl_next_i.list);
     free(rl_prev_i.list);
     free(rl_final.list);
     free(rl_anchors.list);
     free(rl_anchors_r.list);

     free(code_seq);
     free(seq_strand);
     //free(k0);
     //free(l0);
     //free(k1);
     //free(l1);
     //free(ki0);
     //free(li0);
     //free(ki1);
     //free(li1);

     free(seq_dup);
     free(quality_clipping);
     array_list_free(tmp_mapping_list, NULL);

     return array_list_size(mapping_list);

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
/*
size_t bwt_map_inexact_read(fastq_read_t *read, 
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    array_list_t *mapping_list) {
  
  size_t padding = array_list_size(mapping_list);
  
  size_t num_mappings = bwt_map_inexact_seq(read->sequence, 
					    bwt_optarg, index, 
					    mapping_list);

  alignment_t *mapping;
  unsigned int header_len, quality_len;
  unsigned int quality_type;
  
  for (int i = 0; i < num_mappings; i++) {
    mapping = (alignment_t *) array_list_get(i + padding, mapping_list);
    
    if(mapping != NULL){
      header_len = strlen(read->id) + 1;
      mapping->query_name = (char *)malloc(sizeof(char)*header_len);
      get_to_first_blank(read->id, header_len, mapping->query_name);
  
      quality_len = strlen(read->quality) + 1;
      quality_type = atoi(mapping->quality);
      //printf("Quality len from batch: %i\n", quality_len);
      free(mapping->quality);
      mapping->quality = (char *)malloc(sizeof(char)*(quality_len + 1));
  
      //printf("Read:%s\n", mapping->sequence);
      if (mapping->seq_strand == 0) {
	   memcpy(mapping->quality, 
		  read->quality, 
		  quality_len);
      } else {
	   reverse_str(read->quality,
		       mapping->quality, quality_len);
      }
	//mapping->quality[quality_len] = '\0';
	//printf("(%i)NORMAL : %s\n", mapping->seq_strand, mapping->quality);
      //}      
      //mapping->quality = read->quality;
    }else{
      printf("Error to extract item\n");
    }
  }
  
  return num_mappings;
}
*/
//-----------------------------------------------------------------------------
/*
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
*/

//-----------------------------------------------------------------------------
/*
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
  size_t quality_len;
  //struct timeval start_time, end_time;
  //double timer, parallel_t, sequential_t;
  unsigned int j, header_len;
  size_t read_id = 0;

  for (unsigned int i = 0; i < num_reads; i++) {
    individual_mapping_list_p[i] = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
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
	  //printf("Header len %i\n", header_len);
	  mapping->query_name = (char *)malloc(sizeof(char)*header_len + 1);
	  get_to_first_blank(&(batch->header[batch->header_indices[i]]), header_len, mapping->query_name);
	  //printf("--->header id %s to %s\n",read->id, header_id);	  
	  bwt_cigar_cpy_batch(mapping, i, batch);
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
      //memcpy(&(unmapped_batch->seq[unmapped_batch->data_indices[read_pos - 1]]),      &(batch->seq[batch->data_indices[i]]),      quality_len);
      strcpy_capitalize(&(unmapped_batch->seq[unmapped_batch->data_indices[read_pos - 1]]),      &(batch->seq[batch->data_indices[i]]),      quality_len);
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
*/
//-----------------------------------------------------------------------------				  
/*
size_t bwt_map_inexact_array_list(array_list_t *reads,
				  bwt_optarg_t *bwt_optarg, 
				  bwt_index_t *index,
				  array_list_t **lists,
				  size_t *num_unmapped, 
				  size_t *unmapped_indices) {

  size_t header_len, total_mappings;
  size_t num_threads = bwt_optarg->num_threads;
  size_t num_reads = array_list_size(reads);
  size_t chunk = MAX(1, num_reads/(num_threads*10));
  fastq_read_t* read;
  fastq_read_pe_t* read_pe;
  /*
  *num_mapped = 0;
  *num_unmapped = 0;
  //printf("%i reads\n", num_reads);
  if (single_end) {
    #pragma omp parallel for private(read) schedule(dynamic, chunk)
    for (size_t i = 0; i < num_reads; i++) {
      read_p = (fastq_read_t *)array_list_get(i, reads);
      //printf("Extract...\n");
      //printf("%s\n", read_p->sequence);
      bwt_map_inexact_seq(read_p->sequence, 
			  bwt_optarg, index, 
			  lists[i]);
      
    }
  } else {
    int j = 0;
    #pragma omp parallel for private(read_pe) schedule(dynamic, chunk)
    for (size_t i = 0; i < num_reads; i++) {
      read_pe = (fastq_read_pe_t *)array_list_get(i, reads);
      //printf("Extract...\n");
      //printf("%s\n", read_p->sequence);
      bwt_map_inexact_seq(read_pe->sequence1, 
			  bwt_optarg, index, 
			  lists[j++]);
  
      bwt_map_inexact_seq(read_pe->sequence2, 
			  bwt_optarg, index, 
			  lists[j++]);
      
    }

    if (single_end) {
      size_t num_mappings = 0;
      for (size_t i = 0; i < num_reads; i++) {
	num_mappings = array_list_size(lists[i]);
	if (num_mappings > 0) {
	  read_p = array_list_get(i, reads);
	  total_mappings += num_mappings;
	  array_list_set_flag(1, lists[i]);
	  for (size_t j = 0; j < num_mappings; j++) {
	    alignment = (alignment_t *) array_list_get(j, lists[i]);
	    header_len = strlen(read->id);
	    alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	    get_to_first_blank(read->id, header_len, alignment->query_name);
	    
	    bwt_cigar_cpy(alignment, read->quality);
	    //	  free(alignment->quality);
	    //alignment->quality = strdup(&(batch->quality[batch->data_indices[i]]));
	  }
	} else {
	  unmapped_indices[(*num_unmapped)++] = i;
	  array_list_set_flag(0, lists[i]);
	  //	printf("\tbwt.c: bwt_map_inexact_batch_by_filter, setting flag to 0 for list %i\n", i);
	}
      }
    } else {
      int num_mappings[2];
      char *id, *quality;
      int type, read_index, num_lists = array_list_size(lists);
      for (size_t i = 0; i < num_lists; i+=2) {
      //for (size_t i = 0; i < num_reads; i++) {
	//num_mappings[0] = array_list_size(lists[i]);
	//num_mappings[1] = array_list_size(lists[i+1]);
	
	read_index = i / 2;
	read_pe = array_list_get(read_index, reads);

	type = 0;
	for (int j = 0; j < 2; j++) {
	  num_mappings[j] = array_list_size(lists[i+j]);
	  array_list_set_flag(0, lists[i+j]);
	  	  
	  if (num_mappings[j] > 0) {
	    if (j == 0) {
	      id = read_pe->id1;
	      quality = read_pe->quality1;
	    } else {
	      id = read_pe->id2;
	      quality = read_pe->quality2;
	    }
	    total_mappings += num_mappings[j];
	    array_list_set_flag(1, lists[i+j]);	    
	    for (size_t k = 0; k < num_mappings[j]; k++) {
	      alignment = (alignment_t *) array_list_get(k, lists[i+j]);	    
	      header_len = strlen(id);
	      alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	      //printf("Process %s\n", read_p->id);
	      //TODO: Extract the next function out 
	      get_to_first_blank(id, header_len, alignment->query_name);
	      bwt_cigar_cpy(alignment, quality);
	    }
	    } else {*/
	    /**********************************/
	    /* TYPE 1: PAIR 1, NOT MAPPED     */
	    /* TYPE 2: PAIR 2, NOT MAPPED     */
	    /* TYPE 3: PAIR 1 & 2, NOT MAPPED */
	    /**********************************/
	    /*if (j == 0) {
	      id = read_pe->id1;
	      quality = read_pe->quality1;
	      type += 1;
	    } else {
	      id = read_pe->id2;
	      quality = read_pe->quality2;
	      type += 2;
	    }
	  }
	}
	if (type) {
	  unmapped_indices[(*num_unmapped)] = read_index;
	  unmapped_pairs[(*num_unmapped)++] = type;
	}
      }
    }
  }
  return total_mappings;
}
*/
//-----------------------------------------------------------------------------
alignment_t* add_optional_fields(alignment_t *alignment, size_t n_mappings, int read_length) {
  char *p, *optional_fields;
  size_t optional_fields_length = 100;
  int distance;
  int AS;
  int cigar_len;

  distance = alignment->map_quality;    
  AS = (read_length * 1) - (distance * 0.4);
  AS = AS * 100 / read_length;

  optional_fields = (char *)calloc(optional_fields_length, sizeof(char));
  p = optional_fields;

  sprintf(p, "ASi");
  p += 3;
  memcpy(p, &AS, sizeof(int));
  p += sizeof(int);
  
  sprintf(p, "NHi");
  p += 3;
  memcpy(p, &n_mappings, sizeof(int));
  p += sizeof(int);
  
  sprintf(p, "NMi");
  p += 3;
  cigar_len = strlen(alignment->cigar);
  
  alignment->map_quality = 255;
  
  memcpy(p, &distance, sizeof(int));
  p += sizeof(int);
  alignment->optional_fields_length = p - optional_fields;
  alignment->optional_fields = optional_fields;
  
  return alignment;

}

//-----------------------------------------------------------------------------
/*
size_t bwt_map_inexact_array_list(array_list_t *reads,
				  bwt_optarg_t *bwt_optarg, 
				  bwt_index_t *index,
				  array_list_t **lists,
				  size_t *num_unmapped, 
				  size_t *unmapped_indices) {

  alignment_t *alignment;
  size_t header_len, num_mappings, total_mappings;
  size_t num_threads = bwt_optarg->num_threads;
  size_t num_reads = array_list_size(reads);
  size_t chunk = MAX(1, num_reads/(num_threads*10));
  fastq_read_t* fq_read;
 
  *num_unmapped = 0;
  //printf("%i reads\n", num_reads);
  #pragma omp parallel for private(fq_read) schedule(dynamic, chunk)
  for (size_t i = 0; i < num_reads; i++) {
    fq_read = (fastq_read_t *) array_list_get(i, reads);
    //printf("Extract...\n");
    //printf("%s\n", read_p->sequence);
    bwt_map_inexact_seq(fq_read->sequence, 
			bwt_optarg, index, 
			lists[i]);
  }

  for (size_t i = 0; i < num_reads; i++) {
    num_mappings = array_list_size(lists[i]);
    total_mappings += num_mappings;
    fq_read = (fastq_read_t *) array_list_get(i, reads);

    if (num_mappings > 0) {
      array_list_set_flag(1, lists[i]);
      for (size_t j = 0; j < num_mappings; j++) {
	alignment = (alignment_t *) array_list_get(j, lists[i]);
	header_len = strlen(fq_read->id);
	alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	//printf("Process %s\n", fq_read->id);                                                                                                                         
	get_to_first_blank(fq_read->id, header_len, alignment->query_name);
	bwt_cigar_cpy(alignment, fq_read->quality);
	//alignment->quality = strdup(&(batch->quality[batch->data_indices[i]]));                                                                                     

	alignment = add_optional_fields(alignment, num_mappings);

      }
    } else {
      unmapped_indices[(*num_unmapped)++] = i;
      array_list_set_flag(0, lists[i]);
      //      printf("\tbwt.c: bwt_map_inexact_batch_by_filter, setting flag to 0 for list %i\n", i);                                                       
    }
  }
}
*/
//-----------------------------------------------------------------------------
/*
void bwt_map_inexact_array_list_by_filter(array_list_t *reads,
					  bwt_optarg_t *bwt_optarg, 
					  bwt_index_t *index,
					  array_list_t **lists,
					  size_t *num_unmapped, 
					  size_t *unmapped_indices) {
  alignment_t *alignment;
  size_t header_len, num_mappings, total_mappings;
  size_t num_threads = bwt_optarg->num_threads;
  size_t num_reads = array_list_size(reads);
  size_t chunk = MAX(1, num_reads/(num_threads*10));
  fastq_read_t* fq_read;
  *num_unmapped = 0;

  for (size_t i = 0; i < num_reads; i++) {
    fq_read = (fastq_read_t *) array_list_get(i, reads);
    array_list_set_flag(1, lists[i]);
    num_mappings = 0;

    //if (fq_read->length < 200) {
    printf("BWT %s : \n", fq_read->id);
    num_mappings = bwt_map_inexact_seq(fq_read->sequence, 
				       bwt_optarg, index, 
				       lists[i]);
    //}
    if (num_mappings > 0) {
      array_list_set_flag(1, lists[i]);
      for (size_t j = 0; j < num_mappings; j++) {
	alignment = (alignment_t *) array_list_get(j, lists[i]);
	header_len = strlen(fq_read->id);
	alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	
	get_to_first_blank(fq_read->id, header_len, alignment->query_name);
	bwt_cigar_cpy(alignment, fq_read->quality);
	

	alignment = add_optional_fields(alignment, num_mappings);

      }
    } else if (array_list_get_flag(lists[i]) != 2) {
	unmapped_indices[(*num_unmapped)++] = i;
	array_list_set_flag(0, lists[i]);
    }
  } 
}
*/
//-----------------------------------------------------------------------------
/*
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
      //      printf("bwt, read #%i of %i (len = %i) : %s\n", i, num_reads, strlen(&(batch->seq[batch->data_indices[i]])), &(batch->seq[batch->data_indices[i]]));

      num_mappings = bwt_map_inexact_seq(&(batch->seq[batch->data_indices[i]]), 
					 bwt_optarg, index, 
					 lists[i]);

      if (num_mappings > 0) {
	mapped[(*num_mapped)++] = i;
	array_list_set_flag(1, lists[i]);

	//	printf("\tbwt.c: bwt_map_inexact_batch_by_filter, setting flag to 1 for list %i (num_mappings  %d)\n", i, num_mappings);
	for (size_t j = 0; j < num_mappings; j++) {
	  alignment = (alignment_t *) array_list_get(j, lists[i]);

	  header_len = batch->header_indices[i+1] - batch->header_indices[i];
	  alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	  get_to_first_blank(&(batch->header[batch->header_indices[i]]), header_len, alignment->query_name);

	  bwt_cigar_cpy_batch(alignment, i, batch);
	  //	  free(alignment->quality);
	  //alignment->quality = strdup(&(batch->quality[batch->data_indices[i]]));
	}
      } else {
	unmapped[(*num_unmapped)++] = i;
	array_list_set_flag(0, lists[i]);
	//	printf("\tbwt.c: bwt_map_inexact_batch_by_filter, setting flag to 0 for list %i\n", i);
      }
    }
  } else {
    printf("not yet implemented, true filter\n");
    exit(-1);
  }
}
*/

//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seeds_seq(char *seq, size_t seed_size, size_t min_seed_size,
				 bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				 array_list_t *mapping_list) {
  
  size_t len = strlen(seq);
  size_t offset, num_seeds = len / seed_size;
  uint8_t *code_seq = (uint8_t *) malloc((len + 1) * sizeof(uint8_t));
  int seq_id = 0;
  //char aux_seq[50];


  encode_bases(code_seq, seq, len, index->bwt_config.table);

  // first 'pasada'
  offset = 0;
  for (size_t i = 0; i < num_seeds; i++) {
    //memcpy(aux_seq, seq + offset, seed_size);
    //aux_seq[offset + seed_size] = '\0';
    //printf("Seq(%i-%i)\n", offset, offset + seed_size - 1);
    
    bwt_map_inexact_seed(code_seq, len, offset, offset + seed_size - 1,
			 bwt_optarg, index, mapping_list, seq_id);

    //printf("\tregions found: %i\n", array_list_size(mapping_list));
    offset += seed_size;
    seq_id++;

  }


  // special processing for the last seed !!
  if (len % seed_size >= min_seed_size) {
    bwt_map_inexact_seed(code_seq, len, offset, len - 1,
			 bwt_optarg, index, mapping_list, seq_id);
  }

  free(code_seq);

  return array_list_size(mapping_list);
  
}

//-----------------------------------------------------------------------------
size_t bwt_map_inexact_seeds_by_region(int start_position, int end_position, 
				       int strand_target, 
				       int chromosome_target, int start_target,
				       int end_target,
				       char *seq, size_t seed_size,
				       size_t min_seed_size,
				       bwt_optarg_t *bwt_optarg,
				       bwt_index_t *index, 
				       array_list_t *mapping_list) {
  //array_list_t *aux_list = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  //array_list_t *regions_list = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  if (seed_size <= 10) { seed_size = 16; }
  int len_region = end_position - start_position;
  size_t len = strlen(seq);
  size_t offset, num_seeds = len_region / seed_size;
  uint8_t *code_seq = (uint8_t *) malloc((len + 1) * sizeof(uint8_t));
  int seq_id = 1;
  char aux_seq[50];

  encode_bases(code_seq, seq, len, index->bwt_config.table);

  int padding_regions = 0;
  int num_regions;
  int distance_merge = 0;
  int merge_distance = 5;
  int found, found_id;



  region_t *found_fusion = NULL;
  int region_search;
  if (start_position <= 10) {
    region_search = 0;
  } else {
    region_search = 1;
  }

  //seed_size = ;
  //printf(" REGION INTERVAL (%i)[%i:%lu-%lu]\n", strand_target, chromosome_target, start_target, end_target);
  found_fusion = 0;
  // first 'pasada'
  offset = start_position;
  for (size_t i = 0; i < num_seeds; i++) {
    memcpy(aux_seq, seq + offset, seed_size);
    aux_seq[seed_size] = '\0';


    bwt_map_inexact_seed_by_region(code_seq, len, offset, offset + seed_size - 1,
				   bwt_optarg, index, mapping_list, seq_id++, 
				   chromosome_target, start_target, end_target, 
				   strand_target);
    //num_regions  = array_list_size(aux_list);

    //printf("Seq(%i-%i) Tot Regions %i: %s\n", offset, offset + seed_size - 1, num_regions - padding_regions, aux_seq);

    //found = 0;
    /*    
    for (int j1 = 0; j1 < num_regions; j1++) {
      region_t *region = array_list_get(j1, aux_list);
      printf("\t%i:%lu-%lu\n", region->chromosome_id, region->start, region->end);
      for (int j2 = 0; j2 < array_list_size(regions_list); j2++) {
	region_t *region_prev = array_list_get(j2, regions_list);
	if (region_prev->end + merge_distance >= region->start) {
	  region_prev->end = region->end;
	  region_prev->seq_end = region->seq_end;
	  found = 1;
	  found_id = j1;	  
	}
      }      
    }
    
    if (!found) {
      for (int j1 = array_list_size(aux_list) - 1; j1 >= 0; j1--) {
	region_t *region = array_list_remove_at(j1, aux_list);
	array_list_insert(region, regions_list);
      }
    } else {
      for (int j1 = 0; j1 < array_list_size(aux_list); j1++) {
	region_t *region = array_list_get(j1, aux_list);
	if (j1 != found_id) {
	  array_list_insert(region, regions_list);
	} else {
	  region_bwt_free(region);
	}
      }
      array_list_clear(aux_list, NULL);
    }
    */
    //merge_distance += seed_size;
	   
    offset += seed_size;
  }

  
  if (len_region % seed_size >= seed_size / 4) {
    bwt_map_inexact_seed_by_region(code_seq, len, end_position - seed_size, end_position - 1,
				   bwt_optarg, index, mapping_list, seq_id++, 
				   chromosome_target, start_target, end_target, 
				   strand_target);
  }
  

  /* found = 0; */
  /* int max_region = 0; */
  /* printf("totallllllllll %i\n", array_list_size(regions_list) ); */

  /* for (int i = 0; i < array_list_size(regions_list); i++) { */
  /*   region_t *region = array_list_get(i, regions_list);       */

  /*   if (region->seq_end - region->seq_start > seed_size) { */
  /*     found = 1; */
  /*     if (region->seq_end - region->seq_start > max_region) { */
  /* 	max_region = region->seq_end - region->seq_start; */
  /* 	found_fusion = region; */
  /*     } */
  /*   } */
  /* } */

  /* if (found) { */
  /*   array_list_insert(found_fusion, mapping_list); */
  /*   for (int i = 0; i < array_list_size(regions_list); i++) { */
  /*     region_t *region = array_list_get(i, regions_list);     */
  /*     if (region != found_fusion) { region_bwt_free(region); } */
  /*   } */
  /* } else { */
  /*   found_fusion = array_list_get(0, regions_list); */
  /*   for (int i = 1; i < array_list_size(regions_list); i++) { */
  /*     region_t *region = array_list_get(i, regions_list);     */
  /*     if (region_search == 0) { */
  /* 	if (found_fusion->start < region->start) { */
  /* 	  found_fusion = region; */
  /* 	} */
  /*     } else { */
  /* 	if (found_fusion->start > region->start) { */
  /* 	  found_fusion = region; */
  /* 	} */
  /*     }     */

  /*     array_list_insert(found_fusion, mapping_list); */
  /*     for (int i = 0; i < array_list_size(regions_list); i++) { */
  /* 	region_t *region = array_list_get(i, regions_list);     */
  /* 	if (region != found_fusion) { region_bwt_free(region); } */
  /*     } */
  /*   } */
  /* } */

  //array_list_free(regions_list, NULL);
  //array_list_free(aux_list, NULL);

  free(code_seq);

  return array_list_size(mapping_list);
  
}

//-----------------------------------------------------------------------------
/*
size_t bwt_map_exact_seeds_seq2(int padding_left, int padding_right, 
			       char *seq, size_t seed_size, size_t min_seed_size,
			       bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
			       array_list_t *mapping_list, unsigned char step_id) {
  size_t len = strlen(seq);  
  //size_t extra_seed_size = 15;
  int seed_id = 0;

  size_t offset, num_seeds;
  
  //printf(" len=%i, seed_size=%i, num_seeds=%i\n", len, seed_size, num_seeds);
  char *code_seq = (char *) calloc(len, sizeof(char));

  replaceBases(seq, code_seq, len);

  //Second first seed
  padding_left = seed_size / 2;
  //char aux_seq[50];

  //memcpy(aux_seq, seq + padding_left, seed_size);
  //aux_seq[seed_size] = '\0';
  //printf("Seq(%i-%i): %s\n", padding_left, padding_left + seed_size, aux_seq);

  int extra_seed_size = 16;//seed_size;
  int first_seed = 16;

  bwt_map_exact_seed(code_seq, len, 0, first_seed - 1,
  		     bwt_optarg, index, mapping_list, seed_id++);

  bwt_map_exact_seed(code_seq, len, padding_left, padding_left + extra_seed_size - 1,
  		     bwt_optarg, index, mapping_list, seed_id++);

  num_seeds = (len - first_seed*2) / seed_size;

  // first 'pasada'
  offset = extra_seed_size;
  for (size_t i = 0; i < num_seeds; i++) {
    //memcpy(aux_seq, seq + offset, seed_size);
    //aux_seq[seed_size] = '\0';
    //printf("Seq(%i-%i): %s\n", offset, offset + seed_size - 1, aux_seq);
    bwt_map_exact_seed(code_seq, len, offset, offset + seed_size - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    offset += seed_size;
  }
  
  //Last seed
  if ((len - first_seed*2) % seed_size >= 0) {
    bwt_map_exact_seed(code_seq, len, len - first_seed - seed_size, len - first_seed - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    offset += seed_size;
    //} else {
    //padding_right = seed_size / 2;
  }

  bwt_map_exact_seed(code_seq, len, len - first_seed, len - 1,
		     bwt_optarg, index, mapping_list, seed_id++);
  //memcpy(aux_seq, seq + (len - padding_right - seed_size), seed_size);
  //aux_seq[seed_size] = '\0';
  //printf("Seq+(%i-%i): %s\n", len - padding_right - seed_size, len - 1 - padding_right, aux_seq);
  padding_right = padding_left;
  bwt_map_exact_seed(code_seq, len, len - padding_right - extra_seed_size, len - padding_right - 1,
  		     bwt_optarg, index, mapping_list, seed_id);

  free(code_seq);

  return array_list_size(mapping_list);
  
}
*/

/*
size_t bwt_map_seeds_IA(int padding_left,
			int padding_right,
			char *seq, 
			size_t seed_size, size_t min_seed_size,
			bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
			array_list_t *mapping_list, unsigned char step_id) {
  size_t len = strlen(seq);
  int seed_id = 0;
  size_t offset, num_seeds;
  //printf(" len=%i, seed_size=%i, num_seeds=%i\n", len, seed_size, num_seeds);
  uint8_t *code_seq = (uint8_t *) malloc(len * sizeof(uint8_t));
  int num_regions_start, num_regions_end;


  encode_bases(code_seq, seq, len, index->bwt_config.table);

  //First Exact seed Start read and End
  num_regions_start = bwt_map_exact_seed(code_seq, len, 0, 20 - 1,
					 bwt_optarg, index, mapping_list, seed_id++);
  
  printf("START REGIONS:\n");
  for (int i = 0; i < num_regions_start; i++) {
    region_t *region = array_list_get(i, mapping_list);
    printf("\t Region (%i)[%i:%lu-%lu]\n", region->strand, region->chromosome_id, region->start, region->end);
  }

  //First Exact seed Start read and End
  num_regions_end = bwt_map_exact_seed(code_seq, len, len - 20, len - 1,
				       bwt_optarg, index, mapping_list, seed_id++);
  printf("START REGIONS:\n");
  for (int i = num_regions_start; i < num_regions_end; i++) {
    region_t *region = array_list_get(i, mapping_list);
    printf("\t Region (%i)[%i:%lu-%lu]\n", region->strand, region->chromosome_id, region->start, region->end);
  }

  free(code_seq);
  return array_list_size(mapping_list);
}

*/


size_t bwt_map_exact_seeds_seq(int padding_left, int padding_right, 
			       char *seq, size_t seed_size, size_t min_seed_size,
			       bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
			       array_list_t *mapping_list, unsigned char step_id) {
  size_t len = strlen(seq);
  int seed_id = 0;

  size_t offset, num_seeds;
  if (seed_size <= 10) { seed_size = 16; }

  //printf(" len=%i, seed_size=%i, num_seeds=%i\n", len, seed_size, num_seeds);
  uint8_t *code_seq = (uint8_t *) malloc(len * sizeof(uint8_t));

  encode_bases(code_seq, seq, len, index->bwt_config.table);  

  //Second first seed
  padding_left = seed_size / 2;
  padding_right = padding_left;

  int extra_seed_size = 16;//seed_size;
  //int first_seed = 16;

  //Extra seed for splice junctions
  bwt_map_exact_seed(code_seq, len, padding_left, padding_left + extra_seed_size - 1,
  		     bwt_optarg, index, mapping_list, seed_id++);

  num_seeds = len / seed_size;

  // first 'pasada'
  offset = 0;
  for (size_t i = 0; i < num_seeds; i++) {
    bwt_map_exact_seed(code_seq, len, offset, offset + seed_size - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    offset += seed_size;
  }
  
  //Last seed
  if (len % seed_size > 0) {
    bwt_map_exact_seed(code_seq, len, len - seed_size, len - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
  }

  if (len % seed_size != padding_right) {
    //Extra seed for splice junctions
    bwt_map_exact_seed(code_seq, len, len - extra_seed_size - padding_right, len - padding_right - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
  }

  free(code_seq);

  return array_list_size(mapping_list);
  
}

void insert_seeds_and_merge(array_list_t *mapping_list, linked_list_t ***cals_list,  size_t max_cal_distance) {  
  for (int m = array_list_size(mapping_list) - 1; m >= 0; m--) {
    region_t *region = array_list_remove_at(m, mapping_list);
    int  chromosome_id = region->chromosome_id;
    int strand = region->strand;    
    //printf("Insert Region (%i)[%i:%lu|%i-%i|%lu]\n", region->strand, region->chromosome_id, 
    //	   region->start, region->seq_start, region->seq_end, region->end);
    my_cp_list_append_linked_list(cals_list[strand][chromosome_id], region, max_cal_distance, 1000);

    region_bwt_free(region);
  }
}    

/*
size_t bwt_generate_cals(char *seq, size_t seed_size, 
			 bwt_optarg_t *bwt_optarg,
			 cal_optarg_t *cal_optarg,
			 bwt_index_t *index, 
			 array_list_t *cal_list, 
			 unsigned int nchromosomes) {

  //printf("New function generate CALs\n");
  //NEW! Function! Fusion between Seeding and Caling
  array_list_t *mapping_list = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  const int nstrands = 2;
  nchromosomes = 30; //TODO: Parameter
  
  linked_list_t ***cals_list = (linked_list_t ***)malloc(sizeof(linked_list_t **)*nstrands);

  for (unsigned int i = 0; i < nstrands; i++) {
    cals_list[i] = (linked_list_t **)malloc(sizeof(linked_list_t *)*nchromosomes);
    for (unsigned int j = 0; j < nchromosomes; j++) {
      cals_list[i][j] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    }
  }

  size_t len = strlen(seq);
  int seed_id = 0;

  size_t offset, num_seeds;
  
  //printf(" len=%i, seed_size=%i, num_seeds=%i, seq=%s\n", len, seed_size, num_seeds, seq);
  char *code_seq = (char *) calloc(len, sizeof(char));

  encode_bases(code_seq, seq, len, index->bwt_config.table);
  //Second first seed
  if (seed_size <= 10) { seed_size = 16; }

  int padding_left = seed_size / 2;
  int padding_right = padding_left;

  int extra_seed_size = 16;//seed_size;
  //int first_seed = 16;

  //Extra seed for splice junctions
  bwt_map_exact_seed(code_seq, len, padding_left, padding_left + extra_seed_size - 1,
  		     bwt_optarg, index, mapping_list, seed_id++);

  insert_seeds_and_merge(mapping_list, cals_list,  len);
  num_seeds = len / seed_size;

  //Extra seed for splice junctions
  // first 'pasada'
  offset = 0;
  for (size_t i = 0; i < num_seeds; i++) {
    bwt_map_exact_seed(code_seq, len, offset, offset + seed_size - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    insert_seeds_and_merge(mapping_list, cals_list,  len);
    offset += seed_size;
  }
  
  //Last seed
  if (len % seed_size > 0) {
    bwt_map_exact_seed(code_seq, len, len - seed_size, len - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    insert_seeds_and_merge(mapping_list, cals_list,  len);
  }

  if (len % seed_size != padding_right) {
    //Extra seed for splice junctions
    bwt_map_exact_seed(code_seq, len, len - extra_seed_size - padding_right, len - padding_right - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    insert_seeds_and_merge(mapping_list, cals_list,  len);
  }

  free(code_seq);

  //Store CALs in Array List for return results                                                                                                   
  size_t start_cal, end_cal;
  size_t seq_start, seq_end;
  seed_region_t *s, *s_first, *s_last, *s_aux, *seed_region;
  cal_t *cal;
  short_cal_t *short_cal, *short_cal_aux;
  linked_list_iterator_t itr, itr2;
  linked_list_item_t *list_item_cal, *list_item_aux;
  size_t min_cal_size = 20;//cal_optarg->min_cal_size; TODO:PARAMETER
  const int max_intron_size = 500000; //TODO: Parameter
  size_t read_length = len;

  for (unsigned int j = 0; j < nchromosomes; j++) {
    for (unsigned int i = 0; i < nstrands; i++) {
      linked_list_iterator_init(cals_list[i][j], &itr);
      list_item_cal = linked_list_iterator_list_item_curr(&itr);
      while ((list_item_cal != NULL )) {
	short_cal = (short_cal_t *)list_item_cal->item;
	if (short_cal->end - short_cal->start + 1 >= min_cal_size) {
	  linked_list_t *list_aux = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_list)) {
	    //TODO: Change all parameters to seed_region_t
	    append_seed_region_linked_list(list_aux,
					   s->read_start, s->read_end, 
					   s->genome_start, s->genome_end, 
					   s->id);	    
	    seed_region_free(s);	    
	  }
	  cal = cal_new(j, i, short_cal->start, 
			short_cal->end, short_cal->num_seeds, 
			list_aux, short_cal->sr_duplicate_list);
 	  array_list_insert(cal, cal_list);
	  
	  short_cal->sr_duplicate_list = NULL;
	  
	  //============= Search if one seed map near to this CAL ================//
	  s_first = (seed_region_t *)linked_list_get_first(list_aux);	  
	  s_last = (seed_region_t *)linked_list_get_last(list_aux);
	  
	  if (s_first && s_first->read_start > seed_size) {
	    //Search start seed <-----
	    list_item_aux = list_item_cal->prev;
	    while (list_item_aux) {
	      short_cal_aux = list_item_aux->item;
	      if (!short_cal_aux ||
		  short_cal_aux->end - short_cal_aux->start + 1 >= min_cal_size || 
		  short_cal->start - short_cal_aux->end >= max_intron_size) {
		break;
	      }
	   
	      linked_list_iterator_init(short_cal_aux->sr_list, &itr2);
	      s = linked_list_iterator_curr(&itr2);
	      while (s != NULL) {
		if (s->read_end < s_first->read_start) {
		  seed_region = seed_region_new(s->read_start, s->read_end,
						s->genome_start, s->genome_end, 0);	
		  array_list_insert(seed_region, cal->candidates_seeds_start);
		}
		s = linked_list_iterator_next(&itr2);

	      }
	      list_item_aux = list_item_aux->prev;	      
	    }
	  }	

	  if (s_last && s_last->read_end < read_length - seed_size) {
	    // Search start seed ----->
	    list_item_aux = list_item_cal->next;
	    while (list_item_aux) {
	      short_cal_aux = list_item_aux->item;
	      if (!short_cal_aux ||
		  short_cal_aux->end - short_cal_aux->start + 1 >= min_cal_size || 
		  short_cal_aux->start - short_cal->end >= max_intron_size) {
		break;
	      }

	      linked_list_iterator_init(short_cal_aux->sr_list, &itr2);
	      s = linked_list_iterator_curr(&itr2);
	      while (s != NULL) {
		if (s->read_start > s_last->read_end) {
		  seed_region = seed_region_new(s->read_start, s->read_end,
						s->genome_start, s->genome_end, 1);	
		  array_list_insert(seed_region, cal->candidates_seeds_end);
		}
		s = linked_list_iterator_next(&itr2);
	      }
	      list_item_aux = list_item_aux->next;
	    }	    
	  }
	  //======================================================================//
        } 
	linked_list_iterator_next(&itr);
	list_item_cal = linked_list_iterator_list_item_curr(&itr);
      }
    }
  }

  for (unsigned int i = 0; i < nstrands; i++) {
    for (unsigned int j = 0; j < nchromosomes; j++) {
      linked_list_free(cals_list[i][j], short_cal_free);
    }
    free(cals_list[i]);
  }
  
  array_list_free(mapping_list, NULL);
  free(cals_list);


  return array_list_size(cal_list);
  
}
*/


size_t bwt_generate_cals(char *seq, 
			 size_t seed_size, 
			 bwt_optarg_t *bwt_optarg,
			 cal_optarg_t *cal_optarg,
			 bwt_index_t *index, 
			 array_list_t *cal_list, 
			 unsigned int nchromosomes) {
  nchromosomes += 1;
  //printf("%i\n", nchromosomes);
  //printf("============New function generate CALs================\n");
  //NEW! Function! Fusion between Seeding and Caling
  array_list_t *mapping_list = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  const int nstrands = 2;
  
  linked_list_t ***cals_list = (linked_list_t ***)malloc(sizeof(linked_list_t **)*nstrands);

  for (unsigned int i = 0; i < nstrands; i++) {
    cals_list[i] = (linked_list_t **)malloc(sizeof(linked_list_t *)*nchromosomes);
    for (unsigned int j = 0; j < nchromosomes; j++) {
      cals_list[i][j] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    }
  }

  size_t len = strlen(seq);
  int seed_id = 0;

  size_t offset, num_seeds;
  
  //printf(" len=%i, seed_size=%i, num_seeds=%i, seq=%s\n", len, seed_size, num_seeds, seq);
  uint8_t *code_seq = (uint8_t *) malloc(len * sizeof(uint8_t));

  encode_bases(code_seq, seq, len, index->bwt_config.table);

  //Second first seed
  if (seed_size <= 10) { seed_size = 16; }

  int padding_left = seed_size / 2;
  int padding_right = padding_left;

  int extra_seed_size = 16;//seed_size;
  //int first_seed = 16;

  //Extra seed for splice junctions
  bwt_map_exact_seed(code_seq, len, padding_left, padding_left + extra_seed_size - 1,
  		     bwt_optarg, index, mapping_list, seed_id++);

  insert_seeds_and_merge(mapping_list, cals_list,  len);
  num_seeds = len / seed_size;

  //Extra seed for splice junctions
    // first 'pasada'
  offset = 0;
  for (size_t i = 0; i < num_seeds; i++) {
    //printf("SEED %i-%i\n", offset, offset + seed_size - 1);
    bwt_map_exact_seed(code_seq, len, offset, offset + seed_size - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    //printf("\tMapping list %i\n", array_list_size(mapping_list));
    insert_seeds_and_merge(mapping_list, cals_list,  len);
    offset += seed_size;
  }
  
  //Last seed
  if (len % seed_size > 0) {
    bwt_map_exact_seed(code_seq, len, len - seed_size, len - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    insert_seeds_and_merge(mapping_list, cals_list,  len);
  }

  if (len % seed_size != padding_right) {
    //Extra seed for splice junctions
    bwt_map_exact_seed(code_seq, len, len - extra_seed_size - padding_right, len - padding_right - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    insert_seeds_and_merge(mapping_list, cals_list,  len);
  }

  free(code_seq);

  //Store CALs in Array List for return results                                                                                                   
  size_t start_cal, end_cal;
  size_t seq_start, seq_end;
  seed_region_t *s, *s_first, *s_last, *s_aux, *seed_region;
  cal_t *cal;
  short_cal_t *short_cal, *short_cal_aux;
  linked_list_iterator_t itr, itr2;
  linked_list_item_t *list_item_cal, *list_item_aux;

  size_t min_cal_size = cal_optarg->min_cal_size;
  const int max_intron_size = cal_optarg->max_intron_size;

  
  //if (min_cal_size != 20) { printf("%i\n", min_cal_size); exit(-1); }
  //printf("%i - %i\n", min_cal_size, max_intron_size);
  size_t read_length = len;
  //min_cal_size = 20;
  //cal_optarg->min_cal_size;// TODO:PARAMETER
  //const int max_intron_size = 500000; //TODO: Parameter  

  for (unsigned int j = 0; j < nchromosomes; j++) {
    for (unsigned int i = 0; i < nstrands; i++) {
      linked_list_iterator_init(cals_list[i][j], &itr);
      list_item_cal = linked_list_iterator_list_item_curr(&itr);
      while ((list_item_cal != NULL )) {
	short_cal = (short_cal_t *)list_item_cal->item;
	//printf("%i:%i-%i\n", j + 1, short_cal->start, short_cal->end);
	if (short_cal->end - short_cal->start + 1 >= min_cal_size) {
	  //printf("Yes\n");
	  linked_list_t *list_aux = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_list)) {
	    //TODO: Change all parameters to seed_region_t
	    append_seed_region_linked_list(list_aux,
					   s->read_start, s->read_end, 
					   s->genome_start, s->genome_end, 
					   s->id, 0, 0, 0);	    
	    seed_region_free(s);	    
	  }
	  cal = cal_new(j, i, short_cal->start, 
			short_cal->end, short_cal->num_seeds, 
			list_aux, short_cal->sr_duplicate_list);
 	  array_list_insert(cal, cal_list);
	  
	  short_cal->sr_duplicate_list = NULL;
	  
	  //============= Search if one seed map near to this CAL ================//
	  s_first = (seed_region_t *)linked_list_get_first(list_aux);	  
	  s_last = (seed_region_t *)linked_list_get_last(list_aux);
	  
	  //printf("%i > %i", s_first->read_start, seed_size);
	  if (s_first && s_first->read_start >= seed_size) {
	    //Search start seed <-----
	    //printf("Search <-------");
	    list_item_aux = list_item_cal->prev;
	    while (list_item_aux) {
	      short_cal_aux = list_item_aux->item;
	      if (!short_cal_aux ||
		  short_cal_aux->end - short_cal_aux->start + 1 >= min_cal_size || 
		  short_cal->start - short_cal_aux->end >= max_intron_size) {
		break;
	      }
	   
	      linked_list_iterator_init(short_cal_aux->sr_list, &itr2);
	      s = linked_list_iterator_curr(&itr2);
	      while (s != NULL) {
		if (s->read_end < s_first->read_start) {
		  seed_region = seed_region_new(s->read_start, s->read_end,
						s->genome_start, s->genome_end, 0, 0, 0);	
		  array_list_insert(seed_region, cal->candidates_seeds_start);
		}
		s = linked_list_iterator_next(&itr2);

	      }
	      list_item_aux = list_item_aux->prev;	      
	    }
	  }	

	  if (s_last && s_last->read_end < read_length - seed_size) {
	    // Search start seed ----->
	    list_item_aux = list_item_cal->next;
	    while (list_item_aux) {
	      short_cal_aux = list_item_aux->item;
	      if (!short_cal_aux ||
		  short_cal_aux->end - short_cal_aux->start + 1 >= min_cal_size || 
		  short_cal_aux->start - short_cal->end >= max_intron_size) {
		break;
	      }

	      linked_list_iterator_init(short_cal_aux->sr_list, &itr2);
	      s = linked_list_iterator_curr(&itr2);
	      while (s != NULL) {
		if (s->read_start > s_last->read_end) {
		  seed_region = seed_region_new(s->read_start, s->read_end,
						s->genome_start, s->genome_end, 1, 0, 0);	
		  array_list_insert(seed_region, cal->candidates_seeds_end);
		}
		s = linked_list_iterator_next(&itr2);
	      }
	      list_item_aux = list_item_aux->next;
	    }	    
	  }
	  //======================================================================//
        } 
	linked_list_iterator_next(&itr);
	list_item_cal = linked_list_iterator_list_item_curr(&itr);
      }
    }
  }

  for (unsigned int i = 0; i < nstrands; i++) {
    for (unsigned int j = 0; j < nchromosomes; j++) {
      linked_list_free(cals_list[i][j], (void *)short_cal_free);
    }
    free(cals_list[i]);
  }
  
  array_list_free(mapping_list, NULL);
  free(cals_list);


  return array_list_size(cal_list);
  
}

size_t bwt_map_exact_seeds_by_region(int start_position, int end_position, 
				     char *seq, int seed_size, int min_seed_size,
				     bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
				     array_list_t *mapping_list) {
  //printf("%s\n", seq);
  int len = end_position - start_position;
  int seed_id = 1;
  if (seed_size <= 10) { seed_size = 16; }

  int offset, num_seeds = len / seed_size;
  int extra_seed_offset = seed_size / 2;
  int extra_seed_size = 16;


  uint8_t *code_seq = (uint8_t *) malloc((strlen(seq) + 1) * sizeof(uint8_t));
  encode_bases(code_seq, seq, strlen(seq), index->bwt_config.table);
  /*
  printf("len = %i\n",len);

  for (uintmax_t i_=0; i_<((uintmax_t) (len)); i_++) { 
    printf("%ju vs %c,", (uintmax_t) (code_seq)[i_], seq[i_]);
  }
  printf("\n");
*/
  num_seeds = len / seed_size;

  offset = start_position;
  char aux_seq[1000];
  for (size_t i = 0; i < num_seeds; i++) {
    //memcpy(aux_seq, seq + offset, seed_size); aux_seq[seed_size] = '\0';
    //LOG_DEBUG_F("\t Seed [%i-%i]: %s\n", start_position + offset, start_position + offset + seed_size - 1, aux_seq);

    bwt_map_exact_seed(code_seq, strlen(seq), offset, offset + seed_size - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    offset += seed_size;
  }
  
  if (len % seed_size >  0) {
    //memcpy(aux_seq, seq + (end_position - seed_size), seed_size); aux_seq[seed_size] = '\0';
    //LOG_DEBUG_F("\t Seed [%i-%i]: %s\n", end_position - seed_size, end_position - 1, aux_seq);
    bwt_map_exact_seed(code_seq, strlen(seq), end_position - seed_size, end_position - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
  }

  free(code_seq);

  return array_list_size(mapping_list);
}

size_t bwt_generate_cals_between_coords(int strand_target, int chromosome_target,
					size_t start_target, size_t end_target, 
					int start_position, int end_position, 
					char *seq, int seed_size, int min_seed_size,
					bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
					array_list_t *init_list, array_list_t *cal_list) {

  //printf("SEEDS CAL BETWEEN %lu - %lu\n", start_position, end_position);

  //NEW! Function! Fusion between Seeding and Caling
  array_list_t *mapping_list = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  const int nstrands = 2;
  const int nchromosomes = 30; //TODO: Parameter

  linked_list_t ***cals_list = (linked_list_t ***)malloc(sizeof(linked_list_t **)*nstrands);

  for (unsigned int i = 0; i < nstrands; i++) {
    cals_list[i] = (linked_list_t **)malloc(sizeof(linked_list_t *)*nchromosomes);
    for (unsigned int j = 0; j < nchromosomes; j++) {
      cals_list[i][j] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    }
  }
  
  //TODO: ¿¿We need filter mappings?? ¿¿Call bwt_map_exact_seed_with_filter??
  if (seed_size <= 10) { seed_size = 16; }

  int len = end_position - start_position;
  int seed_id = 1;
  int offset, num_seeds = len / seed_size;
  int extra_seed_offset = seed_size / 2;
  int extra_seed_size = 16;


  uint8_t *code_seq = (uint8_t *) malloc(len * sizeof(uint8_t));
  encode_bases(code_seq, seq, len, index->bwt_config.table);

  num_seeds = len / seed_size;

  insert_seeds_and_merge(init_list, cals_list,  strlen(seq));

  offset = 0;
  char aux_seq[1000];
  printf("Nt without process %i\n", len % seed_size);
  for (size_t i = 0; i < num_seeds; i++) {
    memcpy(aux_seq, seq + start_position + offset, seed_size); aux_seq[seed_size] = '\0';
    LOG_DEBUG_F("\t Seed [%i-%i]: %s\n", start_position + offset, start_position + offset + seed_size - 1, aux_seq);

    bwt_map_exact_seed_by_region(code_seq, strlen(seq), start_position + offset, start_position + offset + seed_size - 1,
				 bwt_optarg, index, mapping_list, seed_id++, strand_target, chromosome_target,
				 start_target, end_target);

    insert_seeds_and_merge(mapping_list, cals_list,  strlen(seq));
    offset += seed_size;
  }
  
  if (len % seed_size >  0) {
    memcpy(aux_seq, seq + (end_position - seed_size), seed_size); aux_seq[seed_size] = '\0';
    LOG_DEBUG_F("\t Seed [%i-%i]: %s\n", end_position - seed_size, end_position - 1, aux_seq);

    bwt_map_exact_seed_by_region(code_seq, strlen(seq), end_position - seed_size, end_position - 1,
				 bwt_optarg, index, mapping_list, seed_id++, strand_target, chromosome_target,
				 start_target, end_target);

    insert_seeds_and_merge(mapping_list, cals_list,  strlen(seq));        
  }

  //Store CALs in Array List for return results                            
  size_t start_cal, end_cal;
  size_t seq_start, seq_end;
  seed_region_t *s, *s_first, *s_last, *s_aux, *seed_region;
  cal_t *cal;
  short_cal_t *short_cal;
  linked_list_iterator_t itr;
  linked_list_item_t *list_item_cal;
  size_t min_cal_size = 16;//seed_size;//cal_optarg->min_cal_size; TODO:PARAMETER

  for (unsigned int j = 0; j < nchromosomes; j++) {
    for (unsigned int i = 0; i < nstrands; i++) {
      linked_list_iterator_init(cals_list[i][j], &itr);
      list_item_cal = linked_list_iterator_list_item_curr(&itr);
      while ((list_item_cal != NULL )) {
	short_cal = (short_cal_t *)list_item_cal->item;
	//printf("Short CAL %i:%lu-%lu \n", short_cal->start, short_cal->end);
	if (short_cal->end - short_cal->start + 1 >= min_cal_size) {
	  linked_list_t *list_aux = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_list)) {
	    printf("\t Seed Region [%lu|%i-%i|%lu] %i\n", s->genome_start, 
		   s->read_start, s->read_end, s->genome_end, s->id);
	    //TODO: Change all parameters to seed_region_t
	    append_seed_region_linked_list(list_aux,
					   s->read_start, s->read_end, 
					   s->genome_start, s->genome_end, 
					   s->id, 0, 0, 0);  
	    seed_region_free(s);
	  }
	  cal = cal_new(j, i, short_cal->start, 
			short_cal->end, short_cal->num_seeds, 
			list_aux, short_cal->sr_duplicate_list);

 	  array_list_insert(cal, cal_list);
	  short_cal->sr_duplicate_list = NULL;
        }
	linked_list_iterator_next(&itr);
	list_item_cal = linked_list_iterator_list_item_curr(&itr);
      }
    }
  }

  if (array_list_size(cal_list) == 0) {
    assert(array_list_size(cal_list) != 0);
  }

  for (unsigned int i = 0; i < nstrands; i++) {
    for (unsigned int j = 0; j < nchromosomes; j++) {
      linked_list_free(cals_list[i][j], (void *)short_cal_free);
    }
    free(cals_list[i]);
  }

  free(cals_list);
  free(code_seq);
  
  //printf("Partial Seeding End\n");
  //return array_list_size(mapping_list);
  return array_list_size(cal_list);
}

//-----------------------------------------------------------------------------
// 11 -> 9 seeds (first pass), and 22 (second pass)
inline size_t seedingOK(char *code_seq, size_t seq_len, size_t num_seeds,
		      size_t seed_size, size_t min_seed_size,
		      bwt_optarg_t *bwt_optarg, bwt_index_t *index,
		      array_list_t *mapping_list) {

  size_t num_mappings = 0;
  size_t offset, offset_inc, offset_end = seq_len - min_seed_size;
  int seed_id = 0;

  if (seed_size * num_seeds > seq_len) {
    offset_inc = seed_size - (((seed_size * num_seeds) - seq_len) / num_seeds);
  } else {
    offset_inc = seq_len / num_seeds;
  }

  for (offset = 0; offset < offset_end; offset += offset_inc) {
    //printf("\nseed [%i - %i]\n", offset, offset + seed_size - 1);
    num_mappings += bwt_map_exact_seed(code_seq, seq_len, offset, offset + seed_size - 1,
				       bwt_optarg, index, mapping_list, seed_id);
    seed_id++;
  }

  return num_mappings;
}

//-----------------------------------------------------------------------------

inline size_t seeding(char *code_seq, size_t seq_len, size_t num_seeds,
		      size_t seed_size, size_t min_seed_size,
		      bwt_optarg_t *bwt_optarg, bwt_index_t *index,
		      array_list_t *mapping_list) {

  size_t n_seeds, total_mappings = 0, num_mappings = 0;
  //  size_t offset, offset_inc, offset_end = seq_len - min_seed_size;
  size_t start, end;
  size_t offset = 0, offset_inc, offset_end = seq_len - min_seed_size;
  int seed_id = 0;

  n_seeds = num_seeds;
  offset_inc = ceil(1.0f * seq_len / (num_seeds + 1));
  if (offset_inc <= 0) offset_inc = 1;
  /*
  if (seed_size * num_seeds > seq_len) {
    size_t max_seeds = seq_len - min_seed_size;
    if (num_seeds >= max_seeds) {
      n_seeds = max_seeds;
      offset_inc = 1;
    } else {
      n_seeds = num_seeds;
      offset_inc = seq_len / num_seeds;
    }
  } else {
    n_seeds = num_seeds;
    offset_inc = seq_len / num_seeds;
  }
  */

  start = 0;
  for (size_t i = 0; i < n_seeds; i++) {
    end = start + seed_size;
    if (end >= seq_len) end = seq_len;
    num_mappings = bwt_map_exact_seed(code_seq, seq_len, start, end - 1,
				      bwt_optarg, index, mapping_list, seed_id);
    seed_id++;
    total_mappings += num_mappings;
    //    LOG_DEBUG_F("\tseed %i\t[%i - %i], length read = %i, num. mappings = %i\n", 
    //		i + 1, start, end, seq_len, num_mappings);
    start += offset_inc;
    if (start > offset_end) {
      if (offset_inc == 1) break;
      start = offset_inc / 2;
    }
    /*
    if (start > offset_end) {
      offset++;
      start = offset;
    }
    */
  }

  //  LOG_DEBUG_F("\t\ttotal mappings = %i\n", total_mappings);

  return total_mappings;
}

//-----------------------------------------------------------------------------

/* size_t bwt_map_exact_seeds_seq_by_num(char *seq, size_t num_seeds,  */
/* 				      size_t seed_size, size_t min_seed_size, */
/* 				      bwt_optarg_t *bwt_optarg, bwt_index_t *index, */
/* 				      array_list_t *mapping_list) { */
/*   size_t seq_len = strlen(seq); */
/*   size_t num_mappings = 0; */
/*   if (seed_size <= 10) { seed_size = 16; } */
/*   char *code_seq = (char *) calloc(seq_len, sizeof(char)); */

/*   replaceBases(seq, code_seq, seq_len); */

/*   num_mappings = seeding(code_seq, seq_len, num_seeds, seed_size, min_seed_size, */
/* 			 bwt_optarg, index, mapping_list); */
/*   //  printf("\tfirst, num_mappings = %d\n", num_mappings); */
/*   /\* */
/*   if (num_mappings < 10) { */
/*     num_mappings += seeding(code_seq, seq_len, max_num_seeds, seed_size - 4, min_seed_size - 4, */
/* 			    bwt_optarg, index, mapping_list); */
/*     //    printf("\tsecond -4, num_mappings (include first) = %d\n", num_mappings); */
/*     //  } else if (num_mappings >= 10000) { */
/*   } else if (num_mappings >= bwt_optarg->filter_read_mappings) { */
/*     array_list_clear(mapping_list, (void *) region_bwt_free); */
/*     num_mappings = seeding(code_seq, seq_len, max_num_seeds, seed_size + 2, min_seed_size + 2, */
/* 			   bwt_optarg, index, mapping_list);  */
/*     //    printf("\tthird +2, num_mappings = %d\n", num_mappings); */
/*  } */
/*   *\/ */
/*   free(code_seq); */
/*   return num_mappings; */
/* } */

//-----------------------------------------------------------------------------
// CAL functions
//-----------------------------------------------------------------------------

int print_short_cal(void *item, void *dummy){
        short_cal_t *coordenate_p = (short_cal_t *)item;
	printf("[%lu-%lu]-> ",  coordenate_p->start, coordenate_p->end);

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

void print_se_region(short_cal_t *short_cal){
  printf("(%lu-%lu)->", short_cal->start, short_cal->end);
}


void print_se_region_cp(short_cal_t *short_cal, void* dummy){
  printf("(%lu-%lu)->", short_cal->start, short_cal->end);
}



void seed_region_select_linked_list(linked_list_t* sr_list, linked_list_t* sr_duplicate_list, 
				    size_t read_start, size_t read_end,
				    size_t genome_start, size_t genome_end,
				    int seed_id, int pos_err, int type_err,
				    unsigned char *seeds_ids_array) {
  //printf("\tInsert [Seed:=%lu-%lu](%i): ", read_start, read_end, seed_id);
  seed_region_t *item;
  if (!seeds_ids_array[seed_id]) { 
    //printf(" Not in!!\n");
    seeds_ids_array[seed_id]++;
    item = seed_region_new(read_start, read_end, genome_start, genome_end, seed_id, pos_err, type_err);
    linked_list_insert(item, sr_list);
  } else {
    //printf(" In!!\n");
    if (seeds_ids_array[seed_id] == 1) { 
      linked_list_iterator_t* itr = linked_list_iterator_new(sr_list);
      item = (seed_region_t *)linked_list_iterator_curr(itr);
      while (item != NULL) {
	if (item->id == seed_id) {
	  item = linked_list_iterator_remove(itr);
	  linked_list_insert(item, sr_duplicate_list);
	  //printf("\tRemove [Seed:=%lu-%lu](%i)\n", item->read_start, item->read_end, item->id);
	  break;
	}
	linked_list_iterator_next(itr);
	item = linked_list_iterator_curr(itr);
      }
      linked_list_iterator_free(itr);
    }
    item = seed_region_new(read_start, read_end, genome_start, genome_end, seed_id, pos_err, type_err);
    linked_list_insert(item, sr_duplicate_list);
    seeds_ids_array[seed_id]++;
  }


  //printf("\t sr_size = %i, sr_duplicate = %i\n", linked_list_size(sr_list), linked_list_size(sr_duplicate_list));

  /*
  unsigned char actualization = 0;
  seed_region_t *item, *item_aux, *new_item, *item_free;
  linked_list_iterator_t* itr = linked_list_iterator_new(list_p);

  printf("\tInsert [Seed:=%lu-%lu]\n", read_start, read_end);
  if (linked_list_size(list_p) <= 0) {
    new_item = seed_region_new(read_start, read_end, genome_start, genome_end);
    linked_list_insert(new_item, list_p);
  } else {
    item = (seed_region_t *)linked_list_iterator_curr(itr);
    while (item != NULL) {
      //printf("\t compare with %lu\n", item->start);
      if (read_start <= item->read_start) {
	new_item = seed_region_new(read_start, read_end, genome_start, genome_end);
	linked_list_iterator_insert(new_item, itr);
	break;
      }      
      //continue loop...
      linked_list_iterator_next(itr);
      item = linked_list_iterator_curr(itr);      
    }// end while

    if (item == NULL) {
      new_item = seed_region_new(read_start, read_end, genome_start, genome_end);
      //printf("\tInsert at END\n");
      linked_list_insert_last(new_item, list_p);
    }
    //printf("Insert OK! and now actualization\n");
  }

  linked_list_iterator_free(itr);
  */
}


//-----------------------------------------------------------------------------

void append_seed_region_linked_list(linked_list_t* sr_list,
				    size_t read_start, size_t read_end, 
				    size_t genome_start, size_t genome_end, 
				    int seed_id, int pos_err, int type_err, int type_seeds) {
  unsigned char actualization = 0;
  seed_region_t *item, *item_aux, *new_item, *item_free;
  linked_list_iterator_t* itr = linked_list_iterator_new(sr_list);
  
  if (linked_list_size(sr_list) <= 0) {
    new_item = seed_region_new(read_start, read_end, genome_start, genome_end, 
			       seed_id, pos_err, type_err);

    if (type_seeds) {
      new_item->errors_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      if (pos_err != -1) {
	bwt_err_t *err = bwt_err_new(pos_err, type_err == SEED_INSERTION ? 'I' : 'D');
	array_list_insert(err, new_item->errors_list);
      }
    }
    linked_list_insert(new_item, sr_list);
    //printf("Call linked_list_insert\n");
  } else {
    item = (seed_region_t *)linked_list_iterator_curr(itr);
    while (item != NULL) {
      //printf("\t compare with %lu\n", item->start);
      if (genome_start < item->genome_start) {
	if (genome_end + 1 < item->genome_start) {
	  /*********************************************
	   *    Case 1: New item insert before item.   *
           *                                           *
           *        new item     item                  *
           *       |-------| |--------|                *
           ********************************************/
	  //printf("\t Insert now before %lu\n", item->start);
	  new_item = seed_region_new(read_start, read_end, genome_start, 
				     genome_end, seed_id, pos_err, type_err);
	  if (type_seeds) {
	    new_item->errors_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	    if (pos_err != -1) {
	      bwt_err_t *err = bwt_err_new(pos_err, type_err == SEED_INSERTION ? 'I' : 'D');
	      array_list_insert(err, new_item->errors_list);
	    }
	  }
	  linked_list_iterator_insert(new_item, itr);
	  //printf("Call linked_list_iterator_insert\n");
	  linked_list_iterator_prev(itr);
	} else {
	  /********************************************
           *  Case 2: Actualization item start        *
           *           new item                       *
           *          |-------|   item                *
           *                   |--------|             *                            
           ********************************************/
	  /*printf("\tFusion Case 2! [%lu|%i - %i|%lu], [%lu|%i - %i|%lu]\n", 
		 genome_start, read_start, genome_end, read_end, 
		 item->genome_start, item->read_start, item->genome_end, item->read_end);*/
	  item->read_start = read_start;
	  item->genome_start = genome_start;
	  if (genome_end > item->genome_end) {
	    /**************************************************
             *  Case 3: Actualization item start and item end *
             *          new item                              *
             *         |------------|                         *
             *              item                              *    
             *           |--------|                           *                                    
             **************************************************/
	    //printf("\tFusion!\n");
	    item->read_end = read_end;
	    item->genome_end = genome_end;
	    actualization = 1;
	  }
	  if (pos_err != -1 && type_seeds) {
	    bwt_err_t *err = bwt_err_new(pos_err, type_err == SEED_INSERTION ? 'I' : 'D');
	    array_list_insert(err, item->errors_list);	  
	  }
	}
	break;
      } else {
	if (genome_end <= item->genome_end) {
	  /**************************************************                                       
           *  Case 4: The new item don't insert in the list *                             
           *              item                              * 
           *         |-------------|                        * 
           *             new item                           * 
           *            |--------|                          * 
           **************************************************/
	  //printf("\tFusion!\n");
	  break;
	} else if (item->genome_end + 1 >= genome_start) {
	  /********************************************                                              
           *  Case 5: Actualization item end          *
           *            item                          *                                              
           *          |-------| new item              *                                            
           *                 |--------|               *                                              
           ********************************************/
	  //printf("\tFusion!\n");
	  /*printf("\tFusion Case 5! New Item:[%lu|%i - %i|%lu], Item:[%lu|%i - %i|%lu]\n", 
		 genome_start, read_start, read_end, genome_end,  
		 item->genome_start, item->read_start, item->read_end, item->genome_end);
	  */
	  item->read_end = read_end;
	  item->genome_end = genome_end;
	  actualization = 1;
	  if (pos_err != -1 && type_seeds) {
	    bwt_err_t *err = bwt_err_new(pos_err, type_err == SEED_INSERTION ? 'I' : 'D');
	    array_list_insert(err, item->errors_list);	  
	  }
	  break;
	}
      } // end else

      //continue loop...
      linked_list_iterator_next(itr);
      item = linked_list_iterator_curr(itr);
      
    } // end while

    if (item == NULL) {
      /******************************************************* 
       * Case 5: Insert new item at the end of the list      * 
       *                 item    new item                    * 
       *              |-------| |--------|                   *    
       *******************************************************/
      new_item = seed_region_new(read_start, read_end, genome_start, 
				 genome_end, seed_id, pos_err, type_err);
      if (type_seeds) {
	new_item->errors_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	if (pos_err != -1) {
	  bwt_err_t *err = bwt_err_new(pos_err, type_err == SEED_INSERTION ? 'I' : 'D');
	  array_list_insert(err, new_item->errors_list);
	}
      }
      //printf("\tInsert at END\n");
      linked_list_insert_last(new_item, sr_list);
      //printf("Call linked_list_insert_last\n");
    }
    //printf("Insert OK! and now actualization\n");
    if (actualization == 1) {
      //printf("\tActualization RIGHT items (Next). Current item [%d-%d]\n", item->start, item->end);
      //printf("List before actualization\n");
      //linked_list_print(list_p, print_se_region);
      //printf("\n\n");
 
      linked_list_iterator_next(itr);
      item_aux = linked_list_iterator_curr(itr);

      while (item_aux != NULL) {
	//printf("\t\tFusion right items. item->end=%d < item_aux->start=%d?\n", item->end, item_aux->start);
	//printf("\tIterator are in [%lu-%lu] and compare with item [%lu-%lu]\n", item_aux->start, item_aux->end, item->start, item->end);
	if (item->genome_end + 1 < item_aux->genome_start) {
	  //printf("\t\tSTOP Actualization\n");
	  break;
	} else {
	  //printf("\t\tCONTINUE Actualization. item->end=%d < item_aux->end=%d?\n", item->end, item_aux->end);
	  if (item->genome_end + 1 < item_aux->genome_end) {
	    //printf("\t\tActualization end value %d\n", item_aux->end);
	    item->read_end = item_aux->read_end;
	    item->genome_end = item_aux->genome_end;
	    if (type_seeds && item_aux->errors_list) {
	      for (int i = 0; i < array_list_size(item_aux->errors_list); i++) {
		array_list_insert(array_list_get(i, item_aux->errors_list), item->errors_list);
	      }
	    }
	  }
          //printf("\t\tDelete item %d-%d\n", item_aux->start, item_aux->end);
	  item_free = linked_list_iterator_remove(itr);
	  if (item_free) { 
	    //if (type_seeds && item_free->errors_list) { array_list_free(item_free->errors_list, (void *)NULL); }
	    seed_region_free(item_free);
	  }
	  //printf("\t\tDelete OK!\n");
	  item_aux = linked_list_iterator_curr(itr);
	}                                                                                             
      }
      /*printf("List after actualization\n");
      linked_list_print(list_p, print_se_region);
      printf("\n\n");*/
    }
  }//end first else
  //printf("End insert and actualization\n");
  linked_list_iterator_free(itr);

  /*printf("Status Seed Region list %lu [%i-%i][%i-%i]:\n", linked_list_size(sr_list), 
	 ((seed_region_t *)sr_list->first->item)->read_start, ((seed_region_t *)sr_list->first->item)->read_end,
	 ((seed_region_t *)sr_list->last->item)->read_start,  ((seed_region_t *)sr_list->last->item)->read_end);
  */
  /*for (linked_list_item_t *list_item = sr_list->first; list_item != NULL; list_item = list_item->next) {
    printf("Status Seed Region list %lu:\n", linked_list_size(sr_list));
    seed_region_t *s = list_item->item;
    printf("[%i|%i - %i|%i]  ", s->genome_start, s->read_start, s->read_end, s->genome_end);
  }
  printf("\n");
  */
}



//-----------------------------------------------------------------------------

void my_cp_list_append_linked_list(linked_list_t* list_p, region_t *region, size_t max_cal_distance, int max_seeds) {  
  unsigned char actualization = 0;
  short_cal_t *item, *item_aux, *new_item_p, *item_free;
  
  //int strand = region->strand;
  size_t start = region->start;
  size_t end = region->end;
  size_t seq_start = region->seq_start;
  size_t seq_end = region->seq_end;
  size_t seq_len = region->seq_len;
  linked_list_iterator_t* itr = linked_list_iterator_new(list_p);
  //printf("LINKED LIST INSERT %lu|%i-%i|%lu:\n", start, seq_start, seq_end, end);

  if (linked_list_size(list_p) <= 0) {
    new_item_p = short_cal_new(start, end, seq_start, seq_end, 
			       seq_len, max_seeds, region->id,
			       region->pos_err,
			       region->type_err);
    linked_list_insert(new_item_p, list_p);
  } else {
    item = (short_cal_t *)linked_list_iterator_curr(itr);
    while (item != NULL) {
      if (start < item->start) {
	if (end + max_cal_distance < item->start) {
	  /*********************************************
	   *    Case 1: New item insert before item.   *
           *                                           *
           *        new item     item                  *
           *       |-------| |--------|                *
           ********************************************/

	  //printf("\t Insert now before %lu\n", item->start);
	  new_item_p = short_cal_new(start, end, seq_start, seq_end,
				     seq_len, max_seeds, region->id,
				     region->pos_err,
				     region->type_err);

	  linked_list_iterator_insert(new_item_p, itr);
	  linked_list_iterator_prev(itr);
	} else {
	  /********************************************
           *  Case 2: Actualization item start        *
           *           new item                       *
           *          |-------|   item                *
           *                   |--------|             *                            
           ********************************************/
	  //printf("\tFusion!\n");
	  item->start = start;
	  item->num_seeds++;
	  //item->seq_start = seq_start;
	  seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list, 
					 seq_start, seq_end, start, end, region->id,
					 region->pos_err,
					 region->type_err,
					 item->seeds_ids_array);
	  if (end > item->end) {
	    /**************************************************
             *  Case 3: Actualization item start and item end *
             *          new item                              *
             *         |------------|                         *
             *              item                              *    
             *           |--------|                           *                                    
             **************************************************/
	    //printf("\tFusion!\n");
	    item->end = end;
	    //item->seq_end = seq_end;
	    actualization = 1;
	  }
	}
	break;
      } else {
	if (end <= item->end) {
	  /**************************************************                                       
           *  Case 4: The new item don't insert in the list *                             
           *              item                              * 
           *         |-------------|                        * 
           *             new item                           * 
           *            |--------|                          * 
           **************************************************/
	  //printf("\tFusion!\n");
	  item->num_seeds++;
	  seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list, seq_start,
					 seq_end, start, end, region->id, 
					 region->pos_err,
					 region->type_err,
					 item->seeds_ids_array);
	  break;
	} else if (item->end + max_cal_distance >= start) {
	  /********************************************                                              
           *  Case 5: Actualization item end          *
           *            item                          *                                              
           *          |-------| new item              *                                            
           *                 |--------|               *                                              
           ********************************************/
	  //printf("\tFusion!\n");
	  item->end = end;
	  //item->seq_end = seq_end;
	  seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list, seq_start,
					 seq_end, start, end, region->id,
					 region->pos_err,
					 region->type_err,
					 item->seeds_ids_array);
	  actualization = 1;
	  item->num_seeds++;
	  break;
	}
      } // end else

      //continue loop...
      linked_list_iterator_next(itr);
      item = linked_list_iterator_curr(itr);
      
    } // end while

    if (item == NULL) {
      /******************************************************* 
       * Case 5: Insert new item at the end of the list      * 
       *                 item    new item                    * 
       *              |-------| |--------|                   *    
       *******************************************************/
      new_item_p = short_cal_new(start, end, seq_start, seq_end, seq_len, max_seeds, region->id,
				 region->pos_err, region->type_err);
      //printf("\tInsert at END\n");
      linked_list_insert_last(new_item_p, list_p);
    }
    //printf("Insert OK! and now actualization\n");
    if (actualization == 1) {
      //printf("\tActualization RIGHT items (Next). Current item [%d-%d]\n", item->start, item->end);
      //printf("List before actualization\n");
      //linked_list_print(list_p, print_se_region);
      //printf("\n\n");
 
      linked_list_iterator_next(itr);
      item_aux = linked_list_iterator_curr(itr);

      while (item_aux != NULL) {
	//printf("\t\tFusion right items. item->end=%d < item_aux->start=%d?\n", item->end, item_aux->start);
	//printf("\tIterator are in [%lu-%lu] and compare with item [%lu-%lu]\n", item_aux->start, item_aux->end, item->start, item->end);
	if (item->end + max_cal_distance < item_aux->start) {
	  //printf("\t\tSTOP Actualization\n");
	  break;
	} else {
	  //printf("\t\tCONTINUE Actualization. item->end=%d < item_aux->end=%d?\n", item->end, item_aux->end);
	  if (item->end < item_aux->end) {
	    //printf("\t\tActualization end value %d\n", item_aux->end);
	    item->end = item_aux->end;
	    //item->seq_end = item_aux->seq_end;	    
	    seed_region_t *seed_region_aux;
	    while (seed_region_aux = linked_list_remove_first(item_aux->sr_list)) {
	      seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list, seed_region_aux->read_start, 
					     seed_region_aux->read_end, seed_region_aux->genome_start,
					     seed_region_aux->genome_end, seed_region_aux->id, 
					     seed_region_aux->pos_err, seed_region_aux->type_err,
					     item->seeds_ids_array);
	      seed_region_free(seed_region_aux);
	    }

	  }
          //printf("\t\tDelete item %d-%d\n", item_aux->start, item_aux->end);
	  item_free = linked_list_iterator_remove(itr);
	  if (item_free) { short_cal_free(item_free); }
	  //printf("\t\tDelete OK!\n");
	  item_aux = linked_list_iterator_curr(itr);
	}                                                                                             
      }
      /*printf("List after actualization\n");
      linked_list_print(list_p, print_se_region);
      printf("\n\n");*/
    }
  }//end first else
  //printf("End insert and actualization\n");
  linked_list_iterator_free(itr);
}

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
/*
size_t bwt_generate_cal_list_linked_list(array_list_t *mapping_list,
					 cal_optarg_t *cal_optarg,
					 int *min_seeds, int *max_seeds,
					 size_t nchromosomes,
					 array_list_t *cal_list,
					 size_t read_length) {
  //  printf("::: CALS PROCESS with max seeds %i:::\n", *max_seeds);
  short_cal_t *short_cal, *short_cal_aux;
  linked_list_item_t *list_item_cal, *list_item_aux, *list_item_s;
  region_t *region;
  size_t min_cal_size = cal_optarg->min_cal_size;
  size_t max_cal_distance  = read_length;
  size_t num_mappings = array_list_size(mapping_list);
  size_t chromosome_id;
  short int strand;
  size_t start, end;
  int seed_size = cal_optarg->seed_size; //Change parameter
  int min_intron_size = 500000; //Change parameter
  linked_list_iterator_t itr, itr2;
  const unsigned char nstrands = 2;
  linked_list_t ***cals_list = (linked_list_t ***)malloc(sizeof(linked_list_t **)*nstrands);
  nchromosomes += 1;

  for (unsigned int i = 0; i < nstrands; i++) {
    cals_list[i] = (linked_list_t **)malloc(sizeof(linked_list_t *)*nchromosomes);
    for (unsigned int j = 0; j < nchromosomes; j++) {
      cals_list[i][j] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    }
  }

  // from the results mappings, generates CALs
  for (unsigned int m = 0; m < num_mappings; m++) {
    region = array_list_get(m, mapping_list);
    chromosome_id = region->chromosome_id;
    strand = region->strand;
    
    //printf("Inserting [Region(%i):=%i:%lu-%lu] [Seed:=%lu-%lu]\n", strand, chromosome_id, 
    //	   region->start, region->end, region->seq_start, region->seq_end);
    
    my_cp_list_append_linked_list(cals_list[strand][chromosome_id], region, max_cal_distance, *max_seeds);    
  }

  //Store CALs in Array List for return results                                                                                                                        
  size_t start_cal, end_cal, len;
  size_t seq_start, seq_end;
  seed_region_t *s, *s_first, *s_last, *s_aux, *seed_region;
  cal_t *cal;

  for (unsigned int j = 0; j < nchromosomes; j++) {
    for (unsigned int i = 0; i < nstrands; i++) {
      linked_list_iterator_init(cals_list[i][j], &itr);
      list_item_cal = linked_list_iterator_list_item_curr(&itr);
      //short_cal_p = linked_list_iterator_curr(&itr);
      while ((list_item_cal != NULL )) {
	short_cal = (short_cal_t *)list_item_cal->item;
	//printf("Min CAL [%i:%lu-%lu](%i):", j, short_cal->start, short_cal->end, short_cal->end - short_cal->start);
	if (short_cal->end - short_cal->start + 1 >= min_cal_size) {
	  //printf("Min CAL [%i:%lu-%lu](%i):", j, short_cal->start, short_cal->end, short_cal->end - short_cal->start);
	  //printf(" In!\n");
	  linked_list_t *list_aux = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  //	  printf("Removing %i items\n", linked_list_size(short_cal->sr_list));
	  while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_list)) {
	    //TODO: Change all parameters to seed_region_t
	    //printf("Extract[%i|%i - %i|%i] and Insert:\n", s->genome_start, s->read_start, s->read_end, s->genome_end);
	    append_seed_region_linked_list(list_aux,
					   s->read_start, s->read_end, 
					   s->genome_start, s->genome_end, 
					   s->id);	    
	    seed_region_free(s);	    
	  }
	  cal = cal_new(j, i, short_cal->start, 
			short_cal->end, short_cal->num_seeds, 
			list_aux, short_cal->sr_duplicate_list);
 	  array_list_insert(cal, cal_list);

	  short_cal->sr_duplicate_list = NULL;

	  //============= Search if one seed map near to this CAL ================//
	  s_first = (seed_region_t *)linked_list_get_first(list_aux);	  
	  s_last = (seed_region_t *)linked_list_get_last(list_aux);

	  if (s_first && s_first->read_start > seed_size) {
	    //Search start seed <-----
	    list_item_aux = list_item_cal->prev;
	    while (list_item_aux) {
	      short_cal_aux = list_item_aux->item;
	      if (!short_cal_aux ||
		  short_cal_aux->end - short_cal_aux->start + 1 >= min_cal_size || 
		  short_cal->start - short_cal_aux->end >= min_intron_size) {
		break;
	      }
	   
	      linked_list_iterator_init(short_cal_aux->sr_list, &itr2);
	      s = linked_list_iterator_curr(&itr2);
	      while (s != NULL) {
		if (s->read_end < s_first->read_start) {
		  seed_region = seed_region_new(s->read_start, s->read_end,
						s->genome_start, s->genome_end, 0);	
		  array_list_insert(seed_region, cal->candidates_seeds_start);
		}
		s = linked_list_iterator_next(&itr2);

	      }
	      list_item_aux = list_item_aux->prev;	      
	    }
	  }	

	  if (s_last && s_last->read_end < read_length - seed_size) {
	    // Search start seed ----->
	    list_item_aux = list_item_cal->next;
	    while (list_item_aux) {
	      short_cal_aux = list_item_aux->item;
	      if (!short_cal_aux ||
		  short_cal_aux->end - short_cal_aux->start + 1 >= min_cal_size || 
		  short_cal_aux->start - short_cal->end >= min_intron_size) {
		break;
	      }

	      linked_list_iterator_init(short_cal_aux->sr_list, &itr2);
	      s = linked_list_iterator_curr(&itr2);
	      while (s != NULL) {
		if (s->read_start > s_last->read_end) {
		  seed_region = seed_region_new(s->read_start, s->read_end,
						s->genome_start, s->genome_end, 1);	
		  array_list_insert(seed_region, cal->candidates_seeds_end);
		}
		s = linked_list_iterator_next(&itr2);
	      }
	      list_item_aux = list_item_aux->next;
	    }	    
	  }
	  	  
	  //printf("\t First Seed [%i-%i]\n", s_first->read_start, s_first->read_end);
	  //printf("\t Last Seed [%i-%i]\n", s_last->read_start, s_last->read_end);

	  //======================================================================//

        } //else { printf("\n"); }
	linked_list_iterator_next(&itr);
	//linked_list_item_free(list_item_cal, short_cal_free);
	list_item_cal = linked_list_iterator_list_item_curr(&itr);
      }
      //free(cals_list[i][j]);
    }
  }

  for (unsigned int i = 0; i < nstrands; i++) {
    for (unsigned int j = 0; j < nchromosomes; j++) {
      linked_list_free(cals_list[i][j], short_cal_free);
    }
    free(cals_list[i]);
  }
  

  //for (unsigned int i = 0; i < nstrands; i++) {
  //free(cals_list[i]);
  //}

  free(cals_list);

  *min_seeds = 10000;
  *max_seeds = 0;
  for (int i = 0; i < array_list_size(cal_list); i++) {
    cal_t *cal = array_list_get(i, cal_list);
    if (*min_seeds > cal->num_seeds) *min_seeds = cal->num_seeds;
    if (*max_seeds < cal->num_seeds) *max_seeds = cal->num_seeds;
  }
  return array_list_size(cal_list);
}
*/
//-----------------------------------------------------------------------------

size_t bwt_generate_cal_list_linked_list(array_list_t *mapping_list,
					 cal_optarg_t *cal_optarg,
					 int *min_seeds,
					 int *max_seeds,
					 size_t nchromosomes,
					 array_list_t *cal_list,
					 size_t read_length,
					 size_t min_cal_size,
					 int type_seeds) {
  //  printf("::: CALS PROCESS with max seeds %i:::\n", *max_seeds);
  short_cal_t *short_cal, *short_cal_aux;
  linked_list_item_t *list_item_cal, *list_item_aux, *list_item_s;
  region_t *region;

  //size_t min_cal_size = cal_optarg->min_cal_size;
  int seed_size = cal_optarg->seed_size; //Change parameter
  const int max_intron_size = cal_optarg->max_intron_size;  

  size_t max_cal_distance  = read_length;
  size_t num_mappings = array_list_size(mapping_list);
  size_t chromosome_id;
  short int strand;
  size_t start, end;
  //if (min_cal_size != 20) { printf("%i\n", min_cal_size); exit(-1); }

  //int max_intron_size = 500000; //Change parameter
  linked_list_iterator_t itr, itr2;
  const unsigned char nstrands = 2;
  linked_list_t ***cals_list = (linked_list_t ***)malloc(sizeof(linked_list_t **)*nstrands);
  nchromosomes += 1;

  for (unsigned int i = 0; i < nstrands; i++) {
    cals_list[i] = (linked_list_t **)malloc(sizeof(linked_list_t *)*nchromosomes);
    for (unsigned int j = 0; j < nchromosomes; j++) {
      cals_list[i][j] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    }
  }

  // from the results mappings, generates CALs
  for (unsigned int m = 0; m < num_mappings; m++) {
    region = array_list_get(m, mapping_list);
    chromosome_id = region->chromosome_id;
    strand = region->strand;
    
    //printf("Inserting [Region(%i):=%i:%lu-%lu] [Seed:=%lu-%lu]\n", strand, chromosome_id, 
    //	   region->start, region->end, region->seq_start, region->seq_end);
    
    my_cp_list_append_linked_list(cals_list[strand][chromosome_id], region, max_cal_distance, *max_seeds);    
  }

  //Store CALs in Array List for return results                                                                                                                        
  size_t start_cal, end_cal, len;
  size_t seq_start, seq_end;
  seed_region_t *s, *s_first, *s_last, *s_aux, *seed_region;
  cal_t *cal;

  for (unsigned int j = 0; j < nchromosomes; j++) {
    for (unsigned int i = 0; i < nstrands; i++) {
      linked_list_iterator_init(cals_list[i][j], &itr);
      list_item_cal = linked_list_iterator_list_item_curr(&itr);
      //short_cal_p = linked_list_iterator_curr(&itr);
      while ((list_item_cal != NULL )) {
	short_cal = (short_cal_t *)list_item_cal->item;
	//printf("Min CAL [%i:%lu-%lu](%i):", j, short_cal->start, short_cal->end, short_cal->end - short_cal->start);
	if (short_cal->end - short_cal->start + 1 >= min_cal_size) {
	  //printf("Min CAL [%i:%lu-%lu](%i):", j, short_cal->start, short_cal->end, short_cal->end - short_cal->start);
	  //printf(" In!\n");
	  linked_list_t *list_aux = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  //	  printf("Removing %i items\n", linked_list_size(short_cal->sr_list));
	  while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_list)) {
	    //TODO: Change all parameters to seed_region_t
	    //printf("Extract[%i|%i - %i|%i] and Insert:\n", s->genome_start, s->read_start, s->read_end, s->genome_end);
	    append_seed_region_linked_list(list_aux,
					   s->read_start, s->read_end, 
					   s->genome_start, s->genome_end, 
					   s->id, s->pos_err, s->type_err, type_seeds);	    
	    seed_region_free(s);	    
	  }
	  cal = cal_new(j, i, short_cal->start, 
			short_cal->end, short_cal->num_seeds, 
			list_aux, short_cal->sr_duplicate_list);
	  cal->type_seeds = type_seeds;
 	  array_list_insert(cal, cal_list);

	  short_cal->sr_duplicate_list = NULL;

	  if (!type_seeds) {
	    //============= Search if one seed map near to this CAL ================//
	    s_first = (seed_region_t *)linked_list_get_first(list_aux);	  
	    s_last = (seed_region_t *)linked_list_get_last(list_aux);

	    if (s_first && s_first->read_start >= seed_size) {
	      //Search start seed <-----
	      list_item_aux = list_item_cal->prev;
	      while (list_item_aux) {
		short_cal_aux = list_item_aux->item;
		if (!short_cal_aux ||
		    short_cal_aux->end - short_cal_aux->start + 1 >= min_cal_size || 
		    short_cal->start - short_cal_aux->end >= max_intron_size) {
		  break;
		}
	   
		linked_list_iterator_init(short_cal_aux->sr_list, &itr2);
		s = linked_list_iterator_curr(&itr2);
		while (s != NULL) {
		  if (s->read_end < s_first->read_start) {
		    seed_region = seed_region_new(s->read_start, s->read_end,
						  s->genome_start, s->genome_end, 0, 0, 0);	
		    array_list_insert(seed_region, cal->candidates_seeds_start);
		  }
		  s = linked_list_iterator_next(&itr2);

		}
		list_item_aux = list_item_aux->prev;	      
	      }
	    }	

	    if (s_last && s_last->read_end < read_length - seed_size) {
	      // Search start seed ----->
	      list_item_aux = list_item_cal->next;
	      while (list_item_aux) {
		short_cal_aux = list_item_aux->item;
		if (!short_cal_aux ||
		    short_cal_aux->end - short_cal_aux->start + 1 >= min_cal_size || 
		    short_cal_aux->start - short_cal->end >= max_intron_size) {
		  break;
		}

		linked_list_iterator_init(short_cal_aux->sr_list, &itr2);
		s = linked_list_iterator_curr(&itr2);
		while (s != NULL) {
		  if (s->read_start > s_last->read_end) {
		    seed_region = seed_region_new(s->read_start, s->read_end,
						  s->genome_start, s->genome_end, 1, 0, 0);	
		    array_list_insert(seed_region, cal->candidates_seeds_end);
		  }
		  s = linked_list_iterator_next(&itr2);
		}
		list_item_aux = list_item_aux->next;
	      }	    
	    }
	  	  
	    //printf("\t First Seed [%i-%i]\n", s_first->read_start, s_first->read_end);
	    //printf("\t Last Seed [%i-%i]\n", s_last->read_start, s_last->read_end);
	    //======================================================================//
	  }
        } //else { printf("\n"); }
	linked_list_iterator_next(&itr);
	//linked_list_item_free(list_item_cal, short_cal_free);
	list_item_cal = linked_list_iterator_list_item_curr(&itr);
      }
      //free(cals_list[i][j]);
    }
  }

  for (unsigned int i = 0; i < nstrands; i++) {
    for (unsigned int j = 0; j < nchromosomes; j++) {
      linked_list_free(cals_list[i][j], (void *) short_cal_free);
    }
    free(cals_list[i]);
  }
  

  //for (unsigned int i = 0; i < nstrands; i++) {
  //free(cals_list[i]);
  //}

  free(cals_list);

  *min_seeds = 10000;
  *max_seeds = 0;
  for (int i = 0; i < array_list_size(cal_list); i++) {
    cal_t *cal = array_list_get(i, cal_list);
    if (*min_seeds > cal->num_seeds) *min_seeds = cal->num_seeds;
    if (*max_seeds < cal->num_seeds) *max_seeds = cal->num_seeds;
  }

    return array_list_size(cal_list);
}

//-----------------------------------------------------------------------------

void order_cal_linked_list(int num_order, array_list_t *cal_list) {
  int cal_list_size = array_list_size(cal_list);
  
  for (int i = 0; i < num_order; i++) {
    cal_t *cal_prev = array_list_get(i, cal_list);
    for (int j = i + 1; j < cal_list_size; j++) {
      cal_t *cal_next = array_list_get(j, cal_list);
      if (cal_next->read_area > cal_prev->read_area) {
	array_list_swap(i, j, cal_list);
	cal_prev = cal_next;
	//printf("Swap between %i(%i) vs %i(%i),", i, cal_prev->read_area, j, cal_next->read_area);
      }
    }
  }

}
//-----------------------------------------------------------------------------
/*
size_t bwt_generate_cal_rna_list_linked_list(array_list_t *mapping_list,
					     cal_optarg_t *cal_optarg,
					     size_t *min_seeds, 
					     size_t *max_seeds,
					     size_t nchromosomes,
					     array_list_t *cal_list,
					     size_t read_length) {
  //printf("::: CALS PROCESS with max seeds %i:::\n", *max_seeds);

  short_cal_t *short_cal;
  linked_list_item_t *list_item_cal, list_item_seed;
  region_t *region;
  size_t min_cal_size = cal_optarg->min_cal_size;
  size_t max_cal_distance  = cal_optarg->max_cal_distance;
  size_t num_mappings = array_list_size(mapping_list);
  size_t chromosome_id;
  short int strand;
  size_t start, end;
  linked_list_iterator_t itr;
  const unsigned char nstrands = 2;
  
  linked_list_t ***cals_list = (linked_list_t ***)malloc(sizeof(linked_list_t **)*nstrands);

  for (unsigned int i = 0; i < nstrands; i++) {
    cals_list[i] = (linked_list_t **)malloc(sizeof(linked_list_t *)*nchromosomes);
    for (unsigned int j = 0; j < nchromosomes; j++) {
      cals_list[i][j] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    }
  }

  // from the results mappings, generates CALs
  for (unsigned int m = 0; m < num_mappings; m++) {
    region = array_list_get(m, mapping_list);
    chromosome_id = region->chromosome_id;
    strand = region->strand;
    my_cp_list_append_linked_list(cals_list[strand][chromosome_id], region, max_cal_distance, *max_seeds);     
  }


  //  printf("End inserts. Select CALs\n");
  //Store CALs in Array List for return results
  size_t start_cal, end_cal, len;
  size_t seq_start, seq_end;
  seed_region_t *s;
  array_list_t *short_cals_list = array_list_new(256,
                                                 1.25f,
                                                 COLLECTION_MODE_ASYNCHRONIZED);
  short_cal_t *short_cal_prev = NULL, *short_cal_last = NULL, *short_cal_first = NULL;
  size_t max_intron_size = 1000000;
  int pending_insert;

  //seed_region_t *s;
  for (unsigned int j = 0; j < nchromosomes; j++) {
    for (unsigned int i = 0; i < nstrands; i++) {
      linked_list_iterator_init(cals_list[i][j], &itr);
      list_item_cal = linked_list_iterator_list_item_curr(&itr);
      short_cal_first = NULL;
      pending_insert = 0;
      while (list_item_cal != NULL ) {
	short_cal = (short_cal_t *)list_item_cal->item;
	if (short_cal->end - short_cal->start + 1 >= min_cal_size) {
	  //	  printf("short_cal = [%i:%lu - %lu]\n", j, short_cal->start, short_cal->end);
	  if (!short_cal_first) {
	    //	    printf("************* SELECT FIRST ***************\n");
	    short_cal_first = short_cal;
	    pending_insert = 1;
	  } else {
	    if (short_cal->start <= (short_cal_first->end + max_intron_size)) {
	      short_cal_first->end = short_cal->end;
	      short_cal_first->num_seeds += short_cal->num_seeds;
	      while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_list)) {
		seed_region_select_linked_list(short_cal_first->sr_list, short_cal_first->sr_duplicate_list, 
					       s->read_start, s->read_end, s->genome_start, s->genome_end, s->id,
					       short_cal_first->seeds_ids_array);
	      }
	      linked_list_free(short_cal->sr_list, NULL);
	      while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_duplicate_list)) {
		linked_list_insert_first(s, short_cal_first->sr_duplicate_list);
	      }
	      linked_list_free(short_cal->sr_duplicate_list, NULL);
	      short_cal_free(short_cal);
	    } else {
	      //Insertar en la lista de CALs short_cal_first
	      linked_list_t *list_aux = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	      while (s = (seed_region_t *)linked_list_remove_last(short_cal_first->sr_list)) {
		//TODO: Change all parameters to seed_region_t
		append_seed_region_linked_list(list_aux,
					       s->read_start, s->read_end,
					       s->genome_start, s->genome_end, 
					       s->id);	    
		seed_region_free(s);
	      }
	      array_list_insert(cal_new(j, i, short_cal_first->start, 
					short_cal_first->end, short_cal_first->num_seeds, 
					list_aux, short_cal_first->sr_duplicate_list), cal_list);
	      short_cal_first = short_cal;
	    }	    
	  }
	  //TODO: Free short cal
	  array_list_insert(short_cal, short_cals_list);
	} else {
	  //Not min CAL size free short_cal
	  linked_list_free(short_cal->sr_list, free);
	  linked_list_free(short_cal->sr_duplicate_list, free);
	  short_cal_free(short_cal);
	}
	
	linked_list_iterator_next(&itr);
	//linked_list_item_free(list_item_cal, short_cal_free);

	list_item_cal = linked_list_iterator_list_item_curr(&itr);

      }
      
      if (pending_insert) {
	//Insertar en la lista de CALs short_cal_first
	linked_list_t *list_aux = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	while (s = (seed_region_t *)linked_list_remove_last(short_cal_first->sr_list)) {
	  //TODO: Change all parameters to seed_region_t
	  append_seed_region_linked_list(list_aux,
					 s->read_start, s->read_end, 
					 s->genome_start, s->genome_end, 
					 s->id);	    
	  seed_region_free(s);		
	}
	array_list_insert(cal_new(j, i, short_cal_first->start, 
				  short_cal_first->end, short_cal_first->num_seeds, 
				  list_aux, short_cal_first->sr_duplicate_list), cal_list);
      }

      free(cals_list[i][j]);
      //linked_list_free(cals_list[i][j], (void *)short_cal_free);
    }
  }

  for (unsigned int i = 0; i < nstrands; i++) {
    free(cals_list[i]);
  }

  free(cals_list);
  
  order_cal_linked_list(array_list_size(cal_list), cal_list);

  /*printf(":: CALS RESULT: \n");

  /*
  for (int i = 0; i < array_list_size(cal_list); i++) {
    cal_t *cal = array_list_get(i, cal_list);
    printf("\tCAL%i:= Num Seeds: %i, chr %i:(%i)[%lu-%lu], long: %i\n",i, cal->num_seeds,
	   cal->chromosome_id, cal->strand, cal->start, cal->end, cal->read_area);

    printf("\tTotal Seeds Regions Uniq %lu: \n", linked_list_size(cal->sr_list));
    for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
      seed_region_t *s = list_item->item;
      printf("\t\t[%i|%i - %i|%i]\n", s->genome_start, s->read_start, s->read_end, s->genome_end);
    }
    printf("\n");

    printf("\tTotal Seeds Regions Duplicate %lu: \n", linked_list_size(cal->sr_duplicate_list));
    for (linked_list_item_t *list_item = cal->sr_duplicate_list->first; list_item != NULL; list_item = list_item->next) {
      seed_region_t *s = list_item->item;
      printf("\t\t[%i|%i - %i|%i]\n", s->genome_start, s->read_start, s->read_end, s->genome_end);
    }
    printf("\n");
  }


  return array_list_size(cal_list);
}
*/

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
	      printf("mapping: start = %lu end = %lu\n", start, end);

	    if (start < cal->start) {
	      if (print) 
		printf("\tChange start: chrm %lu CAL %d [%lu-%lu] -> new start = %lu\n", cal->chromosome_id, c, cal->start, cal->end, start);
	      cal->start = start; 
	    }
	    if (end > cal->end) {
	      if (print) 
		printf("\tChange end: chrm %lu CAL %d [%lu-%lu] -> new end = %lu\n", cal->chromosome_id, c, cal->start, cal->end, end);
	      cal->end = end; 
	    }
	    if (print) 
	      printf("\t\tExtend CAL %d Result :: %lu-%lu\n", c, cal->start, cal->end);
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
		  printf("\tFound CAL in range. Start End Actualization cal:[%lu-%lu] extend cal1:[%lu-%lu]\n", cal->start, cal->end, cal1->start, cal1->end);
		if(start > cal1->start){
		  cal->start = cal1->start;
		  extend = 1;
		}
		if(end < cal1->end){
		  cal->end = cal1->end;
		  extend = 1;
		}
		if (print) 
		  printf("\t\tResult cal [%lu-%lu] and Delete %d\n", cal->start, cal->end, extend);
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
      array_list_insert(cal_new(chromosome_id, strand, start, end, 0, NULL, NULL),
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
  unsigned int num_mappings = bwt_map_exact_seeds_seq(5, 5, seq, cal_optarg->seed_size,
						      cal_optarg->min_seed_size,
						      bwt_optarg, index, mapping_list, 0);
  

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
/*
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
*/
//-----------------------------------------------------------------------------
/*
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
*/
//-----------------------------------------------------------------------------
size_t bwt_generate_cal_list_linkedlist(array_list_t *mapping_list,
					cal_optarg_t *cal_optarg,
					int *min_seeds, int *max_seeds,
					size_t nchromosomes,
					array_list_t *cal_list) {
  /*
  short_cal_t *short_cal_p;
  region_t *region;
  size_t min_cal_size = cal_optarg->min_cal_size;
  size_t max_cal_distance  = cal_optarg->max_cal_distance;
  size_t num_mappings = array_list_size(mapping_list);
  size_t chromosome_id;
  short int strand;
  size_t start, end;  
  linked_list_item_t *list_item;
  linked_list_iterator_t itr;

  *max_seeds = 0;
  *min_seeds = 1000;

  const unsigned char nstrands = 2;
  linked_list_t ***cals_list = (linked_list_t ***)malloc(sizeof(linked_list_t **)*nstrands);

  for (unsigned int i = 0; i < nstrands; i++) {
    cals_list[i] = (linked_list_t **)malloc(sizeof(linked_list_t *)*nchromosomes);
    for (unsigned int j = 0; j < nchromosomes; j++) {
      cals_list[i][j] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    }
  }
  
  // from the results mappings, generates CALs
  for (unsigned int m = 0; m < num_mappings; m++) {
    region = array_list_get(m, mapping_list);

    chromosome_id = region->chromosome_id;
    strand = region->strand;

    my_cp_list_append_linked_list(cals_list[strand][chromosome_id], region, max_cal_distance);
  }

  // for debugging
  //  printf("num. seed mappings: %i\n", array_list_size(mapping_list));

  //Store CALs in Array List for return results
  cal_t *cal;
  for (unsigned int j = 0; j < nchromosomes; j++) {
    for (unsigned int i = 0; i < nstrands; i++) {
      linked_list_iterator_init(cals_list[i][j], &itr);
      list_item = linked_list_iterator_list_item_curr(&itr);

      // for debugging
      //      if (list_item != NULL) {
      //	printf("\t: strand %c\t chrom. %i: num. cals = %i\n", (i == 0 ? '-' : '+'), j, linked_list_size(cals_list[i][j]));
      //      }
      while (list_item != NULL) {
	short_cal_p = (short_cal_t *)list_item->item;
	if (short_cal_p->end - short_cal_p->start + 1 >= min_cal_size) {

	  if (*min_seeds > short_cal_p->num_seeds) *min_seeds = short_cal_p->num_seeds;
	  if (*max_seeds< short_cal_p->num_seeds) *max_seeds = short_cal_p->num_seeds;

	  if (i) {
	    // strand -
	    cal = cal_new(j, i, 
			  short_cal_p->start - (short_cal_p->seq_len - short_cal_p->seq_start) + 1, 
			  short_cal_p->end + short_cal_p->seq_end, short_cal_p->num_seeds, 0, 0);
	  } else {
	    // strand +
	    cal = cal_new(j, i, 
			  short_cal_p->start - short_cal_p->seq_start, 
			  short_cal_p->end + (short_cal_p->seq_len - short_cal_p->seq_end) + 1,
			  short_cal_p->num_seeds, 0, 0);
	  }

	  array_list_insert(cal, cal_list);
	}

        linked_list_iterator_next(&itr);
	linked_list_item_free(list_item, short_cal_free);
        list_item = linked_list_iterator_list_item_curr(&itr);
      }
      free(cals_list[i][j]);
    }
  }
  
  for (unsigned int i = 0; i < nstrands; i++) {
    free(cals_list[i]);
  }

  free(cals_list);

  // for debugging
  //  printf("--->final num. cals = %i\n", array_list_size(cal_list));

  return array_list_size(cal_list);  */
}
