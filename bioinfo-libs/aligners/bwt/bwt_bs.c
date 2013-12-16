#include "bwt.h"

//-----------------------------------------------------------------------------
// added by JT and PP
size_t __bwt_map_inexact_read_bs(fastq_read_t *read,
				 bwt_optarg_t *bwt_optarg, 
				 bwt_index_t *index, 
				 array_list_t *mapping_list, 
				 int type);

void *__bwt_generate_anchor_list(size_t k_start, size_t l_start, int len_calc, 
				 bwt_optarg_t *bwt_optarg, 
				 bwt_index_t *index, int type, 
				 array_list_t *anchor_list, int type_anchor,
				 int dsp);

// end of added by JT and PP
//-----------------------------------------------------------------------------

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

size_t bwt_map_exact_seed(char *seq, size_t seq_len, 
			  size_t start, size_t end,
			  bwt_optarg_t *bwt_optarg,
			  bwt_index_t *index,
			  array_list_t *mapping_list);

//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seq(char *seq, 
			   bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, 
			   array_list_t *mapping_list);

size_t bwt_map__forward_inexact_seq(char *seq, 
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

//------------------------------------------------------------------------------

size_t bwt_map_inexact_seed(char *seq, size_t seq_len, 
			    size_t start, size_t end,
			    bwt_optarg_t *bwt_optarg,
			    bwt_index_t *index,
			    array_list_t *mapping_list);

//------------------------------------------------------------------------------

char *bwt_error_type(char error_kind);

//-----------------------------------------------------------------------------
/*
void bwt_generate_index_files_bs(char *ref_file, char *output_dir, 
				 unsigned int s_ratio, char *bases) {

  byte_vector X, B, Bi;
  vector C, C1;
  comp_vector S, Si, Scomp, Scompi;
  comp_vector R, Ri, Rcomp, Rcompi;
  comp_matrix O, Oi;

  exome ex;

  //initReplaceTable();
  initReplaceTable_bs(bases);
  saveNucleotide(bases, output_dir, "Nucleotide");

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
*/
//-----------------------------------------------------------------------------
// general bs functions
//-----------------------------------------------------------------------------

size_t bwt_map_exact_seed_bs(char *seq, size_t seq_len,
			     size_t seq_start, size_t seq_end,
			     bwt_optarg_t *bwt_optarg, 
			     bwt_index_t *index, 
			     array_list_t *mapping_list,
			     int id, int type) {
  
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
  //  array_list_t *mappings = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, 
  //						     COLLECTION_MODE_ASYNCHRONIZED);
  int discard_seed = 0;
  int actual_mappings = 0;
  struct timeval t_start, t_end;

  result.k = 0;
  result.l = size_SA(index->backward) - 1;
  result.start = start;
  result.end = end;
  if (type == 1) {
    // strand +
    result.pos = end;
    aux_seq_start = seq_start;
    aux_seq_end = seq_end;        
    BWExactSearchBackward(code_seq, index->backward, &result, NULL);
  } else {
    // strand -
    aux_seq_start = seq_len - seq_end - 1;
    aux_seq_end = seq_len - seq_start - 1;
    result.pos = start; 
    BWExactSearchForward(code_seq, index->backward_rev, &result, NULL);
  }
  
  //start_timer(t_start);
  
  k_aux = result.k;
  l_aux = result.l;
  actual_mappings += (result.l - result.k + 1);
  
  //if (actual_mappings > 150) {//bwt_optarg->filter_seed_mappings) {
  if (actual_mappings > bwt_optarg->filter_seed_mappings) {
    //discard_seed = 1;
    //break;
    //printf("**************** EXCEED mappings: %lu\n", actual_mappings);
    k_aux = result.k;
    l_aux = result.k + 10;
  } else {
    //printf("\tk=%d - l=%d\n", r->k, r->l);      
    k_aux = result.k;
    l_aux = result.l;
  }
      
  for (size_t j = k_aux; j <= l_aux; j++) {
      key = get_SA(j, index->backward);      
      idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
      
      if (key + len <= index->karyotype.offset[idx]) {
	//start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]); 
	start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	printf("\t\tstrand:%c\tchromosome:%s\tstart:%u\tend:%u\n",plusminus[type],
	       index->karyotype.chromosome + (idx-1) * IDMAX,
	       start_mapping, start_mapping + len);

	// save all into one alignment structure and insert to the list
	region = region_bwt_new(idx, !type, start_mapping, start_mapping + len, aux_seq_start, aux_seq_end, seq_len, id);

	assert(region != NULL);

	if (!array_list_insert((void*) region, mapping_list)){
	  printf("Error to insert item into array list\n");
	}

	num_mappings++;

      }
    }

  return num_mappings;  
}

//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------

size_t bwt_map_inexact_read_bs(fastq_read_t *read, 
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       array_list_t *mapping_list,
			       int type) {
  
  return __bwt_map_inexact_read_bs(read, 
				   bwt_optarg, 
				   index, 
				   mapping_list,
				   type);  
}

size_t __bwt_map_inexact_read_bs(fastq_read_t *read,
				 bwt_optarg_t *bwt_optarg, 
				 bwt_index_t *index, 
				 array_list_t *mapping_list,
				 int type) {
  
  alignment_t *alignment;
  char *seq = read->sequence; 
  size_t len = read->length;
     
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

  char *code_seq = (char *) malloc(len * sizeof(char));

  encode_bases(code_seq, seq, len, index->bwt_config.table);

  // calculate vectors k and l
  intmax_t *k1 = (intmax_t *) malloc(len * sizeof(intmax_t));
  intmax_t *l1 = (intmax_t *) malloc(len * sizeof(intmax_t));
  intmax_t *ki1 = (intmax_t *) malloc(len * sizeof(intmax_t));
  intmax_t *li1 = (intmax_t *) malloc(len * sizeof(intmax_t));

  size_t last_k1, last_l1;
  size_t last_ki1, last_li1;

  int back1_nt, forw1_nt;
  back1_nt = BWExactSearchVectorBackward(code_seq, start, end, 0, size_SA(index->backward)-1,
					 k1, l1, index->backward);
  
  forw1_nt = BWExactSearchVectorForward(code_seq, start, end, 0, size_SA(index->forward)-1,
					ki1, li1, index->forward);

  // compare the vectors k and l to get mappings in the genome
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
  //  size_t i, j, z;
  size_t *allocate_pos_alignments;
  size_t k_start, l_start;

  size_t start_mapping;
  size_t tot_alignments = 0;
  const int MAX_BWT_ALIGNMENTS = 10;
  int filter_exceeded = 0;
  //seq_dup = (char *)malloc(sizeof(char)*(len + 1));
  seq_strand = strdup(seq);
  error = MISMATCH;


  char *quality_clipping = (char *) malloc(sizeof(char) * 50);
  seq_dup = (char *) malloc(sizeof(char) * (len + 1));

  new_results_list(&r_list, bwt_optarg->filter_read_mappings);

  array_list_t *tmp_mapping_list = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, 
						  COLLECTION_MODE_ASYNCHRONIZED);

  //for (int type = 1; type >= 0; type--) {
  r_list.num_results = 0;
  r_list.read_index = 0;

  //    printf("*** bwt.c: calling BWSearch1 with type = %d...\n", type);
  //if (type == 1) {
  BWSearch1VectorHelper(code_seq, start, end, k1, l1, ki1, li1, 
			index->backward, index->forward, 
			&r_list, index->bwt_config.nA);

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
	  

    //if (type) {
    direction = r->dir;
    //} else {
    //direction = !r->dir;
    //}

    len_calc = len;
    if (error == DELETION) {
      len_calc--;
    } else if (error == INSERTION) {
      len_calc++;
    }

    // generating cigar
    sprintf(quality_clipping, "%i", NONE_HARD_CLIPPING);
    if (error == 0) {
      sprintf(cigar, "%lu=\0", len);
      num_cigar_ops = 1;
      memcpy(seq_dup, seq_strand, len);
      seq_dup[len] = '\0';
		 
    } else if (error == MISMATCH) {
      if (pos == 0) {
	//Positive strand
	//if(type) { 
	sprintf(cigar, "1S%luM\0", len-1); 
	start_mapping++;
	//}
	//else { 
	//sprintf(cigar, "%luM1S\0", len-1); 
	//}
	num_cigar_ops = 2;
      } else if (pos == len - 1) {
	//Positive strand
	//if(type) { 
	sprintf(cigar, "%luM1S\0", len - 1); 
	//}
	//else{ 
	//sprintf(cigar, "1S%luM\0", len-1); 
	//start_mapping++;
	//}
	num_cigar_ops = 2;
      } else {
	sprintf(cigar, "%luM\0", len);
	num_cigar_ops = 1;
      }
      memcpy(seq_dup, seq_strand, len);
      seq_dup[len] = '\0';
      //printf("MISMATCH\n");
		 
    } else if (error == INSERTION) {
      //printf("INSERTION\n");
      if (pos == 0) {
	//if(type) {
	sprintf(cigar, "1M1D%luM\0", len - 1); 
	//}
	//else{ 
	//sprintf(cigar, "%luM1D1M\0", len - 1); 
	//}	      
      } else if (pos == len - 1) {
	//if(type) { 
	sprintf(cigar, "%luM1D1M\0", len - 1); 
	//}
	//else{ 
	//sprintf(cigar, "1M1D%luM\0", len - 1); 
	//}
      } else {
		   
	//if(type) {
	if(r->dir)
	  sprintf(cigar, "%iM1D%luM\0", pos, len - pos);
	else
	  sprintf(cigar, "%iM1D%luM\0", pos + 1, len - pos - 1);
	/*} else { 
	  if(r->dir)
	  sprintf(cigar, "%luM1D%dM\0", len - pos, pos);
	  else
	  sprintf(cigar, "%luM1D%dM\0", len - pos - 1, pos + 1);
	  }*/
      }
      num_cigar_ops = 3;
      memcpy(seq_dup, seq_strand, len);
      seq_dup[len] = '\0';
    } else if (error == DELETION) {	     
      //printf("DELETION\n");
      if (pos == 0) {
	//if(type) { 
	sprintf(cigar, "1I%luM\0", len -1); 
	//}
	//else{ 
	//sprintf(cigar, "%luM1I\0", len -1); 
	//		   start_mapping++;
	//}
		   
	num_cigar_ops = 2;		
      } else if (pos == len - 1) {
	//if(type) { 
	sprintf(cigar, "%luM1I\0", len -1); 
	//		   start_mapping++;
	//}
	// else{ 
	//sprintf(cigar, "1I%luM\0", len -1); 
	//}
	num_cigar_ops = 2;
      } else {
	//if(type) { 
	sprintf(cigar, "%dM1I%luM\0", pos, len - pos - 1); 
	//}
	//else{ 
	//sprintf(cigar, "%luM1I%dM\0", len - pos - 1, pos); 
	//}
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
    //printf("IN FUNCTION SEQ_DUP %d :: %s\n", strlen(seq_dup), seq_dup);


    //	 printf("\tstart_mapping = %lu, cigar = %s, type = %i, j (k to l) = %lu\n", start_mapping, cigar, type, j);


    //printf("Max alignments per read %i >= %i\n", r->l - r->k + 1, bwt_optarg->max_alignments_per_read);
    //if ( r->l - r->k + 1 > bwt_optarg->filter_read_mappings) {
    //k_start = r->k;
    //l_start = k_start +  bwt_optarg->filter_read_mappings;
    //filter_exceeded = 1;
    //type = -1;
    //break;
    //} else {
    k_start = r->k;
    l_start = r->l;
    //}


    tot_alignments += (l_start - k_start);
    // check filter_read_mappings for bisulfite case
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
      idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
      if(key + len_calc <= index->karyotype.offset[idx]) {
	start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	// save all into one alignment structure and insert to the list
	alignment = alignment_new();

	alignment_init_single_end(NULL, strdup(seq_dup), strdup(quality_clipping), type, 
				  idx - 1, //index->karyotype.chromosome + (idx-1) * IDMAX,
				  start_mapping, //index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]), 
				  strdup(cigar), num_cigar_ops, 254, 1, (num_mappings > 0), 0, NULL, alignment);
	array_list_insert((void*) alignment, tmp_mapping_list);
      }
    }//end for k and l
  }//end for 

  if (filter_exceeded) {
    //printf("Limit Exceeded %d\n", bwt_optarg->filter_read_mappings);
    array_list_clear(tmp_mapping_list, (void *)alignment_free);
    if (!mapping_list->size) {
      array_list_set_flag(2, mapping_list);       
    }
    goto exit;
    // break;
  }
     
  //
  //
  //Search for equal BWT mappings and set the mappings that will be delete
  int n_mappings = array_list_size(tmp_mapping_list);
  printf("**Finish BWT with n_mappings = %i\n", n_mappings);
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
     
  if (array_list_get_flag(mapping_list) == 1 &&
      n_mappings) {
    printf("\tNUM ITEMS TO FREE %i\n", array_list_size(mapping_list));
    array_list_clear(mapping_list, (void *)bwt_anchor_free);
    printf("\tFREE END\n");
  }

  //Delete all set mappings
  int primary_delete = 0;
  int header_len;
  for (int m = n_mappings - 1; m >= 0; m--) {
    alig_1 = array_list_remove_at(m, tmp_mapping_list);
    //alignment_print(alig_1);
    if (delete_mark[m]) {
      if (!is_secondary_alignment(alig_1)) { primary_delete = 1; }
      alignment_free(alig_1);
    } else {
      set_secondary_alignment(num_mappings > 0, alig_1);
      num_mappings++;
      //printf("header : %s\n", read->id);
      header_len = strlen(read->id);
      alig_1->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
      get_to_first_blank(read->id, header_len, alig_1->query_name);
      //printf("header alig: %s\n", alig_1->query_name);
      bwt_cigar_cpy(alig_1, read->quality);
      //************************* OPTIONAL FIELDS ***************************//
      alig_1 = add_optional_fields(alig_1, n_mappings, len);
      array_list_insert((void*) alig_1, mapping_list);
    }
  }

  if (primary_delete) { 
    alig_1 = array_list_get(0, mapping_list); 
    set_secondary_alignment(0, alig_1);
  }
  free(delete_mark);

  if (n_mappings == 0) {     
    if (array_list_get_flag(mapping_list) == 1) {
      //====================================================================================
      //printf("BWT(+): FORWARD(k-l)%i:%lu-%lu BACKWARD(k-l)%i:%lu-%lu\n", forw1_nt, last_ki1, last_li1, back1_nt, last_k1, last_l1);
      //printf("BWT(-): FORWARD(k-l)%i:%lu-%lu BACKWARD(k-l)%i:%lu-%lu\n", forw0_nt, last_ki0, last_li0, back0_nt, last_k0, last_l0);
      array_list_t *forward_anchor_list, *backward_anchor_list;       
      int new_type = !type; //1 == (+)

      __bwt_generate_anchor_list(k1[end - back1_nt + 1], 
				 l1[end - back1_nt + 1], back1_nt, bwt_optarg, 
				 index, type, mapping_list, BACKWARD_ANCHOR, 0);

      //printf("FORWARD (+)\n");
      __bwt_generate_anchor_list(ki1[start + forw1_nt - 1], 
				 li1[start + forw1_nt - 1], forw1_nt, bwt_optarg, 
				 index, type,  mapping_list, FORWARD_ANCHOR, 0);
      
      //printf("**BWT Anchors items = %i\n", mapping_list->size);

      //type = 0;//(-)
      //printf("BACKWARD (-)\n");
      /* __bwt_generate_anchor_list(last_k0, last_l0, back0_nt, bwt_optarg, 
	 index, type,  mapping_list, BACKWARD_ANCHOR, back0_nt - forw0_nt);
	 //printf("FORWARD (-)\n");
	 __bwt_generate_anchor_list(last_ki0, last_li0, forw0_nt, bwt_optarg, 
	 index, type,  mapping_list, FORWARD_ANCHOR, back0_nt - forw0_nt);*/
      /*
	__bwt_generate_anchor_list(last_k0, last_l0, forw0_nt, bwt_optarg, 
	index, type,  mapping_list, BACKWARD_ANCHOR, back0_nt - forw0_nt);
	//printf("FORWARD (-)\n");
	__bwt_generate_anchor_list(last_ki0, last_li0, back0_nt, bwt_optarg, 
	index, type,  mapping_list, FORWARD_ANCHOR, back0_nt - forw0_nt);
      */
    }
  } else {
    array_list_set_flag(0, mapping_list);
    //printf("########## EXACT! #########\n");
  }

 exit:
  free(r_list.list);
  free(code_seq);
  free(seq_strand);
  //free(k0);
  //free(l0);
  free(k1);
  free(l1);
  //free(ki0);
  //free(li0);
  free(ki1);
  free(li1);

  free(seq_dup);
  free(quality_clipping);
  array_list_free(tmp_mapping_list, NULL);

  //printf("\tOut function\n");

  return array_list_size(mapping_list);
}

//-----------------------------------------------------------------------------

void bwt_map_inexact_array_list_by_filter_bs(array_list_t *reads,
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

    //printf("read %lu (of %lu)\t", i, num_reads);
    //printf("%s\n", fq_read->sequence);
    num_mappings = bwt_map_inexact_read_bs(fq_read, 
					   bwt_optarg, index, 
					   lists[i], 0);

    //printf("end read %lu\n\n", i);    
    if (num_mappings > 0) {
      array_list_set_flag(1, lists[i]);
      for (size_t j = 0; j < num_mappings; j++) {
	alignment = (alignment_t *) array_list_get(j, lists[i]);
	header_len = strlen(fq_read->id);
	alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	
	get_to_first_blank(fq_read->id, header_len, alignment->query_name);
	bwt_cigar_cpy(alignment, fq_read->quality);
	
	//************************* OPTIONAL FIELDS ***************************//
	alignment = add_optional_fields(alignment, num_mappings, fq_read->length);
	//*********************** OPTIONAL FIELDS END ************************//
      }
    } else if (array_list_get_flag(lists[i]) != 2) {
	unmapped_indices[(*num_unmapped)++] = i;
	array_list_set_flag(0, lists[i]);
    }
  } 
}

//-----------------------------------------------------------------------------

inline size_t seeding_bs(char *code_seq, size_t seq_len, size_t num_seeds,
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

  start = 0;
  for (size_t i = 0; i < n_seeds; i++) {
    end = start + seed_size;
    if (end >= seq_len) end = seq_len;
    num_mappings = bwt_map_exact_seed_bs(code_seq, seq_len, start, end - 1,
					 bwt_optarg, index, mapping_list, seed_id, 0);
    seed_id++;
    total_mappings += num_mappings;
    //    LOG_DEBUG_F("\tseed %i\t[%i - %i], length read = %i, num. mappings = %i\n", 
    //		i + 1, start, end, seq_len, num_mappings);
    start += offset_inc;
    if (start > offset_end) {
      if (offset_inc == 1) break;
      start = offset_inc / 2;
    }
  }

  //  LOG_DEBUG_F("\t\ttotal mappings = %i\n", total_mappings);

  return total_mappings;
}

//-----------------------------------------------------------------------------

size_t bwt_map_exact_seeds_seq_by_num_bs(char *seq, size_t num_seeds, 
					 size_t seed_size, size_t min_seed_size,
					 bwt_optarg_t *bwt_optarg, bwt_index_t *index,
					 array_list_t *mapping_list) {
  size_t seq_len = strlen(seq);
  size_t num_mappings = 0;

  char *code_seq = (char *) calloc(seq_len, sizeof(char));

  //replaceBases(seq, code_seq, seq_len);
  //encodeBases(seq, code_seq, seq_len);
  encode_bases(code_seq, seq, seq_len, index->bwt_config.table);
  

  num_mappings = seeding_bs(code_seq, seq_len, num_seeds, seed_size, min_seed_size,
  			    bwt_optarg, index, mapping_list);
  //  printf("\tfirst, num_mappings = %d\n", num_mappings);

  free(code_seq);
  return num_mappings;
}

//-----------------------------------------------------------------------------

size_t bwt_generate_cals_bs(char *seq, char *seq2, size_t seed_size, bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, bwt_index_t *index2, array_list_t *cal_list) {

  //printf("New function generate CALs\n");
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

  int seed_id = 0;
  size_t len = strlen(seq);
  size_t len2 = strlen(seq2);

  size_t offset, num_seeds;
  
  //printf(" len=%i, seed_size=%i, num_seeds=%i, seq=%s\n", len, seed_size, num_seeds, seq);
  char *code_seq = (char *) calloc(len, sizeof(char));
  char *code_seq2 = (char *) calloc(len, sizeof(char));

  // be careful, now we use the table from the index
  // replaceBases(seq, code_seq, len);
  encode_bases(code_seq, seq, len, index->bwt_config.table);
  encode_bases(code_seq2, seq2, len2, index->bwt_config.table);

  // second first seed
  if (seed_size <= 10) { seed_size = 16; }

  int padding_left = seed_size / 2;
  int padding_right = padding_left;

  int extra_seed_size = 16;//seed_size;
  //int first_seed = 16;

  /*
  // extra seed for splice junctions
  bwt_map_exact_seed_bs(code_seq, len, padding_left, padding_left + extra_seed_size - 1,
			bwt_optarg, index, mapping_list, seed_id++);

  insert_seeds_and_merge(mapping_list, cals_list,  len);
  */

  {
    size_t num_mappings;
    size_t num_seeds = 10;
    size_t start, end, offset_end = len - 16;
    size_t n_seeds = num_seeds;
    size_t offset_inc = ceil(1.0f * len / (num_seeds + 1));
    if (offset_inc <= 0) offset_inc = 1;

    printf("+++++++++++ read1 %s ++++++++++++++\n", seq);
    printf("+++++++++++ read2 %s ++++++++++++++\n", seq2);

    start = 0;
    for (size_t i = 0; i < n_seeds; i++) {
      end = start + seed_size;
      if (end >= len) end = len;

      //start = 10;
      //end = 33;
      
      num_mappings = bwt_map_exact_seed_bs(code_seq, len, start, end - 1,
					   bwt_optarg, index,  mapping_list, seed_id++, 0);
      printf("----------- seed %lu (%lu - %lu) -> %lu mappings\n", seed_id, start, end, num_mappings);
      //transform_regions(mapping_list);
      insert_seeds_and_merge(mapping_list, cals_list,  len);
      
      num_mappings = bwt_map_exact_seed_bs(code_seq2, len2, start, end - 1,
					   bwt_optarg, index2, mapping_list, seed_id++, 1);
      printf("----------- seed %lu (%lu - %lu) -> %lu mappings\n", seed_id, start, end, num_mappings);
      insert_seeds_and_merge(mapping_list, cals_list,  len);
      
      start += offset_inc;
      if (start > offset_end) {
	if (offset_inc == 1) break;
	start = offset_inc / 2;
      }
    }
  }




  /*
  num_seeds = len / seed_size;

  // first 'pasada'
  offset = 0;
  for (size_t i = 0; i < num_seeds; i++) {
    bwt_map_exact_seed_bs(code_seq, len, offset, offset + seed_size - 1,
			  bwt_optarg, index, mapping_list, seed_id++);
    transform_regions(mapping_list);
    insert_seeds_and_merge(mapping_list, cals_list,  len);

    bwt_map_exact_seed_bs(code_seq2, len2, offset, offset + seed_size - 1,
			  bwt_optarg, index2, mapping_list, seed_id++);
    insert_seeds_and_merge(mapping_list, cals_list,  len);

    offset += seed_size;
  }
  
  // last seed
  if (len % seed_size > 0) {
    bwt_map_exact_seed_bs(code_seq, len, len - seed_size, len - 1,
			  bwt_optarg, index, mapping_list, seed_id++);
    transform_regions(mapping_list);
    insert_seeds_and_merge(mapping_list, cals_list,  len);

    bwt_map_exact_seed_bs(code_seq2, len2, len2 - seed_size, len2 - 1,
			  bwt_optarg, index2, mapping_list, seed_id++);
    insert_seeds_and_merge(mapping_list, cals_list,  len);
  }
  */
  /*
  if (len % seed_size != padding_right) {
    // extra seed for splice junctions
    bwt_map_exact_seed(code_seq, len, len - extra_seed_size - padding_right, len - padding_right - 1,
		       bwt_optarg, index, mapping_list, seed_id++);
    insert_seeds_and_merge(mapping_list, cals_list,  len);
  }
  */
  free(code_seq);
  free(code_seq2);

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

