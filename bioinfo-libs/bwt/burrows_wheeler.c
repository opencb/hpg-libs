#include "burrows_wheeler.h"

/*
//-----------------------------------------------------------------------------
// CAL policy
//-----------------------------------------------------------------------------

cal_optarg_t *cal_optarg_new(const unsigned int num_seeds, const unsigned int max_seed_distance,
			     const unsigned int left_flank_length, const unsigned in right_flank_length);

void cal_optarg_free(cal_optarg_t *optarg);
*/
//-----------------------------------------------------------------------------
// Burrows-Wheeler Transform
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Index for the Burrows-Wheeler transform
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

unsigned int bwt_map_seq_cpu(char *seq, 
			     bwt_optarg_t *bwt_optarg, 
			     bwt_index_t *index, 
			     array_list_t *mapping_list) {

  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_seq_cpu(seq, bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_seq_cpu(seq, bwt_optarg, index, mapping_list);  
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_read_cpu(fastq_read_t *read, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list) {

  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_read_cpu(read, bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_read_cpu(read, bwt_optarg, index, mapping_list);  
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_seqs_cpu(char **seqs, 
			      unsigned int num_reads,
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list) {

  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_seqs_cpu(seqs, num_reads, 
				  bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_seqs_cpu(seqs, num_reads, 
				  bwt_optarg, index, mapping_list);  
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_batch_cpu(fastq_batch_t *batch,
			       bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, 
			       fastq_batch_t *unmapped_batch,
			       array_list_t *mapping_list) {

  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_batch_cpu(batch, bwt_optarg, 
				   index, unmapped_batch, mapping_list);
  }

  return bwt_map_inexact_batch_cpu(batch, bwt_optarg, 
				   index, unmapped_batch, mapping_list);  
}

//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_seq_cpu(char *seq, 
				   bwt_optarg_t *bwt_optarg, 
				   bwt_index_t *index, 
				   array_list_t *mapping_list) {
  return 0;
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_read_cpu(fastq_read_t *read, 
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index, 
				    array_list_t *mapping_list) {
  return 0;
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_seqs_cpu(char **seqs, 
				    unsigned int num_reads,
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index, 
				    array_list_t *mapping_list) {
  return 0;
}


//-----------------------------------------------------------------------------

unsigned int bwt_map_exact_batch_cpu(fastq_batch_t *batch,
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index, 
				     fastq_batch_t *unmapped_batch,
				     array_list_t *mapping_list) {
  return 0;
}

//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seq_cpu(char *seq, 
				     bwt_optarg_t *bwt_optarg, 
				     bwt_index_t *index, 
				     array_list_t *mapping_list) {

  unsigned int len = strlen(seq);
  int start = 0;
  int end = len - 1;

  char *codeSeq = (char *) calloc(len, sizeof(char));

  replaceBases(seq, codeSeq, len);

  // calculate vectors k and l
  size_t *k0, *l0, *k1, *l1;
  size_t *ki0, *li0, *ki1, *li1;
  
  BWExactSearchBackwardVector(codeSeq, start, end, 0, index->h_O.m_count - 2,
			      &k1, &l1, &index->h_C, &index->h_C1, &index->h_O);

  BWExactSearchForwardVector(codeSeq, start, end, 0, index->h_Oi.m_count - 2,
			     &ki1, &li1, &index->h_C, &index->h_C1, &index->h_Oi);

  BWExactSearchForwardVector(codeSeq, start, end, 0, index->h_rO.m_count - 2,
			     &k0, &l0, &index->h_rC, &index->h_rC1, &index->h_rO);

  BWExactSearchBackwardVector(codeSeq, start, end, 0, index->h_rOi.m_count - 2,
			      &ki0, &li0, &index->h_rC, &index->h_rC1, &index->h_rOi);

  // compare the vectors k and l to get mappings in the genome
  unsigned int num_mappings = 0;
  unsigned int found = 0, found2 = 0;
  char plusminus[2] = "-+";
  int idx, key, direction;
  results_list *r_list;
  result *r;

  for (int type = 1; type >= 0; type--) {

    r_list = (results_list *) calloc(1, sizeof(results_list));
    new_results_list(r_list, 2000);
    r_list->n = 0;
    r_list->read_index = 0;

    if (type == 1) {
      BWIterativeSearch1(codeSeq, start, end, k1, l1, ki1, li1, 
			 &index->h_C, &index->h_C1, &index->h_O, &index->h_Oi, r_list);
    } else {
      BWIterativeSearch1(codeSeq, start, end, ki0, li0, k0, l0, 
     			 &index->h_rC, &index->h_rC1, &index->h_rOi, &index->h_rO, r_list);
    }

    found2 = 0;

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
	  found2 = 1;
	  /*
	  alignment_p->flags = CODED_SEQ_FLAG;
	  alignment_p->header_p = &read_batch_p->header[read_batch_p->header_indices[i]];
	  alignment_p->read_p = seq;
	  alignment_p->quality_p = &read_batch_p->quality[read_batch_p->data_indices[i]];

	  mapping_info_p = mapping_info_new(genome_p->chromosome + (index-1) * IDMAX,
					    plusminus[type],
					    genome_p->start[index-1] + (key - genome_p->offset[index-1]));
	  memcpy(&mapping_info_p->errors, r, sizeof(result));
	  mapping_item_p = list_item_new(0, 0, mapping_info_p);
	  list_insert_item(mapping_item_p, &alignment_p->mapping_list);
	  */


	  printf("%s\t%c\t%s %u %s error: %i, pos: %i, base: %i\n",
		 "nothing", plusminus[type],
		 index->karyotype.chromosome + (idx-1) * IDMAX,
		 index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]),
		 seq, r->err_kind[0], r->position[0], r->base[0]);

	  num_mappings++;
	}
      }
    }
    found = found || found2;

    free(r_list);
  } // end for type 
  

  return num_mappings;
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_read_cpu(fastq_read_t *read, 
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list) {

  return bwt_map_inexact_seq_cpu(read->sequence, bwt_optarg, index, mapping_list);
}

//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_seqs_cpu(char **seqs, 
				      unsigned int num_reads,
				      bwt_optarg_t *bwt_optarg, 
				      bwt_index_t *index, 
				      array_list_t *mapping_list) {
  return 0;
}


//-----------------------------------------------------------------------------

unsigned int bwt_map_inexact_batch_cpu(fastq_batch_t *batch,
				       bwt_optarg_t *bwt_optarg, 
				       bwt_index_t *index, 
				       fastq_batch_t *unmapped_batch,
				       array_list_t *mapping_list) {
  return 0;
}

//-----------------------------------------------------------------------------
// cal functions
//-----------------------------------------------------------------------------

unsigned int bwt_find_cals_from_seq_cpu(char *seq, 
					bwt_optarg_t *bwt_optarg, 
					bwt_index_t *index, 
					cal_optarg_t *cal_optarg, 
					array_list_t *cal_list) {
  return 0;
}

//-----------------------------------------------------------------------------

unsigned int bwt_find_cals_from_read_cpu(fastq_read_t *read, 
					 bwt_optarg_t *bwt_optarg, 
					 bwt_index_t *index, 
					 cal_optarg_t *cal_optarg, 
					 array_list_t *cal_list) {
  return 0;
}

//-----------------------------------------------------------------------------

unsigned int bwt_find_cals_from_seqs_cpu(char **seqs, 
					 unsigned int num_reads,
					 bwt_optarg_t *bwt_optarg, 
					 bwt_index_t *index, 
					 unsigned int *num_unmapped,
					 int *unmapped_array,
					 cal_optarg_t *cal_optarg, 
					 array_list_t *cal_list) {
  return 0;
}

//-----------------------------------------------------------------------------

unsigned int bwt_find_cals_from_batch_cpu(fastq_batch_t *batch,
					  bwt_optarg_t *bwt_optarg, 
					  bwt_index_t *index, 
					  fastq_batch_t *unmapped_batch,
					  cal_optarg_t *cal_optarg, 
					  array_list_t *cal_list) {
  return 0;
}

//-----------------------------------------------------------------------------
/*
unsigned int bwt_single_read_cpu(fastq_read_t *read, bwt_optarg_t *optarg, 
				 bwt_index_t *index, 
				 void *mappings) {

  unsigned int len = strlen(read->sequence);
  int start = 0;
  int end = len - 1;

  char *seq = (char *) calloc(len, sizeof(char));

  replaceBases(read->sequence, seq, len);

  // calculate vectors k and l
  size_t *k0, *l0, *k1, *l1;
  size_t *ki0, *li0, *ki1, *li1;
  
  BWExactSearchBackwardVector(seq, start, end, 0, index->h_O.m_count - 2,
			      &k1, &l1, &index->h_C, &index->h_C1, &index->h_O);

  BWExactSearchForwardVector(seq, start, end, 0, index->h_Oi.m_count - 2,
			     &ki1, &li1, &index->h_C, &index->h_C1, &index->h_Oi);

  BWExactSearchForwardVector(seq, start, end, 0, index->h_rO.m_count - 2,
			     &k0, &l0, &index->h_rC, &index->h_rC1, &index->h_rO);

  BWExactSearchBackwardVector(seq, start, end, 0, index->h_rOi.m_count - 2,
			      &ki0, &li0, &index->h_rC, &index->h_rC1, &index->h_rOi);

  // compare the vectors k and l to get mappings in the genome
  unsigned int num_mappings = 0;
  unsigned int found = 0, found2 = 0;
  char plusminus[2] = "-+";
  int idx, key, direction;
  results_list *r_list;
  result *r;

  for (int type = 1; type >= 0; type--) {

    r_list = (results_list *) calloc(1, sizeof(results_list));
    new_results_list(r_list, 2000);
    r_list->n = 0;
    r_list->read_index = 0;

    if (type == 1) {
      BWIterativeSearch1(seq, start, end, k1, l1, ki1, li1, 
			 &index->h_C, &index->h_C1, &index->h_O, &index->h_Oi, r_list);
    } else {
      BWIterativeSearch1(seq, start, end, ki0, li0, k0, l0, 
     			 &index->h_rC, &index->h_rC1, &index->h_rOi, &index->h_rO, r_list);
    }

    found2 = 0;

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
	  found2 = 1;
	  /*
	  alignment_p->flags = CODED_SEQ_FLAG;
	  alignment_p->header_p = &read_batch_p->header[read_batch_p->header_indices[i]];
	  alignment_p->read_p = seq;
	  alignment_p->quality_p = &read_batch_p->quality[read_batch_p->data_indices[i]];

	  mapping_info_p = mapping_info_new(genome_p->chromosome + (index-1) * IDMAX,
					    plusminus[type],
					    genome_p->start[index-1] + (key - genome_p->offset[index-1]));
	  memcpy(&mapping_info_p->errors, r, sizeof(result));
	  mapping_item_p = list_item_new(0, 0, mapping_info_p);
	  list_insert_item(mapping_item_p, &alignment_p->mapping_list);
	  */

/*
	  printf("%s\t%c\t%s %u %s error: %i, pos: %i, base: %i\n",
		 read->id, plusminus[type],
		 index->karyotype.chromosome + (idx-1) * IDMAX,
		 index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]),
		 read->sequence, r->err_kind[0], r->position[0], r->base[0]);

	  num_mappings++;
	}
      }
    }
    found = found || found2;

    free(r_list);
  } // end for type 
  
  if (!found) {
*/
    /*    
    alignment_p->flags = UNMAPPED_FLAG | CODED_SEQ_FLAG;
    alignment_p->header_p = &read_batch_p->header[read_batch_p->header_indices[i]];
    alignment_p->read_p = seq;
    alignment_p->quality_p = &read_batch_p->quality[read_batch_p->data_indices[i]];
    */
    //fprintf(stdout, "%s\n", orig_seq);                                                                                /*                    
//  }

//  return num_mappings;
//}
//*/
/*
//------------------------------------------------------------------------------------

void bwt_seeding_gpu(qreads_t *reads, bwt_optarg_t *bwt_optarg, bwt_context_t *context, 
		     seeding_t *seeding, cal_optarg_t *cal_optarg, 
		     qreads_t *unmapped_reads, alignments_t *mapped_reads);

void bwt_seeding_cpu(qreads_t *reads, bwt_optarg_t *bwt_optarg, bwt_context_t *context, 
		     seeding_t *seeding, cal_optarg_t *cal_optarg,
		     qreads_t *unmapped_reads, alignments_t *mapped_reads);

//------------------------------------------------------------------------------------

void bwt_search(qreads_t *reads, bwt_optarg_t *optarg, bwt_context_t *context, 
		qreads_t *unmapped_reads, alignments_t *mapped_reads);

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

*/
