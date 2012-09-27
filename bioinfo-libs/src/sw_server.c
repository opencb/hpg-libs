#include "sw_server.h"

//====================================================================================
//  Input structure for Smith-Waterman server
//====================================================================================

void sw_server_input_init(list_t *sw_list, list_t* write_list, unsigned int write_size, 
			  float match, float mismatch, float gap_open, float gap_extend, 
			  float min_score, unsigned int flank_length, genome_t *genome, 
			  sw_server_input_t *input) {
  input->sw_list_p = sw_list;
  input->write_list_p = write_list;
  input->write_size = write_size;
  input->genome_p = genome;

  // Smith-Waterman parameters
  input->match = match;
  input->mismatch = mismatch;
  input->gap_open = gap_open;
  input->gap_extend = gap_extend;
  input->min_score = min_score;

  // CAL
  input->flank_length = flank_length;
}

//====================================================================================
// apply_sw
//====================================================================================
int unmapped_by_score_counter[100];

void apply_sw(sw_server_input_t* input, aligner_batch_t *batch) {


  //  printf("START: apply_sw\n"); 
  int tid = omp_get_thread_num();

  cal_t *cal = NULL;
  array_list_t *cal_list = NULL, *mapping_list = NULL;//, *old_list = NULL, *new_list = NULL;
  fastq_batch_t *fq_batch = batch->fq_batch;

  size_t start, end;
  genome_t *genome = input->genome_p;
     
  size_t flank_length = input->flank_length;

  // SIMD support for Smith-Waterman
  float score, min_score = input->min_score;
  //  size_t curr_depth = 0;
  sw_output_t *sw_output;
  //  sw_simd_input_t *sw_sinput = sw_simd_input_new(SIMD_DEPTH);
  //  sw_simd_output_t *sw_soutput = sw_simd_output_new(SIMD_DEPTH);
  //sw_simd_context_t *context = sw_simd_context_new(input->match, input->mismatch, 
  //						    input->gap_open, input->gap_extend); 

  // for tracking the current read, cal being processed using sw_channel_t
  //sw_channel_t *channel;
  //sw_channel_t sw_channels[SIMD_DEPTH];
  //memset(sw_channels, 0, sizeof(sw_channels));
  
  //size_t header_len, read_len;
  //size_t strands[SIMD_DEPTH], chromosomes[SIMD_DEPTH], starts[SIMD_DEPTH];
  
  size_t index, num_cals;
  size_t total = 0, valids = 0;

  size_t num_seqs = batch->num_targets;

  // set to zero
  batch->num_done = batch->num_to_do;
  batch->num_to_do = 0;

  size_t sw_total = batch->num_done;
  /*
  // for all seqs pending to process !!
  size_t sw_total = 0;
  for (size_t i = 0; i < num_seqs; i++) {
    sw_total += array_list_size(batch->mapping_lists[batch->targets[i]]);
  }
  printf("number of sw to run: %d (vs num_done = %d)\n", sw_total, batch->num_done);
  */

  sw_optarg_t sw_optarg; //= sw_optarg_new(gap_open, gap_extend, matrix_filename);
  sw_optarg.gap_open = input->gap_open;
  sw_optarg.gap_extend = input->gap_extend;
  sw_optarg.subst_matrix['A']['A'] = input->match;    sw_optarg.subst_matrix['C']['A'] = input->mismatch; sw_optarg.subst_matrix['T']['A'] = input->mismatch; sw_optarg.subst_matrix['G']['A'] = input->mismatch;
  sw_optarg.subst_matrix['A']['C'] = input->mismatch; sw_optarg.subst_matrix['C']['C'] = input->match;    sw_optarg.subst_matrix['T']['C'] = input->mismatch; sw_optarg.subst_matrix['G']['C'] = input->mismatch;
  sw_optarg.subst_matrix['A']['G'] = input->mismatch; sw_optarg.subst_matrix['C']['T'] = input->mismatch; sw_optarg.subst_matrix['T']['T'] = input->match;    sw_optarg.subst_matrix['G']['T'] = input->mismatch;
  sw_optarg.subst_matrix['A']['T'] = input->mismatch; sw_optarg.subst_matrix['C']['G'] = input->mismatch; sw_optarg.subst_matrix['T']['G'] = input->mismatch; sw_optarg.subst_matrix['G']['G'] = input->match;

  sw_multi_output_t *output = sw_multi_output_new(sw_total);
  char *q[sw_total], *r[sw_total];
  uint8_t strands[sw_total], chromosomes[sw_total];
  size_t starts[sw_total];
  size_t sw_count = 0, read_indices[sw_total];

  // debugging: to kown how many reads are not mapped by SW score
  int unmapped_by_score[batch->num_mapping_lists];

  // initialize query and reference sequences to Smith-Waterman
  for (size_t i = 0; i < num_seqs; i++) {
    index = batch->targets[i];

    // debugging
    unmapped_by_score[index] = 1;

    cal_list = batch->mapping_lists[index];
    num_cals = array_list_size(cal_list);

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // read index
      read_indices[sw_count] = index;

      // query sequence
      q[sw_count] = &(fq_batch->seq[fq_batch->data_indices[index]]);

      // reference sequence
      cal = array_list_get(j, cal_list);
      start = cal->start - flank_length;
      end = cal->end + flank_length;
      r[sw_count] = calloc(1, end - start + 2);
      genome_read_sequence_by_chr_index(r[sw_count], cal->strand,
					cal->chromosome_id - 1, &start, &end, genome);

      // save some stuff, we'll use them after...
      strands[sw_count] = cal->strand;
      chromosomes[sw_count] = cal->chromosome_id;
      starts[sw_count] = start;

      // increase counter
      sw_count++;
    }

    // free cal_list
    array_list_free(cal_list, cal_free);
    batch->mapping_lists[index] = NULL;
  }

  // run Smith-Waterman
  smith_waterman_mqmr(q, r, sw_total, &sw_optarg, 1, output);

  // debugging
  {
    FILE *fd = fopen("sw.out", "w");
    sw_multi_output_save(sw_total, output, fd);
    fclose(fd);
  }

  // filter alignments by min_score
  for (size_t i = 0; i < sw_total; i++) {

    //    score = output->score_p[i] / (strlen(output->query_map_p[i]) * input->match);
    //    if (score >= min_score) {
    /*
    printf("--------------------------------------------------------------\n");
    printf("Smith-Waterman results:\n");
    printf("id\t%s\n", &(batch->fq_batch->header[batch->fq_batch->header_indices[read_indices[i]]]));
    printf("ref\n%s\n", r[i]);
    printf("query\n%s\n", q[i]);
    printf("map\n%s\n", output->ref_map_p[i]);
    printf("ref: chr = %d, strand = %d, start = %d, len = %d\n", chromosomes[i], strands[i], starts[i], strlen(r[i]));
    printf("query-map-start = %d, ref-map-start = %d\n", 
	   output->query_start_p[i], output->ref_start_p[i]);
    printf("score = %0.2f (min. score = %0.2f)\n", output->score_p[i], min_score);
    printf("--------------------------------------------------------------\n");
    */
    if (output->score_p[i] >= min_score) {
      // valid mappings, 
      //insert in the list for further processing
      index = read_indices[i];
      if (batch->mapping_lists[index] == NULL) {
	mapping_list = array_list_new(1000, 
				      1.25f, 
				      COLLECTION_MODE_ASYNCHRONIZED);
	array_list_set_flag(0, mapping_list);
	
	batch->mapping_lists[index] = mapping_list;
      }


      sw_output = sw_output_new(strands[i],
				chromosomes[i],
				starts[i],
				strlen(r[i]),
				strlen(output->query_map_p[i]),
				output->query_start_p[i],
				output->ref_start_p[i],
				output->score_p[i],
				score,
				output->query_map_p[i],
				output->ref_map_p[i]);
      array_list_insert(sw_output, mapping_list);

      batch->num_to_do++;

      // debugging
      unmapped_by_score[index] = 0;
    }

    // free reference
    free(r[i]);
  }
  /*
  // debugging
  for (size_t i = 0; i < num_seqs; i++) {
    index = batch->targets[i];

    if (unmapped_by_score[index] && (strncmp("@rand", &(batch->fq_batch->header[batch->fq_batch->header_indices[index]]), 5))) {
	unmapped_by_score_counter[tid]++;
	printf("by score: %s\n", &(batch->fq_batch->header[batch->fq_batch->header_indices[index]]));
      }
  }
  */

  // update counter
  thr_sw_items[tid] += batch->num_done;

  // free
  sw_multi_output_free(output);

  //  printf("END: apply_sw, (%d Smith-Waterman, %d valids)\n", total, valids);
}

//--------------------------------------------------------------------------------------

sw_output_t *sw_output_new(int strand, size_t chrom, size_t ref_start, size_t ref_len,
			   size_t mref_len, size_t mquery_start, size_t mref_start,
			   float score, float norm_score, char* mquery, char* mref) {

  sw_output_t *p = (sw_output_t *) calloc(1, sizeof(sw_output_t));

  p->strand = strand;
  p->chromosome = chrom;
  p->ref_start = ref_start;
  p->ref_len = ref_len;
  p->mref_len = mref_len;
  p->mquery_start = mquery_start;
  p->mref_start = mref_start;
  p->score = score;
  p->norm_score = norm_score;
  p->mquery = strdup(mquery);
  p->mref = strdup(mref);
  //p->mquery = NULL;
  //p->mref = NULL;

  return p;
}

//--------------------------------------------------------------------------------------

void sw_output_free(sw_output_t *p) {
  if (p == NULL) return;

  if (p->mquery != NULL) free(p->mquery);
  if (p->mref != NULL) free(p->mref);

  free(p);
}

//--------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

