#include <omp.h>

#include "rna_server.h"

#define MINIMUN_CAL_LENGTH 10
#define MINIMUM_INTRON_LENGTH 10
#define ERRORS_ZONE 8

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

//======== SPLICE JUNCTION TYPE ============
#define NOT_SPLICE	0
#define GT_AG_SPLICE  	1
#define CT_AC_SPLICE  	2
//==========================================

#define SEARCH_START 0
#define SEARCH_END 1
#define POSSIBLES_MARKS 2

#define STRANDS_NUMBER 2


void seq_reverse_complementary(char *seq, unsigned int len){
  unsigned int j = 0;
  
  char *seq_tmp = (char *)malloc(len*sizeof(char));
  memcpy(seq_tmp, seq, len);
  
  for (int i = len - 1; i >= 0; i--){
    if (seq_tmp[i] == 'A' || seq_tmp[i] == 'a'){
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
    }
    j++;  
  }

  free(seq_tmp);

}

void cal_fusion_data_init(unsigned int id, unsigned int start, unsigned int end, unsigned int strand, unsigned int chromosome, unsigned int fusion_start, unsigned int fusion_end, cal_fusion_data_t *cal_fusion_data_p){
  cal_fusion_data_p->genome_strand = strand;
  cal_fusion_data_p->id = id;
  cal_fusion_data_p->genome_start = start;
  cal_fusion_data_p->genome_chromosome = chromosome;
  cal_fusion_data_p->genome_end = end;
  cal_fusion_data_p->fusion_start = fusion_start;
  cal_fusion_data_p->fusion_end = fusion_end;
}

char* cigar_automata_status(unsigned char status){
  char *str_status;
  
  switch (status){
    case CIGAR_MATCH_MISMATCH:
      str_status = "MIS/MATCH";
      break;
    case CIGAR_INSERTION:
      str_status = "INSERTION";
      break;
    case CIGAR_DELETION:
      str_status = "DELETION";
      break;
    case CIGAR_BIG_DELETION:
      str_status = "BIG DELETION";
      break;
    case CIGAR_SKIPPED:
      str_status = "SKIPPED";
      break;
    case CIGAR_PADDING:
      str_status = "PADDING";
    default:
      str_status = "UNKNOWN";
      break;
  }
  return str_status;
}

void cigar_generator(cigar_data_t *cigar_p, unsigned int *max_len, unsigned int *pos, unsigned int status, unsigned int *number){
  char cigar_info[20];
  
  if (*number <= 0){
    return;
  }else{
    //printf("**********Allocate %d\n",*number );
    if ((status != CIGAR_DELETION) || ((status == CIGAR_DELETION) && (*number <= MINIMUM_INTRON_LENGTH))){
	  cigar_p[*pos].type = status;
	  cigar_p[*pos].value = *number;
	  *pos += 1;
	  *number = 1;
	  if (*pos >= *max_len){
	    *max_len = *max_len * 2;
	    cigar_p = (cigar_data_t *)realloc(cigar_p, sizeof(cigar_data_t)*(*max_len));
	  }

      }
  }

  return;
  
}

write_batch_t* search_splice_junctions_sw_output(sw_simd_input_t* input_p, sw_simd_output_t* output_p, unsigned int depth, allocate_fusion_data_t *depth_cal_fusion_p, 
						 allocate_splice_elements_t *chromosome_avls_p,  unsigned char * mapping_reads_p, sw_channel_t *sw_channels_p, sw_batch_t *sw_batch_p,
						 list_t* write_list_p, write_batch_t* write_batch_p, 
						 unsigned int write_size, unsigned int sw_id,  size_t *sw_no_valids, float min_score, genome_t *genome_p, unsigned int min_intron_length){
  
  const unsigned char EXTRA_SEARCH = 8;

  char header_id[2048];
  alignment_t *alignment_p;
  
  char *header_match, *read_match, *quality_match;
  
  char extra_nt[EXTRA_SEARCH + 1];
  unsigned int extra_nt_len = 0;
  unsigned char ref_nt_pos;
  unsigned char extra_nt_found;
  int gap_start = -1; int gap_end = -1;
  unsigned char found = 0;
  unsigned char confirmation_type = 0;
  unsigned int l;
  unsigned char count_extra_search;
  unsigned int position;
  unsigned int start_splice, end_splice;
  unsigned char actual_cal = 0;
  unsigned int start_mapped = 0;
  unsigned int n_deletions_seq = 0;
  unsigned int insertion_number;
  unsigned int deletion_number;
  unsigned int deletions_tot;
  unsigned int displacement_start;
  unsigned int mark_dpl;
  unsigned int strand;
  unsigned char strand_splice;
  int cal;
  unsigned int relative_gap_start;
  unsigned int nt_skipped;
  unsigned int str_pos, quality_pos;
  unsigned int reference_position_disp;
  unsigned int start_gap_padding;
  unsigned int start_search_j;

  unsigned int extend_start_sp, extend_end_sp;

  unsigned int splice_number;
  sp_data_t allocate_splice[5];

  //=========CIGAR INFO. && BAM/SAM ========//
  unsigned int cigar_soft;
  unsigned int value;
  int cigar_value;
  unsigned int cigar_max_len = 200;
  unsigned int cigar_pos;
  unsigned char automata_status;
  char flag;
  short int num_cigar_op;
  unsigned int header_len, read_len, mapped_len, bytes;
  short int primary_alignment;
  unsigned int total_cals;
  char end_exceded = 0;
  unsigned int orig_pos = 0;
  unsigned int tmp;
  unsigned int real_sw_len = 0;
  unsigned int j_start;
  unsigned int change_cal = 0, change_cal_number = 0;
  char cigar_segment[10];
  unsigned int change_cal_index[20];
  list_item_t* item_p = NULL;  
  char *cigar_str;
  cigar_data_t *cigar_p = (cigar_data_t *)malloc(sizeof(cigar_data_t)*cigar_max_len);
  if (cigar_p == NULL) { exit(-1); }
  //============================//
  
  
  //============ CANONICAL INTRON MARKS ==============//
  unsigned char canonical_GT_AG[4];
  unsigned char canonical_CT_AC[4];
  
  //START INTRON MARKS 'GT' or 'CT'
  canonical_GT_AG[0] = 'G';
  canonical_GT_AG[1] = 'T';
  canonical_CT_AC[0] = 'C';
  canonical_CT_AC[1] = 'T';
  //END INTRON MARKS 'AG' or 'AC'
  canonical_GT_AG[2] = 'A';
  canonical_GT_AG[3] = 'G';
  canonical_CT_AC[2] = 'A';
  canonical_CT_AC[3] = 'C';
  //================================================//
  
  /*printf(" ======================== Process Output SW %d=========================\n", depth);
  sw_simd_input_display(depth, input_p);
  sw_simd_output_display(depth, output_p);
  printf("======================================================================\n");
  */

  for (int i = 0; i < depth; i++) {    
    splice_number = 0;
    total_cals = depth_cal_fusion_p[i].cal_number;

    //printf("Total cals %d\n", total_cals);

    if (total_cals == 1) {
      if (output_p->score_p[i] < min_score) {
	*sw_no_valids += 1;
	/*free(output_p->mapped_seq_p[i]);
	  free(output_p->mapped_ref_p[i]);*/
	continue;
      }
    } else if (output_p->score_p[i] < (min_score - 50.0)) {
      *sw_no_valids += 1;
      /*free(output_p->mapped_seq_p[i]);
	free(output_p->mapped_ref_p[i]);*/
      continue;
    }
    /*
    printf("=========================== BEFORE [Depth %d - Total CALs %d]==================================\n", i, depth_cal_fusion_p[i].cal_number);
    for(int c = 0; c < depth_cal_fusion_p[i].cal_number; c++){
      printf("CAL %d :\n", c);
      printf("\t->Genome Start : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_start);
      printf("\t->Genome End : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_end);
      printf("\t->Fusion Start : %d\n", depth_cal_fusion_p[i].allocate_data[c].fusion_start);
      printf("\t->Fusion End : %d\n", depth_cal_fusion_p[i].allocate_data[c].fusion_end);
      printf("\t->Chromosome : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_chromosome);
      printf("\t->Strand : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_strand);
    }
    printf("==============================================================================================\n");
    */

    reference_position_disp = 0;

    //printf("Actualization start %d\n", output_p->start_p[i]);
    //====================== Actualization CALS ===========================//
    //Move to init CAL and calculate start padding
    cal = 0;
    displacement_start = 0;
    while (depth_cal_fusion_p[i].allocate_data[cal].fusion_end  <= output_p->start_p[i]) {
      displacement_start = depth_cal_fusion_p[i].allocate_data[cal].fusion_end + 1;
      depth_cal_fusion_p[i].allocate_data[cal].fusion_end = 0;
      depth_cal_fusion_p[i].allocate_data[cal].fusion_start = 0;
      cal++;
      //printf("Delete CAL %d\n", cal);
      if (cal >= total_cals) {
	printf("ERROR: total cals overflow\n");
	continue;
      }
    }

    actual_cal = cal;
    
    displacement_start = output_p->start_p[i] - displacement_start;

    depth_cal_fusion_p[i].allocate_data[cal].genome_start += displacement_start;
        
    start_mapped = depth_cal_fusion_p[i].allocate_data[cal].genome_start;
    start_gap_padding = displacement_start;

    depth_cal_fusion_p[i].allocate_data[cal].fusion_start = 0;
    depth_cal_fusion_p[i].allocate_data[cal].fusion_end = depth_cal_fusion_p[i].allocate_data[cal].fusion_end - output_p->start_p[i];
    cal++;
    while (cal < total_cals) {
      depth_cal_fusion_p[i].allocate_data[cal].fusion_start = depth_cal_fusion_p[i].allocate_data[cal].fusion_start - output_p->start_p[i];
      depth_cal_fusion_p[i].allocate_data[cal].fusion_end = depth_cal_fusion_p[i].allocate_data[cal].fusion_end - output_p->start_p[i];
      cal++;
    }
    
    cal = total_cals - 1;

    while ((depth_cal_fusion_p[i].allocate_data[cal].fusion_start + 1) > output_p->mapped_len_p[i]) {
      depth_cal_fusion_p[i].allocate_data[cal].fusion_end = 0;
      depth_cal_fusion_p[i].allocate_data[cal].fusion_start = 0;
      cal--;
      if (cal < 0) {
	printf("ERROR: total cals overflow\n");
	continue;
      }
    }

    depth_cal_fusion_p[i].allocate_data[cal].genome_end -= (depth_cal_fusion_p[i].allocate_data[cal].fusion_end - output_p->mapped_len_p[i]);
    depth_cal_fusion_p[i].allocate_data[cal].fusion_end  = output_p->mapped_len_p[i] - 1;
    

    
    j_start = 0;
    while (((depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start) < 
	   MINIMUN_CAL_LENGTH)
	   && (actual_cal < total_cals)) {
      j_start += (depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start) + 1;
      depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start = 0;
      depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end = 0;
      actual_cal++;
      if (actual_cal >= total_cals) {
	printf("ERROR: total cals overflow\n");
	continue;
      }
      start_mapped = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start;      
      //printf("Move to next CAL %i, start_j = %i\n", actual_cal, j_start);
    }

    cal = actual_cal;
    change_cal_number = 0;
    real_sw_len = 0;
    memset( change_cal_index, 0, sizeof(unsigned int)*total_cals);

    
    /*printf("Actual CAL %d", actual_cal);
    printf("=========================== AFTER [Depth %d - Total CALs %d]==================================\n", i, depth_cal_fusion_p[i].cal_number);
    for(int c = 0; c < depth_cal_fusion_p[i].cal_number; c++){
      printf("CAL %d :\n", c);
      printf("\t->Genome Start : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_start);
      printf("\t->Genome End : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_end);
      printf("\t->Fusion Start : %d\n", depth_cal_fusion_p[i].allocate_data[c].fusion_start);
      printf("\t->Fusion End : %d\n", depth_cal_fusion_p[i].allocate_data[c].fusion_end);
      printf("\t->Chromosome : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_chromosome);
      printf("\t->Strand : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_strand);
    }
    printf("=============================================================================================================\n");*/
    
    //====================== Actualization CALS End ===========================//

    //Data Info
    position = 0;
    count_extra_search = 0;
    start_splice = 0;
    end_splice = 0;
    insertion_number = 0;
    deletion_number = 0;
    mark_dpl=0;
    relative_gap_start = 0;
    n_deletions_seq = 0;
    orig_pos = 0;
    deletions_tot = 0;    
    str_pos = 0;
    read_len = output_p->mapped_len_p[i];
    read_match = (char *)malloc(sizeof(char)*(read_len + 1));
    if (read_match == NULL) { exit(-1); }
   
    quality_pos = output_p->start_seq_p[i];
    quality_match = (char *)malloc(sizeof(char)*(read_len + 1));
    if (quality_match == NULL) { exit(-1); }
    
    found = NOT_SPLICE;
    //===========Position Automata in initial State==========
    cigar_value = 0;
    cigar_pos = 0;
    
    if (output_p->mapped_ref_p[i][0] != '-') {
      if (output_p->mapped_seq_p[i][0] != '-') {
	automata_status = CIGAR_MATCH_MISMATCH;
      } else {
	automata_status = CIGAR_DELETION;
      }
    } else {
      automata_status = CIGAR_INSERTION;
    }
    //========================================================
    
    change_cal = 0;
    for (int j = j_start; j < output_p->mapped_len_p[i]; j++) {

      if (output_p->mapped_ref_p[i][j] != '-') {
	reference_position_disp++;
	start_gap_padding++;
	relative_gap_start++;
	//TODO:IMPLEMENTS WHEN JUMP TO OTHER CALL WITHOUT GAP
	if (output_p->mapped_seq_p[i][j] != '-') {
	  //Match mismatch area	  
	  if (automata_status == CIGAR_MATCH_MISMATCH) {
	    cigar_value++;
	  } else {
	    cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value);
	    automata_status = CIGAR_MATCH_MISMATCH;
	  }
	}
      } else {
	//Insertion area
	orig_pos++;
	insertion_number++;
	if (automata_status == CIGAR_INSERTION) {
	   cigar_value++;
	} else {
	    cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value);
	    automata_status = CIGAR_INSERTION;
	}
      }//end mapped_ref_p[i][j] != '-'
       

      if (output_p->mapped_seq_p[i][j] == '-') {
	if (gap_start == -1) {
	  gap_start = j;
	}

	n_deletions_seq++;
	deletions_tot++;
	
	if (output_p->mapped_ref_p[i][j] == '-') {
	  //Padding area
	  if (automata_status == CIGAR_PADDING) {
	    cigar_value++;
	  } else {
	    cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value);
	    automata_status = CIGAR_PADDING;    
	  }
	} else {
	  //Deletion area
	  if (automata_status == CIGAR_DELETION) {
	    cigar_value++;
	  } else {
	    cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value);
	    automata_status = CIGAR_DELETION;
	  }
	}
      } else { 
	 read_match[str_pos] = output_p->mapped_seq_p[i][j];
	 quality_match[str_pos++] = sw_batch_p->allocate_reads_p[sw_channels_p[i].read_index]->quality[quality_pos];
	 quality_pos++;
	 
	 if (n_deletions_seq > MINIMUM_INTRON_LENGTH) {
	    
	   gap_end = j - 1;   
	   //printf("gap_start=%d , gap_end=%d\n", gap_start, gap_end);

	   //Canonical marks strand + (GT-AG)
	   if ((output_p->mapped_ref_p[i][gap_start] == canonical_GT_AG[0]) && (output_p->mapped_ref_p[i][gap_start + 1] == canonical_GT_AG[1])) {
	     if ((output_p->mapped_ref_p[i][gap_end - 1] == canonical_GT_AG[2]) && (output_p->mapped_ref_p[i][gap_end] == canonical_GT_AG[3])) {
	       found = GT_AG_SPLICE; 
	     }
	   }
	   //Canonical marks strand - (CT-AC)
	   if ((output_p->mapped_ref_p[i][gap_start] == canonical_CT_AC[0]) && (output_p->mapped_ref_p[i][gap_start + 1] == canonical_CT_AC[1])) {
	     if ((output_p->mapped_ref_p[i][gap_end - 1] == canonical_CT_AC[2]) && (output_p->mapped_ref_p[i][gap_end] == canonical_CT_AC[3])) {
	       found = CT_AC_SPLICE;
	     }
	   }
	   
	   if (!found) {
	     //Search Xnt (+)---->	    
	     confirmation_type = NOT_SPLICE;
	     count_extra_search = 1;
	     position = gap_start + 1;
	     extra_nt[count_extra_search - 1] = output_p->mapped_ref_p[i][position - 1];
	     //printf("Extra Search position %d and gap end %d\n", position, gap_end);
	     while (count_extra_search < EXTRA_SEARCH ) {
	       if ((gap_end + count_extra_search) >= output_p->mapped_len_p[i]) {
		 end_exceded = 1; break;
	       }

	       if ((output_p->mapped_ref_p[i][position] == canonical_GT_AG[0]) && (output_p->mapped_ref_p[i][position + 1] == canonical_GT_AG[1])) {
		 //printf("-Found GT-");
		 if ((output_p->mapped_ref_p[i][gap_end - 1 + count_extra_search] == canonical_GT_AG[2]) && (output_p->mapped_ref_p[i][gap_end + count_extra_search] == canonical_GT_AG[3])) {
		   confirmation_type = GT_AG_SPLICE; 
		   //printf("Found GT_AG Extra Search\n");
		 }
		 
	       } else if ((output_p->mapped_ref_p[i][position] == canonical_CT_AC[0]) && (output_p->mapped_ref_p[i][position + 1] == canonical_CT_AC[1])) {
		 if ((output_p->mapped_ref_p[i][gap_end - 1 + count_extra_search] == canonical_CT_AC[2]) && (output_p->mapped_ref_p[i][gap_end + count_extra_search] == canonical_CT_AC[3])) {
		   confirmation_type = CT_AC_SPLICE; 
		   //printf("Found CT_AC Extra Search\n");
		 }
	       }
	       
	       if (confirmation_type) {
		 extra_nt[count_extra_search] = '\0';
		 extra_nt_len = strlen(extra_nt);
		 ref_nt_pos = gap_end + 1;
		 extra_nt_found = 1;
		 
		 for (int k = 0; k < extra_nt_len; k++) {
		   if (extra_nt[k] != output_p->mapped_ref_p[i][ref_nt_pos]) {
		     extra_nt_found = 0;
		     break;
		   }
		   ref_nt_pos++;
		 }
		
		 if (extra_nt_found) {
		   start_splice = position;
		   end_splice = gap_end + count_extra_search - 1;
		   found= confirmation_type;
		   mark_dpl = count_extra_search;
		   confirmation_type = NOT_SPLICE;
		 }
	       }
	       
	       count_extra_search++;
	       position++;
	       extra_nt[count_extra_search - 1] = output_p->mapped_ref_p[i][position - 1];
	     }//End while extra search
	     //printf("Extra search end\n");
	   }//End if not found
	  
	   //TODO: If not found make new search <-----(-) and consider that in middle of two CALs with splice junction not found other CAL.
	   
	   if (found) {
	    
	     //printf("====================Found Splice Junction %d Splice number========================\n", splice_number);
	
	     relative_gap_start--;
	     
	     relative_gap_start -= n_deletions_seq;//((gap_end - gap_start) + 1);
	     start_splice = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start + (relative_gap_start + mark_dpl);
	     extend_start_sp = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start;
	     
	     allocate_splice[splice_number].start_sp = start_splice;
	
	     //printf("Before.Value %d\n", allocate_splice[0].start_extend_sp);     

	     if (splice_number >= 1) {
	      allocate_splice[splice_number].start_extend_sp = allocate_splice[splice_number - 1].end_sp;
	      allocate_splice[splice_number - 1].end_extend_sp = start_splice;
	     } else {
	       allocate_splice[splice_number].start_extend_sp = extend_start_sp;
	     }

	     //printf("After.Value %d\n", allocate_splice[0].start_extend_sp);     

	    //printf("(+)START SPLICE :: %d + (%d + %d) :: %d\n", depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start, relative_gap_start, mark_dpl, start_splice);
	   
	    if (gap_end >  depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
	      //printf("\tChange to Next CAL\n");
	     
	      tmp = ((depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - j_start) - (depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start - j_start) + 1) - relative_gap_start ;
	      
	      relative_gap_start = n_deletions_seq - tmp;
	     
	      actual_cal++;
	      nt_skipped = 0;

	      if (gap_end >  depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
		while (gap_end >  depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
		  nt_skipped += (depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end -depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start);
		  actual_cal++;
		}
		nt_skipped++;
	      }
	      
	      if (actual_cal >= total_cals){ printf("ERROR:CAL overflow\n");end_exceded = 1; break; }
	      
	      relative_gap_start -= nt_skipped;
	      
	      start_gap_padding = 0;
	      
	      end_splice = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start + relative_gap_start + mark_dpl;
	      end_splice--;
	
	      relative_gap_start++;
	      insertion_number = 0;
	    } else {
	      //printf("\tStay in same CAL\n");
	      relative_gap_start += n_deletions_seq;
	      
	      end_splice = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start + (relative_gap_start - 1) + mark_dpl;
	    }

	    extend_end_sp = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_end;
	    allocate_splice[splice_number].end_sp = end_splice;

	    if (found == CT_AC_SPLICE) {
	      strand_splice = 1;
	    } else {
	      strand_splice = 0;
	    }

	    allocate_splice[splice_number].strand_sp = strand_splice;
	    allocate_splice[splice_number].end_extend_sp = extend_end_sp;
	    cigar_value = (int)end_splice - (int)start_splice + 1;
	    
	    if (cigar_value > min_intron_length) {
	      splice_number++;
	      cigar_p[cigar_pos - 1].value += mark_dpl;	    
	      cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, CIGAR_SKIPPED, &cigar_value);
	      cigar_value -= mark_dpl;
	      /*if ((start_splice == 89074710) && ( end_splice == 89075610)){
		printf("%d.Seq: %s\n", i, input_p->seq_p[i]);
		exit(-1);
		}*/
	    } else {
	      cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, CIGAR_BIG_DELETION, &cigar_value);
	    }
	    
	    found = NOT_SPLICE;
	    deletion_number = 0;
	    //printf("==================== Exact->( %d:%d-%d ) Extend->( %d:%d-%d )  ========================\n", depth_cal_fusion_p[i].allocate_data->genome_chromosome, start_splice, end_splice, depth_cal_fusion_p[i].allocate_data->genome_chromosome, extend_start_sp, extend_end_sp);
	  } else { //IF not found
	    if (reference_position_disp > depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
	      actual_cal++;
	      
	      if (actual_cal >= total_cals) {
		//Cals overflow, read no valid		
		end_exceded = 1;
		break;
	      } else {
		relative_gap_start = gap_end   - (depth_cal_fusion_p[i].allocate_data[actual_cal - 1].fusion_end + insertion_number);
		
		cigar_value = (depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start + relative_gap_start) -  (depth_cal_fusion_p[i].allocate_data[actual_cal - 1].genome_end - ((gap_end - gap_start) - relative_gap_start));
	      }
	      cigar_generator(cigar_p, &cigar_max_len, &cigar_pos,  CIGAR_BIG_DELETION, &cigar_value);
	    } else {
	     cigar_value = gap_end - gap_start + 1;
	     cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, CIGAR_BIG_DELETION, &cigar_value);
	    }

	    deletion_number = 0;
	    insertion_number = 0;
	    n_deletions_seq = 0;

	   //CIGAR_BIG_DELETION
	  }

	 } else {
	   //Deletion
	   deletion_number += n_deletions_seq;
	   n_deletions_seq = 0;
	 }//if minimum deletions mark 
	 gap_start = -1;
	 n_deletions_seq = 0;
      }//Not found
    }//End loop process smith-waterman
    //printf("End process sw %i\n", i);
    if (end_exceded) { end_exceded = 0; continue; }
    //printf("\tEnd Search Splice\n");
    //=============================== MAKE CIGAR & SAM/BAM =====================================
    cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value);
    if (cigar_pos < 1) { continue; }

    //Generate cigar string
    num_cigar_op = 0;
    cigar_str = (char *)malloc(sizeof(char)*cigar_max_len);
    if (cigar_str == NULL) { exit(-1); }
    cigar_str[0] = '\0';
    
    //Hard and Soft clipped start
    if (output_p->start_seq_p[i] > 0) {
      sprintf(cigar_segment, "%iH", output_p->start_seq_p[i]);
      cigar_str = strcat(cigar_str, cigar_segment);
      num_cigar_op++;
    }
    
    if (cigar_p[0].type == CIGAR_MATCH_MISMATCH) {
      cigar_soft = 0;
      value = 0;
      while ((output_p->mapped_ref_p[i][cigar_soft] != '-') && (output_p->mapped_seq_p[i][cigar_soft] != '-') && 
	      (output_p->mapped_seq_p[i][cigar_soft] != output_p->mapped_ref_p[i][cigar_soft])) {
		cigar_soft++;
		value++;
		//printf("(Soft %c!=%c)", output_p->mapped_ref_p[i][cigar_soft], output_p->mapped_seq_p[i][cigar_soft]);
      }
      
      if (value > 0) {
	cigar_p[0].value -= cigar_soft;
	sprintf(cigar_segment, "%iS", value);
	cigar_str = strcat(cigar_str, cigar_segment);
	num_cigar_op++;
      }
    }
    
    //Hard and Soft clipped end
    //printf("Cigar value %d\n", cigar_pos);
    value = 0;
    if (cigar_p[cigar_pos - 1].type == CIGAR_MATCH_MISMATCH) {
      cigar_soft = output_p->mapped_len_p[i] - 1;
      while ((output_p->mapped_ref_p[i][cigar_soft] != '-') && (output_p->mapped_seq_p[i][cigar_soft] != '-') && 
	      (output_p->mapped_seq_p[i][cigar_soft] != output_p->mapped_ref_p[i][cigar_soft])) {
		cigar_soft--;
		value++;
		//printf("(Soft %c!=%c)", output_p->mapped_ref_p[i][cigar_soft], output_p->mapped_seq_p[i][cigar_soft]);
      }
      cigar_p[cigar_pos - 1].value -= value;
    }
    
    for (int cig = 0; cig < cigar_pos; cig++) {
      if (cigar_p[cig].type == CIGAR_MATCH_MISMATCH) {
	flag = 'M';
      } else if (cigar_p[cig].type == CIGAR_DELETION || cigar_p[cig].type == CIGAR_BIG_DELETION) {
	flag = 'D';
      } else if (cigar_p[cig].type == CIGAR_PADDING) {
	flag = 'P';
      } else if (cigar_p[cig].type == CIGAR_INSERTION) {
	flag = 'I';
      } else if (cigar_p[cig].type == CIGAR_SKIPPED) {
	flag = 'N';
      }

      sprintf(cigar_segment, "%i%c", cigar_p[cig].value, flag);
      cigar_str = strcat(cigar_str, cigar_segment);
      num_cigar_op++;
      //printf("%d%c", cigar_p[cig].value, flag);
    }
    
    
    
    if (value > 0) {
      sprintf(cigar_segment, "%iS", value);
      cigar_str = strcat(cigar_str, cigar_segment);
      num_cigar_op++;
     }
    
    if (((output_p->mapped_len_p[i] - deletions_tot) + output_p->start_seq_p[i]) < input_p->seq_len_p[i]) {
      sprintf(cigar_segment, "%iH", input_p->seq_len_p[i] - ((output_p->mapped_len_p[i] - deletions_tot) + output_p->start_seq_p[i]));
      cigar_str = strcat(cigar_str, cigar_segment);
      num_cigar_op++;
    }
    
    //printf("Cigar(%d):%s\n", cigar_pos, cigar_str);
    //end cigar string
    
    if ( write_batch_p->size >= write_batch_p->allocated_size - 1) {
      item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
      if (time_on) { timing_stop(SW_SERVER, sw_id, timing_p); }
      list_insert_item(item_p, write_list_p);
      if (time_on) { timing_start(SW_SERVER, sw_id, timing_p); }
      
      write_batch_p = write_batch_new(write_size, MATCH_FLAG);
      
    }
    
    if (mapping_reads_p[sw_channels_p[i].read_index] == 1) {
	primary_alignment = 1;
    } else {
        primary_alignment = 0;
    }
    
    alignment_p = alignment_new();
    header_len = sw_channels_p[i].header_len;
    //printf("Reallocate %d\n", header_len); 
    get_to_first_blank(sw_batch_p->allocate_reads_p[sw_channels_p[i].read_index]->id, header_len, header_id);
    header_match = (char *)malloc(sizeof(char)*header_len);
    if (header_match == NULL) { exit(-1); }
    memcpy(header_match, header_id, header_len);
    
    //printf("Header len = %d\n", header_len);

    read_match[str_pos] = '\0';
    quality_match[str_pos] = '\0';
    
    //memcpy(quality_match, &(sw_batch_p->quality_p[sw_batch_p->read_indices_p[sw_channels_p[i].read_index]]), read_len );

    alignment_init_single_end(header_match, read_match, quality_match, depth_cal_fusion_p[i].allocate_data->genome_strand, depth_cal_fusion_p[i].allocate_data->genome_chromosome - 1, start_mapped, cigar_str, num_cigar_op, 255, 1, primary_alignment, alignment_p);
    //printf("seq: %s\n", alignment_p->sequence);
    
    ((alignment_t **)write_batch_p->buffer_p)[write_batch_p->size] = alignment_p;
    write_batch_p->size++;
	     
    //============================================================================================
    
    //printf(" (score = %.2f, norm. score = %.2f)\n", output_p->score_p[i], output_p->norm_score_p[i]);
   
    //printf("Free done 1\n");
    /*
      free(output_p->mapped_seq_p[i]);
      free(output_p->mapped_ref_p[i]);
    */
    mapping_reads_p[sw_channels_p[i].read_index] = 1;
    
    //Report Splice Junctions
    for (unsigned int s = 0; s < splice_number; s++) {
      pthread_mutex_lock(&(chromosome_avls_p[depth_cal_fusion_p[i].allocate_data[0].genome_chromosome].mutex));
      allocate_new_splice(depth_cal_fusion_p[i].allocate_data->genome_chromosome, allocate_splice[s].strand_sp, allocate_splice[s].end_sp, allocate_splice[s].start_sp, allocate_splice[s].start_extend_sp, allocate_splice[s].end_extend_sp, chromosome_avls_p);
      pthread_mutex_unlock(&(chromosome_avls_p[depth_cal_fusion_p[i].allocate_data[0].genome_chromosome].mutex));
    }
  }//End loop depth
  
  for (int i = 0; i < depth; i++) {
    free(depth_cal_fusion_p[i].allocate_data);
    free(input_p->seq_p[i]);
  }

  for (unsigned int i = 0; i < 4; i++) {
    free(output_p->mapped_seq_p[i]);
    free(output_p->mapped_ref_p[i]);
  }

  free(cigar_p);

  return write_batch_p;
  
}

void rna_server_omp_smith_waterman(sw_server_input_t* input_p, allocate_splice_elements_t *chromosome_avls_p){
  /**----------------------------------------------------------------------------------**
   **                                    RNA-SEQ                                       **
   **----------------------------------------------------------------------------------**/
  
  printf("rna_server (%i): START\n", omp_get_thread_num());
  unsigned int seed_max_distance = input_p->seed_max_distance;

  unsigned int j, z, k, num_reads, read, bytes, num_cals, cals_positive, cals_negative;
  size_t i;

  //lists and items
  list_t* sw_list_p = input_p->sw_list_p;
  list_item_t *sw_item_p = NULL;
  sw_batch_t* sw_batch_p = NULL;
  
  list_item_t* item_p = NULL;
	
  list_t* write_list_p = input_p->write_list_p;
  unsigned int write_size = input_p->write_size;
  write_batch_t* write_batch_p = write_batch_new(write_size, MATCH_FLAG);
  
  int sw_id = omp_get_thread_num();

  // genome
  char* ref_p;
  unsigned int ref_len;
  
  genome_t* genome_p = input_p->genome_p;
  unsigned int min_intron_length = input_p->min_intron_size;
  
  // SIMD support for Smith-Waterman
  unsigned int depth = 4, curr_depth = 0;
  sw_simd_input_t* sw_input_p = sw_simd_input_new(depth);
  sw_simd_output_t* sw_output_p = sw_simd_output_new(depth);
  sw_simd_context_t *context_p = sw_simd_context_new(input_p->match, input_p->mismatch, input_p->gap_open, input_p->gap_extend);
  float min_score = input_p->min_score;
	
  // for tracking what reads are being mapped successfully
  size_t allocated_mapping_reads = 10000;
  size_t maximum_allocate = 100;
  unsigned char* mapping_reads_p = (unsigned char*) calloc (allocated_mapping_reads, sizeof(unsigned char));
  if (mapping_reads_p == NULL) { exit(-1); }

  // for tracking the current read, cal being processed using sw_channel_t
  sw_channel_t *channel_p, *sw_channels_p = (sw_channel_t*) calloc(depth, sizeof(sw_channel_t));

  unsigned int header_len, read_len;
  unsigned int total = 0, total_valids = 0, total_reads = 0;
  
  cal_t cal_tmp;
  unsigned int *cals_per_chromosome = (unsigned int *) malloc (CHROMOSOME_NUMBER * STRANDS_NUMBER * sizeof(unsigned int));
  if (cals_per_chromosome == NULL) { exit(-1); }

  unsigned int chromosome_tot = 0, cal_chromosome_start = 0, cal_chromosome_end = 0;
  char read_seq[2000];
  char *seq_p;  
  unsigned long int start, end, start_before, end_before;
  unsigned int len;
  char *reference, header_id[1024];
  unsigned int maximum_reference_len = 1024;
  unsigned int len_increment = 512;
  //unsigned int pos_strand;
  char action;
  unsigned int reference_fusion_len = 0;
  char *reference_fusion_p = (char *)malloc(sizeof(char)*maximum_reference_len);
  if (reference_fusion_p == NULL) { exit(-1); }

  unsigned int chromosome;
  unsigned int cal_indice;
  allocate_fusion_data_t *cals_fusion_p, *allocate_cals_fusion_p = (allocate_fusion_data_t *)calloc(depth, sizeof(allocate_fusion_data_t));
  if (allocate_cals_fusion_p == NULL) { exit(-1); }

  char *header_id_match_p, *search_match_p, *quality_match_p, *cigar_p;
  unsigned int len_id_header;
  alignment_t *alignment_p;
  cal_t *cal_prev, *cal_next, *cal_p;
  short int strand;
  unsigned int cals_offset, start_incr, pos_strand;
  unsigned int chromosome_id;
  unsigned int flank_length = input_p->flank_length;
  size_t num_batches = 0, total_reads_process = 0, total_reads_unmapped = 0, num_sw_process = 0;
  size_t sw_no_valids = 0;
  size_t max_reference = 2048;
  unsigned int chr_before;
  unsigned int num_cals_fusion;
  //  int distance;
  char no_store;
  size_t distance;
  unsigned int max_intron_size = input_p->max_intron_size;

  reference = (char *)malloc(sizeof(char)*max_reference);
  
  if (reference == NULL) { exit(-1); }

  // start main loop
  while ( (sw_item_p = list_remove_item(sw_list_p)) != NULL ) {
    num_batches++;
    //printf("RNA Server Processing Batch...\n");
    if (time_on) { timing_start(SW_SERVER, sw_id, timing_p); }
    
    curr_depth = 0;
    sw_batch_p = (sw_batch_t*) sw_item_p->data_p;
    num_reads = sw_batch_p->num_reads;
    total_reads += num_reads;
    
    if (num_reads > allocated_mapping_reads) {
      allocated_mapping_reads = num_reads;
      mapping_reads_p = (unsigned char *) realloc(mapping_reads_p, allocated_mapping_reads * sizeof(unsigned char));
    }

    memset(mapping_reads_p, 0, allocated_mapping_reads * sizeof(unsigned char));

    /** for each read **/
    total_reads_process += num_reads;
    for (i=0; i < num_reads; i++) {
      read_len = strlen(sw_batch_p->allocate_reads_p[i]->sequence);
      header_len = strlen(sw_batch_p->allocate_reads_p[i]->id);
      num_cals = array_list_size(sw_batch_p->allocate_cals_p[i]);
      //printf("RNA SERVER :: %d - %d\n", read_len, header_len);
      total += num_cals;
      //memcpy(read_seq, &(sw_batch_p->read_p[sw_batch_p->read_indices_p[read]]), read_len);
      
      memset(cals_per_chromosome, 0, CHROMOSOME_NUMBER * STRANDS_NUMBER *  sizeof(unsigned int));
      /*
      printf("\n");
      printf("CALS BEFORE ORDER:\n");
      for(int t = 0; t < array_list_size(sw_batch_p->allocate_cals_p[i]); t++){
	cal_p = (cal_t *)array_list_get(t, sw_batch_p->allocate_cals_p[i]);
	printf("\t%d.[chromosome:%d]-[strand:%d]-[start:%d, end:%d]\n", t, cal_p->chromosome_id, cal_p->strand, cal_p->start, cal_p->end);
      }*/
      
      
      // Order CALs and count number of CALs per chromosome and initialize cal type array
      for (j = 0; j < num_cals; j++) {
	/*cal_prev = (cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]);
	//printf("Iter %d\n", j);
	for (z = j + 1; z < num_cals; z++){
	  cal_next = (cal_t *)array_list_get(z, sw_batch_p->allocate_cals_p[i]);
	  //printf("\tCal prev:[%d-%d] <======> Cal next:[%d-%d]\n", cal_prev->start, cal_prev->end, cal_next->start, cal_next->end);
	  if ( cal_prev->chromosome_id > cal_next->chromosome_id ) {
	    //printf("\tSwap chromosome %d-%d\n", j, z);
	    array_list_swap(j, z, sw_batch_p->allocate_cals_p[i]);
	    cal_prev = (cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]);
	  }else if (cal_prev->chromosome_id == cal_next->chromosome_id){
	    if (cal_prev->strand > cal_next->strand){
	      //printf("\tSwap strand %d-%d\n", j, z);
	      array_list_swap(j, z, sw_batch_p->allocate_cals_p[i]);
	      cal_prev = (cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]);
	    }else if(cal_prev->strand == cal_next->strand){
	      if(cal_prev->start > cal_next->start){
		array_list_swap(j, z, sw_batch_p->allocate_cals_p[i]);
		cal_prev = (cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]);
	      }else if (cal_prev->start == cal_next->start) {
		if(cal_prev->end > cal_next->end){
		  //printf("\tSwap end %d-%d\n", j, z);
		  array_list_swap(j, z, sw_batch_p->allocate_cals_p[i]);
		  cal_prev = (cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]);
		}
	      }
	    }
	  }
	  }*/
	cal_p = (cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]);
	cals_per_chromosome[cal_p->chromosome_id * STRANDS_NUMBER + cal_p->strand] += 1;
      }
      /*
      cal_p = (cal_t *)array_list_get(j , sw_batch_p->allocate_cals_p[i]);	
      cals_per_chromosome[cal_p->chromosome_id * STRANDS_NUMBER + cal_p->strand] += 1;
      */
      /*printf("\n");
      printf("CALS AFTER ORDER:\n");
      for(int t = 0; t < array_list_size(sw_batch_p->allocate_cals_p[i]); t++){
	cal_p = (cal_t *)array_list_get(t, sw_batch_p->allocate_cals_p[i]);
	printf("\t%d.[chromosome:%d]-[strand:%d]-[start:%d, end:%d]\n", t, cal_p->chromosome_id, cal_p->strand, cal_p->start, cal_p->end);
	}*/
      

      //CALs actualization and extend 
      j = 0;
      cal_prev = (cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]);
      cal_prev->start -= flank_length;
      cal_prev->end += flank_length;
      j++;
      while(j < num_cals){	  
	cal_next = (cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]);
	cal_next->start -= flank_length;
	cal_next->end += flank_length;
	if((cal_next->chromosome_id == cal_prev->chromosome_id) && 
	   (cal_next->strand == cal_prev->strand) && 
	   (cal_next->start <= (cal_prev->end + seed_max_distance))){
	  //printf("Fusion CALs %d with %d\n", j - 1, j);
	  cal_next->start = cal_prev->start;
	  cal_prev = cal_next;
	  cal_p = array_list_remove_at(j - 1, sw_batch_p->allocate_cals_p[i]);
	  cals_per_chromosome[cal_p->chromosome_id * STRANDS_NUMBER + cal_p->strand] -= 1;
	  cal_free(cal_p);
	  num_cals--;
	}else{
	  cal_prev = cal_next;
	  j++;
	}
      }
      
      /*
      printf("\n");
      printf("CALS AFTER REFUSION:\n");
      for(int t = 0; t < array_list_size(sw_batch_p->allocate_cals_p[i]); t++){
	cal_p = (cal_t *)array_list_get(t, sw_batch_p->allocate_cals_p[i]);
	printf("\t%d.[chromosome:%d]-[strand:%d]-[start:%d, end:%d]\n", t, cal_p->chromosome_id, cal_p->strand, cal_p->start, cal_p->end);
	}	*/
      
      chromosome_tot = 0;
      
      for(int c = 0; c < CHROMOSOME_NUMBER*STRANDS_NUMBER; c+=2){
	if(cals_per_chromosome[c]){ chromosome_tot++; }
	if(cals_per_chromosome[c + 1]){ chromosome_tot++; }
      }
      
      //printf("cals chromosome tot %d\n", chromosome_tot);
      
      // Loop Chromosome CALs
      cal_chromosome_start = 0;
      cal_chromosome_end = 0;
      
      for(int c = 0; c < chromosome_tot; c++){
	cal_p = (cal_t *)array_list_get(cal_chromosome_start, sw_batch_p->allocate_cals_p[i]);
	//printf("Process read %s\n", sw_batch_p->allocate_reads_p[i]->id);
	strand = cal_p->strand;
	num_cals = cals_per_chromosome[cal_p->chromosome_id*STRANDS_NUMBER + cal_p->strand];
	//printf("Num Cals %d\n", num_cals);
	cal_chromosome_end += num_cals;	
	
	reference_fusion_p[0] = '\0';
	reference_fusion_len = 0;
	
	j = cal_chromosome_start;
	no_store = 0;
	//printf("Process from cal %d to cal %d\n", cal_chromosome_start, cal_chromosome_end);
	while( j < cal_chromosome_end ){
	  
	  cals_fusion_p = &allocate_cals_fusion_p[curr_depth];
	  cals_fusion_p->allocate_data = (cal_fusion_data_t *)calloc(num_cals, sizeof(cal_fusion_data_t));
	  
	  cal_p = (cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]);

	  start_before =  cal_p->start;
	  end_before =  cal_p->end;
	  len = (unsigned int)(end_before - start_before + 1);
	  
	  start = start_before;
	  end = end_before;

	  //printf("\tExtract cal %d : chr%d(%d)=[%d-%d]\n", j, cal_p->chromosome_id, strand, start, end);

	  num_cals_fusion = 0;
	  reference_fusion_len = 0;
	  
	  reference_fusion_p[0] = '\0';
	  //cal_indice = 0;
	  distance = 0;
	  
	  while(distance < max_intron_size){
	    
	    cal_fusion_data_init(j, start, end, cal_p->strand, cal_p->chromosome_id, 
				 reference_fusion_len, reference_fusion_len + len -1, &cals_fusion_p->allocate_data[num_cals_fusion]);
	    
	    if(len > max_reference){
	      max_reference += 1.25*len;
	      reference = (char *)malloc(max_reference*sizeof(char));
	    }
	    
	    genome_read_sequence_by_chr_index(reference, 0, 
					      cal_p->chromosome_id - 1, &start, &end, genome_p);
	    
	    reference_fusion_len += len;
	    
	    if(reference_fusion_len >= maximum_reference_len){
	      no_store = 1;
	      break;
	      //maximum_reference_len += 1.25*len;
	      //reference_fusion_p = (char *)realloc(reference_fusion_p, sizeof(char)*maximum_reference_len);
	    }	     
	    strcat(reference_fusion_p, reference);
	    j++;
	    num_cals_fusion++;
	    if( j < cal_chromosome_end){
	      //printf("\tFusion cals...\n");
	      cal_p = (cal_t *)array_list_get(j, sw_batch_p->allocate_cals_p[i]);
	      
	      start_before = start;
	      start = cal_p->start;
	      end = cal_p->end;
		
	      distance = start - start_before;
	      
	      len = (unsigned int)(end - start + 1);
	    }else{
	      break;
	    }
	  }//while to concatenate CALs 
	  
	  if(no_store){
	    no_store = 0;
	    break;
	  }
	  
	  //printf("\tProceses one sw\n");
	  cal_indice = 0;
	  cals_fusion_p->cal_number = num_cals_fusion;   	  
	  channel_p = &sw_channels_p[curr_depth];
	  
	  sw_channel_allocate_ref(reference_fusion_len + 1, channel_p);
	  
	  memcpy(channel_p->ref_p, reference_fusion_p, reference_fusion_len + 1);
	  channel_p->ref_p[reference_fusion_len]= '\0';
	  
	  sw_channel_update(i, strand, read_len, header_len, reference_fusion_len, channel_p);
	  
	  //Reverse read
	  //seq_p = (char *)malloc(sizeof(char *)*(read_len + 1));
	  seq_p = (char *)calloc((read_len + 1), sizeof(char));
	  memcpy(seq_p, sw_batch_p->allocate_reads_p[i]->sequence, read_len);
	  
	  if(strand == 1){
	    //printf("\tReverse complementary\n");
	    seq_reverse_complementary(seq_p, read_len);
	  }
	  seq_p[read_len] = '\0';

	  sw_simd_input_add(seq_p, read_len,
			    channel_p->ref_p, channel_p->ref_len, 
			    curr_depth, sw_input_p);
	  
	  if ((++curr_depth) == depth) {
	  
	    smith_waterman_simd(sw_input_p, sw_output_p, context_p);
	  
	    write_batch_p = search_splice_junctions_sw_output(sw_input_p, sw_output_p, curr_depth, 
							      allocate_cals_fusion_p, chromosome_avls_p,  mapping_reads_p, sw_channels_p, 
							      sw_batch_p, write_list_p, write_batch_p, write_size, sw_id, &sw_no_valids,
							      min_score, genome_p, min_intron_length);
	    num_sw_process += curr_depth;
	    curr_depth = 0;
	  }
	  
	}//while from cal_chromosome_start to cal_chromosome_end
	//printf("End process CALs\n");
	cal_chromosome_start += num_cals;
	
      }//end for start to chromosome_tot
      
    } // end of for 0..num_reads
    
    //printf("%d.Seq(%d): %s\n",cal_p->strand, read_len, sw_batch_p->allocate_reads_p[i]->sequence);
    if (curr_depth > 0) {
      //printf("Current depth %d/%d\n", curr_depth, depth);
      num_sw_process += curr_depth;
      //read = channel_p->read_index;
      //ref_p = channel_p->ref_p;
      //ref_len = channel_p->ref_len;
      //read_len = channel_p->read_len;
      //printf("Depth %d-%d\n", sw_input_p->depth, curr_depth);
      for(k = curr_depth ; k < depth ; k++) {
	sw_simd_input_add(seq_p, read_len,
			  channel_p->ref_p, channel_p->ref_len, 
			  k, sw_input_p);
	//printf("k=%d : %s\n", k , sw_input_p->seq_p[k]);
      }
      //printf("Run smith waterman...\n");
      smith_waterman_simd(sw_input_p, sw_output_p, context_p);
      write_batch_p = search_splice_junctions_sw_output(sw_input_p, sw_output_p, curr_depth, allocate_cals_fusion_p, chromosome_avls_p,  mapping_reads_p, sw_channels_p, sw_batch_p, write_list_p, write_batch_p, write_size, sw_id, &sw_no_valids, min_score, genome_p, min_intron_length);
      
      //found_write_p = process_sw_output(sw_output_p, sw_input_p, min_score, curr_depth, sw_channels_p, sw_batch_p, write_list_p, found_write_p, write_size, sw_id, &total_valids, mapping_reads_p, genome_p);
      curr_depth = 0;
    }
    
    for (i = 0; i < num_reads; i++) {
      if (!mapping_reads_p[i]) {
	//printf("****************** read %i NO MAPPED %i!!!\n", i, mapping_reads_p[i]);
	total_reads_unmapped++;
	read_len = strlen(sw_batch_p->allocate_reads_p[i]->sequence);
	header_len = strlen(sw_batch_p->allocate_reads_p[i]->id);
				
	if ( write_batch_p->size >= write_batch_p->allocated_size - 1) {
	  item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
	  if (time_on) { timing_stop(SW_SERVER, sw_id, timing_p); }
	  list_insert_item(item_p, write_list_p);
	  if (time_on) { timing_start(SW_SERVER, sw_id, timing_p); }
	  
	  write_batch_p = write_batch_new(write_size, MATCH_FLAG);
	}
	
	search_match_p = (char *)malloc(sizeof(char)*(read_len + 1));
	memcpy(search_match_p, sw_batch_p->allocate_reads_p[i]->sequence, read_len);
	search_match_p[read_len] = '\0';
	
	quality_match_p = (char *)malloc(sizeof(char)*(read_len + 1));
	memcpy(quality_match_p, sw_batch_p->allocate_reads_p[i]->quality, read_len);
	quality_match_p[read_len] = '\0';
	
	header_id_match_p = (char *)malloc(sizeof(char)*header_len);
	//memcpy(header_id_match_p, &header_id, len_id_header);
	len_id_header = get_to_first_blank(sw_batch_p->allocate_reads_p[i]->id, header_len, header_id_match_p);
				
	
	alignment_p = alignment_new();
	cigar_p = (char *)malloc(sizeof(char)*10);
	sprintf(cigar_p, "%d%c\0", read_len, 'X');
	//TODO:chromosome 0??
	
	alignment_init_single_end(header_id_match_p, search_match_p, quality_match_p, 0, 0, 0, cigar_p, 1, 255, 0, 0, alignment_p);
	
	//printf("seq: %s\n", alignment_p->sequence);
	((alignment_t **)write_batch_p->buffer_p)[write_batch_p->size] = alignment_p;
	write_batch_p->size++;
	
	mapping_reads_p[i] = 2;
      }else{
	//printf("****************** read %i MAPPED %i!!!\n", i, mapping_reads_p[i]);
      } // end of if mapping_reads_p[i] == 0
    } // end of for 0..num_reads
    
    /*    printf("Reads mapped with global variable at this moment %d\n", reads_mapped_global);
    printf("Reads mapped with local  variable at this moment %d\n", total_reads_unmapped);
    */
    sw_batch_free(sw_batch_p);
    list_item_free(sw_item_p);
    
    if (time_on) { timing_stop(SW_SERVER, sw_id, timing_p); }
    //printf("RNA Server process batch finish!\n");
    } // end of while
  
  // insert or free memory
  if (write_batch_p != NULL) {
    if (write_batch_p->size > 0) {
      item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
      list_insert_item(item_p, write_list_p);
    } else {
      write_batch_free(write_batch_p);
    }
  }
  
  for(k=0 ; k<depth ; k++) {
    free(sw_channels_p[k].ref_p);
  }
  
  free(sw_channels_p);
  sw_simd_input_free(sw_input_p);
  sw_simd_output_free(sw_output_p);
  sw_simd_context_free(context_p);
  
  free(allocate_cals_fusion_p);
  free(cals_per_chromosome);
  free(mapping_reads_p);
  list_decr_writers(write_list_p);
  free(reference_fusion_p);
  free(reference);

  if (statistics_on) { 
    statistics_add(SW_SERVER_ST, 0, num_batches, statistics_p); 
    statistics_add(SW_SERVER_ST, 1, total_reads_process, statistics_p); 
    statistics_add(SW_SERVER_ST, 2, total_reads_process - total_reads_unmapped, statistics_p); 
    statistics_add(SW_SERVER_ST, 3, total_reads_unmapped, statistics_p); 
    statistics_add(SW_SERVER_ST, 4, num_sw_process, statistics_p); 
    statistics_add(SW_SERVER_ST, 5, num_sw_process - sw_no_valids, statistics_p); 
    statistics_add(SW_SERVER_ST, 6, sw_no_valids, statistics_p); 
    
    statistics_add(TOTAL_ST, 1, total_reads_process - total_reads_unmapped, statistics_p); 
    statistics_add(TOTAL_ST, 2, total_reads_unmapped, statistics_p); 
  }


  printf("rna_server: END (%i reads -> unmapped %i, %i smith-waterman -> %i valids)\n", total_reads_process,  total_reads_unmapped, num_sw_process, total_reads_process - total_reads_unmapped);
  
  return;
  }
