#include <omp.h>
#include "rna_server.h"

#define MINIMUN_CAL_LENGTH 10
#define ERRORS_ZONE 8
#define MINIMUM_INTRON_LENGTH 10
#define MAX_CIGAR_OPERATIONS 20

#ifndef MAX
   #define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

//======== SPLICE JUNCTION TYPE ============
#define NOT_SPLICE	0
#define GT_AG_SPLICE  	1
#define CT_AC_SPLICE  	2
//==========================================

#define SEARCH_START 0
#define SEARCH_END 1
#define POSSIBLES_MARKS 2

#define STRANDS_NUMBER 2

int ii = -1;

void cal_fusion_data_init(unsigned int id, size_t start, size_t end, 
			  unsigned int strand, unsigned int chromosome, unsigned int fusion_start, 
			  unsigned int fusion_end, cal_fusion_data_t *cal_fusion_data_p) {
  cal_fusion_data_p->genome_strand = strand;
  cal_fusion_data_p->id = id;
  cal_fusion_data_p->genome_start = start;
  cal_fusion_data_p->genome_chromosome = chromosome;
  cal_fusion_data_p->genome_end = end;
  cal_fusion_data_p->fusion_start = fusion_start;
  cal_fusion_data_p->fusion_end = fusion_end;
}

char* cigar_automata_status(unsigned char status){
  char *str_status = (char *)malloc(sizeof(char)*100);  
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

int cigar_generator(cigar_data_t *cigar_p, unsigned int *max_len, 
		    unsigned int *pos, unsigned int status, int *number, int *distance){
  char cigar_info[20];
  //printf("Store %d in %i position\n", *number, *pos);
  /*if (*pos > MAX_CIGAR_OPERATIONS) {
    return 2;
    }*/
  
  if (*number <= 0) {
    return 1;
  }else{
    //printf("**********Allocate %d\n",*number );
    if ((status != CIGAR_DELETION) || ((status == CIGAR_DELETION) && (*number <= MINIMUM_INTRON_LENGTH))){
	  cigar_p[*pos].type = status;
	  cigar_p[*pos].value = *number;	  
	  *pos += 1;
	  if (status != CIGAR_PADDING && status != CIGAR_MATCH_MISMATCH) {
	    *distance += *number;
	  }
	  *number = 1;
	  if (*pos >= *max_len){
	    *max_len = *max_len * 2;
	    cigar_p = (cigar_data_t *)realloc(cigar_p, sizeof(cigar_data_t)*(*max_len));
	  }
      }
  }
  return 0;
}

//list_t* write_list_p, write_batch_t* write_batch_p, unsigned int write_size,

void search_splice_junctions_sw_output(sw_simd_input_t* input_p, sw_simd_output_t* output_p, 
				       unsigned int depth, allocate_fusion_data_t *depth_cal_fusion_p, 
				       allocate_splice_elements_t *chromosome_avls_p,  
				       sw_channel_t *sw_channels_p, mapping_batch_t *mapping_batch_p,  unsigned int sw_id,
				       size_t *sw_no_valids, float min_score, genome_t *genome_p, int min_intron_length){


  //printf("Process Splice Junctions\n");
  const unsigned char EXTRA_SEARCH = 8;
  char header_id[2048];
  alignment_t *alignment_p;
  
  char *header_match, *read_match, *quality_match, *quality_strand;
  
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
  size_t str_pos;
  unsigned int quality_pos;
  unsigned int reference_position_disp;
  unsigned int start_gap_padding;
  unsigned int start_search_j;

  unsigned int extend_start_sp, extend_end_sp;

  unsigned int splice_number;
  size_t maximum_splice_number = 100;
  sp_data_t *allocate_splice = (sp_data_t *)malloc(sizeof(sp_data_t)*maximum_splice_number);//[maximum_splice_number];

  //=========CIGAR INFO. && BAM/SAM ========//
  int error_cigar;
  unsigned int cigar_soft;
  unsigned int value;
  int cigar_value;
  unsigned int cigar_max_len = 2048;
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
  //unsigned int change_cal_index[20];
  //list_item_t* item_p = NULL;  
  char *cigar_str;
  cigar_data_t *cigar_p = (cigar_data_t *)malloc(sizeof(cigar_data_t)*cigar_max_len);
  char *optional_fields, *p;
  int optional_fields_length;
  int distance;
  int AS;
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

  /*    
  printf("======================== Process Output SW %d=========================\n", depth);
  sw_simd_input_display(depth, input_p);
  sw_simd_output_display(depth, output_p);
  printf("======================================================================\n");
  */

  for (int i = 0; i < depth; i++) {    
    splice_number = 0;
    total_cals = depth_cal_fusion_p[i].cal_number;

    //printf("Total cals %d\n", total_cals);
    //if (total_cals == 1) {
    if (output_p->norm_score_p[i] < min_score) {
      *sw_no_valids += 1;
      continue;
    }
    AS =(int)output_p->score_p[i];
        
    /*printf("=========================== BEFORE [Depth %d - Total CALs %d]==================================\n",
    i, depth_cal_fusion_p[i].cal_number);
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
    fastq_read_t *read_p = array_list_get(sw_channels_p[i].read_index, mapping_batch_p->fq_batch);
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
	//printf("ERROR: total cals overflow\n");
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
	//printf("ERROR: total cals overflow\n");
	continue;
      }
    }

    depth_cal_fusion_p[i].allocate_data[cal].genome_end -= (depth_cal_fusion_p[i].allocate_data[cal].fusion_end - output_p->mapped_len_p[i]);
    depth_cal_fusion_p[i].allocate_data[cal].fusion_end  = output_p->mapped_len_p[i] - 1;
        
    j_start = 0;
    while (((depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start) < 
	   MINIMUN_CAL_LENGTH) && (actual_cal < total_cals)) {
      j_start += (depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start) + 1;
      depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start = 0;
      depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end = 0;
      actual_cal++;
      if (actual_cal >= total_cals) {
	//printf("ERROR: total cals overflow\n");
	continue;
      }
      start_mapped = depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start;      
      //printf("Move to next CAL %i, start_j = %i\n", actual_cal, j_start);
    }

    while (output_p->mapped_ref_p[i][j_start] != output_p->mapped_seq_p[i][j_start]) {
      j_start++;
      if (j_start >= output_p->mapped_len_p[i]){ continue; }
    }

    cal = actual_cal;
    change_cal_number = 0;
    real_sw_len = 0;

    //memset( change_cal_index, 0, sizeof(unsigned int)*total_cals);    
    //printf("Actual CAL %d", actual_cal);
    //printf("=========================== AFTER [Depth %d - Total CALs %d]==================================\n", i, depth_cal_fusion_p[i].cal_number);
    //for(int c = 0; c < depth_cal_fusion_p[i].cal_number; c++){
      //printf("CAL %d :\n", c);
      //printf("\t->Genome Start : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_start);
      //printf("\t->Genome End : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_end);
      //printf("\t->Fusion Start : %d\n", depth_cal_fusion_p[i].allocate_data[c].fusion_start);
      //printf("\t->Fusion End : %d\n", depth_cal_fusion_p[i].allocate_data[c].fusion_end);
      //printf("\t->Chromosome : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_chromosome);
      //printf("\t->Strand : %d\n", depth_cal_fusion_p[i].allocate_data[c].genome_strand);
    //}
    //printf("=============================================================================================================\n");
    
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
    distance = 0;
    optional_fields = (char *)malloc(sizeof(char)*100);

    read_len = output_p->mapped_len_p[i];
    read_match = (char *)malloc(sizeof(char)*(read_len + 1));
    if (read_match == NULL) { printf("Error to allocate memory\n"); exit(-1); }
   
    quality_pos = output_p->start_seq_p[i];
    quality_match = (char *)malloc(sizeof(char)*(read_len + 1));
    if (quality_match == NULL) { printf("Error to allocate memory\n"); exit(-1); }
    
    found = NOT_SPLICE;
    //===========Position Automata in initial State==========
    cigar_value = 0;
    cigar_pos = 0;
    
    if (output_p->mapped_ref_p[i][j_start] != '-') {
      if (output_p->mapped_seq_p[i][j_start] != '-') {
	automata_status = CIGAR_MATCH_MISMATCH;
      } else {
	automata_status = CIGAR_DELETION;
      }
    } else {
      automata_status = CIGAR_INSERTION;
    }
    //========================================================
    
    change_cal = 0;
    //printf("\tStart Search\n");
    for (int j = j_start; j < output_p->mapped_len_p[i]; j++) {
      if (output_p->mapped_ref_p[i][j] != '-') {
	reference_position_disp++;
	start_gap_padding++;
	relative_gap_start++;
	//TODO: IMPLEMENTS WHEN JUMP TO OTHER CALL WITHOUT GAP
	if (output_p->mapped_seq_p[i][j] != '-') {
	  //Match mismatch area	  
	  if (output_p->mapped_seq_p[i][j] != output_p->mapped_ref_p[i][j]) { distance++; /*printf("Dist\n");*/ } 
	  if (automata_status == CIGAR_MATCH_MISMATCH) {
	    cigar_value++;
	  } else {
	    //printf("%s: %d in [j=%d - jstar=%d]\n", cigar_automata_status(automata_status), cigar_value, j, j_start);
	    error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value, &distance);
	    if (error_cigar) { end_exceded = 1; break;}//printf("MATCH-MISMATCH :: %s\n", read_p->id);printf("Error cigar\n"); }
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
	  //printf("%s: %d in [j=%d - jstar=%d]\n", cigar_automata_status(automata_status), cigar_value, j, j_start);
	  //distance += cigar_value;
	  error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value, &distance);
	  if (error_cigar) { end_exceded = 1; break;}//if (error_cigar) { printf("INSERTION :: %s\n", read_p->id);printf("Error cigar"); }
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
	    error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value, &distance);
	     if (error_cigar) { end_exceded = 1; break;}//if (error_cigar) { printf("%s\n", read_p->id);printf("Error cigar"); }
	    automata_status = CIGAR_PADDING;    
	  }
	} else {
	  //Deletion area
	  if (automata_status == CIGAR_DELETION) {
	    cigar_value++;
	  } else {
	    //printf("%s: %d in [j=%d - jstar=%d]\n", cigar_automata_status(automata_status), cigar_value, j, j_start);
	    //distance += cigar_value;
	    error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value, &distance);
	    if (error_cigar) { end_exceded = 1; break;}//if (error_cigar) { printf("DELETION :: %s\n", read_p->id);printf("Error cigar"); }
	    automata_status = CIGAR_DELETION;
	  }
	}
      } else { 
	 read_match[str_pos] = output_p->mapped_seq_p[i][j];
	 quality_match[str_pos++] = read_p->quality[quality_pos];//sw_batch_p->allocate_reads_p[sw_channels_p[i].read_index]->quality[quality_pos];
	 quality_pos++;
	 
	 if (n_deletions_seq > MINIMUM_INTRON_LENGTH) {
	    
	   gap_end = j - 1;   
	   //printf("gap_start=%d , gap_end=%d\n", gap_start, gap_end);

	   //Canonical marks strand + (GT-AG)
	   if ((output_p->mapped_ref_p[i][gap_start] == canonical_GT_AG[0]) && 
	       (output_p->mapped_ref_p[i][gap_start + 1] == canonical_GT_AG[1])) {
	     if ((output_p->mapped_ref_p[i][gap_end - 1] == canonical_GT_AG[2]) && 
		 (output_p->mapped_ref_p[i][gap_end] == canonical_GT_AG[3])) {
	       found = GT_AG_SPLICE; 
	     }
	   }
	   //Canonical marks strand - (CT-AC)
	   if ((output_p->mapped_ref_p[i][gap_start] == canonical_CT_AC[0]) && 
	       (output_p->mapped_ref_p[i][gap_start + 1] == canonical_CT_AC[1])) {
	     if ((output_p->mapped_ref_p[i][gap_end - 1] == canonical_CT_AC[2]) && 
		 (output_p->mapped_ref_p[i][gap_end] == canonical_CT_AC[3])) {
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

	       if ((output_p->mapped_ref_p[i][position] == canonical_GT_AG[0]) && 
		   (output_p->mapped_ref_p[i][position + 1] == canonical_GT_AG[1])) {
		 //printf("-Found GT-");
		 if ((output_p->mapped_ref_p[i][gap_end - 1 + count_extra_search] == canonical_GT_AG[2]) && 
		     (output_p->mapped_ref_p[i][gap_end + count_extra_search] == canonical_GT_AG[3])) {
		   confirmation_type = GT_AG_SPLICE; 
		   //printf("Found GT_AG Extra Search\n");
		 }
		 
	       } else if ((output_p->mapped_ref_p[i][position] == canonical_CT_AC[0]) && 
			  (output_p->mapped_ref_p[i][position + 1] == canonical_CT_AC[1])) {
		 if ((output_p->mapped_ref_p[i][gap_end - 1 + count_extra_search] == canonical_CT_AC[2]) && 
		     (output_p->mapped_ref_p[i][gap_end + count_extra_search] == canonical_CT_AC[3])) {
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
	      tmp = ((depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - j_start) - 
		     (depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start - j_start) + 1) - relative_gap_start ;
	      
	      relative_gap_start = n_deletions_seq - tmp;
	     
	      actual_cal++;
	      nt_skipped = 0;

	      if (gap_end >  depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
		while (gap_end >  depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
		  nt_skipped += (depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end - 
				 depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_start);
		  actual_cal++;
		}
		nt_skipped++;
	      }
	      
	      if (actual_cal >= total_cals){ end_exceded = 1; break; }
	      
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
	      if (splice_number >= maximum_splice_number) {
		maximum_splice_number *= 1.25;
		allocate_splice = (sp_data_t *)realloc(allocate_splice, sizeof(sp_data_t)*maximum_splice_number);
	      }
	      cigar_p[cigar_pos - 1].value += mark_dpl;	    
	      error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, CIGAR_SKIPPED, &cigar_value, &distance);
	      if (error_cigar) { end_exceded = 1; break;}
	      //printf("::: %d || %d :::\n", cigar_value, mark_dpl);	      
	      cigar_value -= mark_dpl;
	      //printf("::: %d || %d :::\n", cigar_value, mark_dpl);
	    } else {
	      //distance += cigar_value;
	      error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, CIGAR_BIG_DELETION, &cigar_value, &distance);
	       if (error_cigar) { end_exceded = 1; break;}
	    }
	    
	    found = NOT_SPLICE;
	    deletion_number = 0;
	    //printf("==================== Exact->( %d:%d-%d ) Extend->( %d:%d-%d )  ========================\n", depth_cal_fusion_p[i].allocate_data->genome_chromosome, start_splice, end_splice, depth_cal_fusion_p[i].allocate_data->genome_chromosome, extend_start_sp, extend_end_sp);
	  } else { //IF not found
	    if (reference_position_disp > depth_cal_fusion_p[i].allocate_data[actual_cal].fusion_end) {
	      actual_cal++;
	      //printf("SPLICE NOT FOUND report deletion\n");
	      if (actual_cal >= total_cals) {
		//Cals overflow, read no valid		
		end_exceded = 1;
		break;
	      } else {
		
		relative_gap_start = gap_end   - (depth_cal_fusion_p[i].allocate_data[actual_cal - 1].fusion_end + insertion_number);
		
		cigar_value = (depth_cal_fusion_p[i].allocate_data[actual_cal].genome_start + relative_gap_start) -  
		  (depth_cal_fusion_p[i].allocate_data[actual_cal - 1].genome_end - ((gap_end - gap_start) - relative_gap_start));
	      }
	      //printf("%s: %d in [j=%d - jstar=%d]\n", cigar_automata_status(CIGAR_BIG_DELETION), cigar_value, j, j_start);
	      //distance += cigar_value;
	      error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos,  CIGAR_BIG_DELETION, &cigar_value, &distance);
	      if (error_cigar) { end_exceded = 1; break;}
	    } else {
	     cigar_value = gap_end - gap_start + 1;
	     //printf("%s: %d in [j=%d - jstar=%d]\n", cigar_automata_status(CIGAR_BIG_DELETION), cigar_value, j, j_start);
	     //distance += cigar_value;
	     error_cigar = cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, CIGAR_BIG_DELETION, &cigar_value, &distance);
	     if (error_cigar) { end_exceded = 1; break;}
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

    //printf("\tEnd process sw %i\n", i);
    //printf("\tStart Search\n");
    if (end_exceded) { free(optional_fields); free(read_match); free(quality_match); end_exceded = 0; continue; }
    //printf("\tEnd Search Splice\n");
    //=============================== MAKE CIGAR & SAM/BAM =====================================
    cigar_generator(cigar_p, &cigar_max_len, &cigar_pos, automata_status, &cigar_value, &distance);

    if (cigar_pos < 1) { free(optional_fields); free(read_match); free(quality_match); continue; }

    //Generate cigar string
    num_cigar_op = 0;
    cigar_str = (char *)malloc(sizeof(char)*cigar_max_len);
    if (cigar_str == NULL) { exit(-1); }
    cigar_str[0] = '\0';
    
    //Hard and Soft clipped start
    if (output_p->start_seq_p[i] > 0) {
      sprintf(cigar_segment, "%iH", output_p->start_seq_p[i] + j_start);
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
		//printf("Soft (%c!=%c)", output_p->mapped_ref_p[i][cigar_soft - 1], output_p->mapped_seq_p[i][cigar_soft - 1]);
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
      //printf("%u\n", cigar_p[cig].value);
      sprintf(cigar_segment, "%u%c", cigar_p[cig].value, flag);
      //printf("%s\n", cigar_segment);
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
      sprintf(cigar_segment, "%iH", input_p->seq_len_p[i] - 
	      ((output_p->mapped_len_p[i] - deletions_tot) + output_p->start_seq_p[i]));
      cigar_str = strcat(cigar_str, cigar_segment);
      num_cigar_op++;
    }
    
    //printf("Cigar(%d):%s\n", cigar_pos, cigar_str);
    //end cigar string
//    if (array_list_size(mapping_batch_p->mapping_lists[sw_channels_p[i].read_index])) {
    if (array_list_size(mapping_batch_p->mapping_lists[sw_channels_p[i].read_index]) > 1) {
      primary_alignment = 1;
    } else {
      primary_alignment = 0;
    }

    alignment_p = alignment_new();
    header_len = sw_channels_p[i].header_len;
    //printf("Reallocate %d\n", header_len); 
    get_to_first_blank(read_p->id, header_len, header_id);
    header_match = (char *)malloc(sizeof(char)*header_len);
    if (header_match == NULL) { exit(-1); }
    memcpy(header_match, header_id, header_len);
    
    //printf("Header len = %d\n", header_len);

    read_match[str_pos] = '\0';
    //quality_match[str_pos] = '\0';
    //memcpy(quality_match, &(sw_batch_p->quality_p[sw_batch_p->read_indices_p[sw_channels_p[i].read_index]]), read_len );
    if (depth_cal_fusion_p[i].allocate_data->genome_strand) {
      quality_strand = (char *)malloc(sizeof(char)*(str_pos + 1));
      reverse_str(quality_match, quality_strand, str_pos);
      //memcpy(quality_strand, quality_match, str_pos);
      free(quality_match);
      quality_strand[str_pos] = '\0';
      quality_match = quality_strand;
    }
    //printf("CIGAR: %s\n", cigar_str);
    p = optional_fields;
    sprintf(p, "ASi");
    p += 3;
    memcpy(p, &AS, sizeof(int));
    p += sizeof(int);

    sprintf(p, "NMi");
    p += 3;
    memcpy(p, &distance, sizeof(int));
    p += sizeof(int);

    optional_fields_length = p - optional_fields;
    //printf("Start mapped %i\n", start_mapped);
    alignment_init_single_end(header_match, read_match, quality_match, depth_cal_fusion_p[i].allocate_data->genome_strand, 
			      depth_cal_fusion_p[i].allocate_data->genome_chromosome - 1, start_mapped - 1, cigar_str, num_cigar_op, 
			      (int)output_p->norm_score_p[i] * 254, 1, primary_alignment, optional_fields_length, optional_fields, alignment_p);

    //printf("Score %i\n", alignment_p->map_quality);
    //printf("seq: %s\n", alignment_p->sequence);
    //printf("Insert mapping\n");
    

    //--------->
    //printf("-->Allocate mapping in %i position\n", sw_channels_p[i].read_index);
    array_list_insert(alignment_p, mapping_batch_p->mapping_lists[sw_channels_p[i].read_index]);
    //if (ii == sw_channels_p[i].read_index){ printf("SW READ list size %i - read position %i\n", array_list_size(allocate_mappings), sw_channels_p[i].read_index); 
    //    alignment_print(alignment_p);
    

    //============================================================================================
    
    //printf(" (score = %.2f, norm. score = %.2f)\n", output_p->score_p[i], output_p->norm_score_p[i]);
   
    //printf("Free done 1\n");
    
    //mapping_reads_p[sw_channels_p[i].read_index] = 1;
    
    //Report Splice Junctions
    for (unsigned int s = 0; s < splice_number; s++) {
      pthread_mutex_lock(&(chromosome_avls_p[depth_cal_fusion_p[i].allocate_data[0].genome_chromosome].mutex));
          allocate_new_splice(depth_cal_fusion_p[i].allocate_data->genome_chromosome, allocate_splice[s].strand_sp, 
			  allocate_splice[s].end_sp, allocate_splice[s].start_sp, allocate_splice[s].start_extend_sp, 
			  allocate_splice[s].end_extend_sp, chromosome_avls_p);
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
  free(allocate_splice);

  //printf("Process Splice Junctions END\n");
  //return write_batch_p;
  
}



void rna_server_omp_smith_waterman(sw_server_input_t* input_p, allocate_splice_elements_t *chromosome_avls_p){
  /**----------------------------------------------------------------------------------**
   **                                    RNA-SEQ                                       **
   **----------------------------------------------------------------------------------**/
  
  LOG_DEBUG_F("rna_server (%i): START\n", omp_get_thread_num());
  unsigned int seed_max_distance = input_p->seed_max_distance;

  size_t j, z, k, num_reads, read, bytes, num_cals, cals_positive, cals_negative;
  size_t i;

  //lists and items
  list_t* sw_list_p = input_p->sw_list_p;
  list_item_t *sw_item_p = NULL;
  mapping_batch_t* mapping_batch_p = NULL;
  
  list_item_t* item_p = NULL;
	
  list_t* alignment_list_p = input_p->alignment_list_p;
  unsigned int write_size = input_p->write_size;
  //write_batch_t* write_batch_p = write_batch_new(write_size, MATCH_FLAG);
  
  int sw_id = omp_get_thread_num();

  // genome
  char* ref_p;
  unsigned int ref_len;
  
  genome_t* genome_p = input_p->genome_p;
  int min_intron_length = input_p->min_intron_size;
  
  // SIMD support for Smith-Waterman
  unsigned int depth = 4, curr_depth = 0;
  sw_simd_input_t* sw_input_p = sw_simd_input_new(depth);
  sw_simd_output_t* sw_output_p = sw_simd_output_new(depth);
  
  sw_simd_context_t *context_p = sw_simd_context_new(input_p->match, input_p->mismatch, input_p->gap_open, input_p->gap_extend);
  

  float min_score = input_p->min_score;
	
  // for tracking what reads are being mapped successfully
  size_t allocated_cals = 10000;
  size_t maximum_allocate = 100;
  array_list_t *cals_list =  array_list_new(MAX_RNA_CALS + 32,
					    1.25f,
					    COLLECTION_MODE_ASYNCHRONIZED);
 //= (unsigned char*) calloc (allocated_mapping_reads, sizeof(unsigned char));
  //if (mapping_reads_p == NULL) { exit(-1); }

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
  size_t start, end, start_before, end_before;
  size_t len;
  char *reference, header_id[1024];
  unsigned int len_increment = 512;
  //unsigned int pos_strand;
  char action;
  unsigned int reference_fusion_len = 0;

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
  unsigned int chr_before;
  unsigned int num_cals_fusion;
  //  int distance;
  char no_store;
  size_t distance;
  unsigned int max_intron_size = input_p->max_intron_size;
  
  int num_mappings;
  alignment_t *alignment_aux;
  array_list_t **allocate_mappings;

  size_t n_alignments_report = 0; 
  size_t n_total_reads_process = 0;
  char *read_reverse;
  unsigned char found_negative;
  size_t reads_no_random = 0;
  size_t num_targets;
  size_t max_reference = 2048;
  size_t maximum_reference_len = 1024;//2048;
  char *reference_fusion_p = (char *)malloc(sizeof(char)*maximum_reference_len);
  if (reference_fusion_p == NULL) { exit(-1); }
  char *optional_fields;
  //int optional_fields_length;

  //allocate_mappings = (array_list_t **)malloc(allocated_mapping_reads*sizeof(array_list_t *));

  reference = (char *)malloc(sizeof(char)*max_reference);
  
  if (reference == NULL) { exit(-1); }

  // start main loop
  while ( (sw_item_p = list_remove_item(sw_list_p)) != NULL ) {
    num_batches++;
    //printf("RNA Server Processing Batch...\n");
    if (time_on) { timing_start(SW_SERVER, sw_id, timing_p); }
    
    curr_depth = 0;
    mapping_batch_p = (mapping_batch_t*) sw_item_p->data_p;
    num_targets = mapping_batch_p->num_targets;//sw_batch_p->num_reads;
    num_reads = array_list_size(mapping_batch_p->fq_batch);
    total_reads += num_targets;
    
    /*    
    if (num_reads >= allocated_mapping_reads) {
      cals_list = num_reads*1.25;
      //mapping_reads_p = (unsigned char *) realloc(mapping_reads_p, allocated_mapping_reads * sizeof(unsigned char));
      allocate_mappings = (array_list_t **)realloc(allocate_mappings, allocated_mapping_reads *sizeof(array_list_t *));
    }
    */
    //memset(mapping_reads_p, 0, allocated_mapping_reads * sizeof(unsigned char));

    /** for each read **/
    //total_reads_process += num_reads;
    for (i = 0; i < num_targets; i++) {      
      //if (allocate_mappings[i] == NULL) {
      /*allocate_mappings[i] = array_list_new(100, 
					    1.25f, 
      				    COLLECTION_MODE_ASYNCHRONIZED);
	//      }
	*/
      fastq_read_t *read_p = array_list_get(mapping_batch_p->targets[i], mapping_batch_p->fq_batch); 
      read_len = strlen(read_p->sequence);
      //printf("Proces(%d) %s\n", read_len, read_p->id);
      header_len = strlen(read_p->id);
      num_cals = array_list_size(mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
      //printf("Process %lu CALs\n", num_cals);
      //printf("RNA SERVER :: %d - %d\n", read_len, header_len);
      total += num_cals;
      //memcpy(read_seq, &(sw_batch_p->read_p[sw_batch_p->read_indices_p[read]]), read_len);
      
      memset(cals_per_chromosome, 0, CHROMOSOME_NUMBER * STRANDS_NUMBER *  sizeof(unsigned int));
      /*   
      printf("\n");
      printf("CALS BEFORE ORDER:\n");
      for(int t = 0; t < array_list_size(mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]); t++){
	cal_p = (cal_t *)array_list_get(t, mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
	printf("\t%d.[chromosome:%d]-[strand:%d]-[start:%d, end:%d]\n", t, cal_p->chromosome_id, cal_p->strand, cal_p->start, cal_p->end);
      }*/
      //exit(-1);
      
      // Order CALs and count number of CALs per chromosome and initialize cal type array
      found_negative = 0;
      for (j = 0; j < num_cals; j++) {
	cal_p = (cal_t *)array_list_get(j, mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
	if (cal_p->strand == 1 && !found_negative) {
	  read_reverse = (char *)calloc((read_len + 1), sizeof(char));
	  memcpy(read_reverse, read_p->sequence, read_len);
	  seq_reverse_complementary(read_reverse, read_len);
	  found_negative = 1;
	}
	cals_per_chromosome[cal_p->chromosome_id * STRANDS_NUMBER + cal_p->strand] += 1;
	array_list_insert(cal_p, cals_list);
      }
      mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]->size = 0;

      //CALs actualization and extend 
      j = 0;
      cal_prev = (cal_t *)array_list_get(j, cals_list);
      
      if (cal_prev->start <= flank_length) {
        cal_prev->start = 0;
      }else {
        cal_prev->start -= flank_length;
      }

      cal_prev->end += flank_length;
      if (cal_prev->end >= genome_p->chr_size[cal_prev->chromosome_id - 1]) {
	cal_prev->end = genome_p->chr_size[cal_prev->chromosome_id - 1] - 1;
      }

      j++;
      while(j < num_cals){	  
	cal_next = (cal_t *)array_list_get(j, cals_list);

	if (cal_next->start <= flank_length) {
          cal_next->start = 0;
        }else {
          cal_next->start -= flank_length;
        }

        cal_next->end += flank_length;
        if (cal_next->end >= genome_p->chr_size[cal_next->chromosome_id - 1]) {
          cal_next->end = genome_p->chr_size[cal_next->chromosome_id - 1] - 1;
        }

	if((cal_next->chromosome_id == cal_prev->chromosome_id) && 
	   (cal_next->strand == cal_prev->strand) && 
	   (cal_next->start <= (cal_prev->end + seed_max_distance))){
	  //printf("Fusion CALs %d with %d\n", j - 1, j);
	  cal_next->start = cal_prev->start;
	  cal_prev = cal_next;
	  cal_p = array_list_remove_at(j - 1, cals_list);
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
	}	
      */
      /*
      printf("\n");
      printf("CALS BEFORE ORDER:\n");
      for(int t = 0; t < array_list_size(cals_list); t++){
	cal_p = (cal_t *)array_list_get(t, cals_list);
	printf("\t%d.[chromosome:%d]-[strand:%d]-[start:%d, end:%d]\n", t, cal_p->chromosome_id, cal_p->strand, cal_p->start, cal_p->end);
	}*/
      //if (array_list_size(sw_batch_p->allocate_cals_p[i]) > 10) { printf("%s\n", sw_batch_p->allocate_reads_p[i]->sequence);exit(0);}
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
	cal_p = (cal_t *)array_list_get(cal_chromosome_start, cals_list);
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
	  
	  cal_p = (cal_t *)array_list_get(j, cals_list);

	  start_before =  cal_p->start;
	  end_before =  cal_p->end;
	  len = end_before - start_before + 1;
	  
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
				 reference_fusion_len, reference_fusion_len + len -1, 
				 &cals_fusion_p->allocate_data[num_cals_fusion]);
	    
	    if(len > max_reference){
	      max_reference += 1.25*len;
	      free(reference);
	      reference = (char *)malloc(max_reference*sizeof(char));
	      if (reference == NULL) { exit(-1); }
	    }
	    
	    genome_read_sequence_by_chr_index(reference, 0, 
					      cal_p->chromosome_id - 1, &start, &end, genome_p);
	    
	    //printf("chr:%i - CAL[%i-%i]\n", cal_p->chromosome_id, start, end);
	    reference_fusion_len += len;
	    
	    if(reference_fusion_len >= maximum_reference_len){
	      no_store = 1;
	      //printf("No process\n");
	      break;
	      // printf("Reallocate\n");
	      //maximum_reference_len += 1.25*len;
	      //reference_fusion_p = (char *)realloc(reference_fusion_p, sizeof(char)*maximum_reference_len);
	    }	     
	    strcat(reference_fusion_p, reference);
	    j++;
	    num_cals_fusion++;
	    if( j < cal_chromosome_end){
	      //printf("\tFusion cals...\n");
	      cal_p = (cal_t *)array_list_get(j, cals_list);
	      
	      start_before = start;
	      start = cal_p->start;
	      end = cal_p->end;
		
	      distance = start - start_before;
	      
	      len = end - start + 1;
	    }else{
	      break;
	    }
	  }//while to concatenate CALs 
	  
	  if(no_store){
	    free(cals_fusion_p->allocate_data);
	    no_store = 0;
	    break;
	  }
	  

	  //printf("Proceses one sw %i\n", reference_fusion_len);
	  
	  cal_indice = 0;
	  cals_fusion_p->cal_number = num_cals_fusion;   	  
	  channel_p = &sw_channels_p[curr_depth];
	  
	  sw_channel_allocate_ref(reference_fusion_len + 1, channel_p);
	  
	  memcpy(channel_p->ref_p, reference_fusion_p, reference_fusion_len + 1);
	  channel_p->ref_p[reference_fusion_len]= '\0';
	  
	  sw_channel_update(mapping_batch_p->targets[i], strand, read_len, header_len, reference_fusion_len, channel_p);
	  //i
	  //Reverse read
	  //seq_p = (char *)malloc(sizeof(char *)*(read_len + 1));
	  seq_p = (char *)calloc((read_len + 1), sizeof(char));
	  	  
	  if(strand == 1){
	    //printf("\tReverse complementary\n");
	    memcpy(seq_p, read_reverse, read_len);
	  }else {
	    memcpy(seq_p, read_p->sequence, read_len);
	  }

	  //seq_p[read_len] = '\0';

	  sw_simd_input_add(seq_p, read_len,
			    channel_p->ref_p, channel_p->ref_len, 
			    curr_depth, sw_input_p);

	  //	  printf("Stored int position %i\n", i); 

	  
	  if ((++curr_depth) == depth) {
	    smith_waterman_simd(sw_input_p, sw_output_p, context_p);
	    search_splice_junctions_sw_output(sw_input_p, sw_output_p, curr_depth, 
					      allocate_cals_fusion_p, chromosome_avls_p, sw_channels_p, 
					      mapping_batch_p, sw_id, &sw_no_valids,
					      min_score, genome_p, min_intron_length);
	    num_sw_process += curr_depth;
	    curr_depth = 0;
	  }
	  
	  //printf("Process End\n");
	}//while from cal_chromosome_start to cal_chromosome_end
	//printf("End process CALs\n");
	cal_chromosome_start += num_cals;
	
      }//end for start to chromosome_tot
      //printf("End process read %s\n", read_p->id);
      if (found_negative) { free(read_reverse); } 
      array_list_clear(cals_list, cal_free);
    } // end of for 0..num_reads
    
    //printf("%d.Seq(%d): %s\n",cal_p->strand, read_len, sw_batch_p->allocate_reads_p[i]->sequence);
    if (curr_depth > 0) {
      //printf("Current depth %d/%d\n", curr_depth, depth);
      num_sw_process += curr_depth;
      for(k = curr_depth ; k < depth ; k++) {
	//printf("STORE : %s\n", sw_input_p->seq_p[0]);
	sw_simd_input_add(sw_input_p->seq_p[0], sw_input_p->seq_len_p[0],
			  channel_p->ref_p, channel_p->ref_len, 
			  k, sw_input_p);
	//printf("k=%d : %s\n", k , sw_input_p->seq_p[k]);
      }
      //printf("Run smith waterman... %i\n", i);
      smith_waterman_simd(sw_input_p, sw_output_p, context_p);
      
      search_splice_junctions_sw_output(sw_input_p, sw_output_p, curr_depth, 
					allocate_cals_fusion_p, chromosome_avls_p,  
					sw_channels_p, mapping_batch_p,
					sw_id, &sw_no_valids, min_score, genome_p, 
					min_intron_length);
      
      //found_write_p = process_sw_output(sw_output_p, sw_input_p, min_score, curr_depth, sw_channels_p, sw_batch_p, write_list_p, found_write_p, write_size, sw_id, &total_valids, mapping_reads_p, genome_p);
      curr_depth = 0;
    }
    
    //**************************************************
    //PROCESS OUTPUT ALIGNMENTS
    for (i = 0; i < num_targets; i++) {
      //num_mappings = array_list_size(mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
      num_mappings = alignments_filter(input_p->bwt_optarg_p->report_all, 
				       input_p->bwt_optarg_p->report_best, 
				       input_p->bwt_optarg_p->report_n_hits, 
				       mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]); 

      if ( num_mappings > 0 ) {
	array_list_set_flag(1, mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
	for (int m = 0; m < num_mappings; m++) {
	  alignment_p = array_list_get(m, mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
	  optional_fields = (char *)alignment_p->optional_fields;
	  optional_fields += alignment_p->optional_fields_length;
	  sprintf(optional_fields, "NHi");
	  alignment_p->optional_fields_length += 3;
	  optional_fields += 3;
	  memcpy(optional_fields, &num_mappings, sizeof(int));
	  alignment_p->optional_fields_length += sizeof(int);
	  //alignment_p->optional_fields = (uint8_t *)optional_fields;
	  //printf("::%s\n", optional_fields);
	}
      } else {
	array_list_set_flag(0, mapping_batch_p->mapping_lists[mapping_batch_p->targets[i]]);
      }
    }

    //printf("Process Alignments %i\n", num_reads);
    /*
    fastq_read_t *read_p;
    for (i = 0; i < num_reads; i++) {
      read_p = array_list_get(i, mapping_batch_p->fq_batch);
      //printf("Read %i : %i\n", i, array_list_size(allocate_mappings[i]));
      n_total_reads_process++;
      //if (ii == i && ii != -1) { printf("1.READ %i -- Mappings found %i\n", ii, array_list_size(allocate_mappings[i]));  }
      num_mappings = alignments_filter(input_p->bwt_optarg_p->report_all, 
				       input_p->bwt_optarg_p->report_best, 
				       input_p->bwt_optarg_p->report_n_hits, 
				       mapping_batch_p->mapping_lists[i]); 
      
      //if (ii == i && ii != -1) { printf("2.Mappings found %i\n", num_mappings); exit(-1); }
      //printf("Read[%i] => %i\n", i, num_mappings);
      for (size_t m = 0; m < num_mappings; m++) {
	alignment_aux = array_list_get(m, mapping_batch_p->mapping_lists[i]);
	//printf("\t%i.%s\n", m, alignment_aux->cigar);
	if ( write_batch_p->size >= write_batch_p->allocated_size - 1) {
	  item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
	  if (time_on) { timing_stop(SW_SERVER, sw_id, timing_p); }
	  list_insert_item(item_p, write_list_p);
	  if (time_on) { timing_start(SW_SERVER, sw_id, timing_p); }
	  
	  write_batch_p = write_batch_new(write_size, MATCH_FLAG);
	  
	}

	((alignment_t **)write_batch_p->buffer_p)[write_batch_p->size] = alignment_aux;
	write_batch_p->size++;
      }

      if (!num_mappings) {
	//printf("****************** read %i NO MAPPED %i!!!\n", i, mapping_reads_p[i]);
	total_reads_unmapped++;
	read_len = strlen(read_p->sequence);
	header_len = strlen(read_p->id);
	
	if ( write_batch_p->size >= write_batch_p->allocated_size - 1) {
	  item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
	  if (time_on) { timing_stop(SW_SERVER, sw_id, timing_p); }
	  list_insert_item(item_p, write_list_p);
	  if (time_on) { timing_start(SW_SERVER, sw_id, timing_p); }
	  
	  write_batch_p = write_batch_new(write_size, MATCH_FLAG);
	}
	
	search_match_p = (char *)malloc(sizeof(char)*(read_len + 1));
	memcpy(search_match_p, read_p->sequence, read_len);
	search_match_p[read_len] = '\0';
	
	quality_match_p = (char *)malloc(sizeof(char)*(read_len + 1));
	memcpy(quality_match_p, read_p->quality, read_len);
	quality_match_p[read_len] = '\0';
	
	header_id_match_p = (char *)malloc(sizeof(char)*header_len);
	//memcpy(header_id_match_p, &header_id, len_id_header);
	len_id_header = get_to_first_blank(read_p->id, header_len, header_id_match_p);


	//	if (strncmp("@rand", sw_batch_p->allocate_reads_p[i]->id, 5)) { reads_no_random++; printf("%s\n", sw_batch_p->allocate_reads_p[i]->sequence); exit(-1);}
	
	alignment_p = alignment_new();
	cigar_p = (char *)malloc(sizeof(char)*10);
	sprintf(cigar_p, "%d%c\0", read_len, 'X');
	//TODO:chromosome 0??
	
	alignment_init_single_end(header_id_match_p, search_match_p, quality_match_p, 
				  0, -1, -1, cigar_p, 1, 0, 0, 0, 0, NULL, alignment_p);
	
	//printf("seq: %s\n", alignment_p->sequence);
	((alignment_t **)write_batch_p->buffer_p)[write_batch_p->size] = alignment_p;
	write_batch_p->size++;


	//mapping_reads_p[i] = 2;
      }
      //allocate_mappings[i]->size = 0;
      array_list_free(mapping_batch_p->mapping_lists[i], NULL);
      //array_list_free(allocate_mappings[i], NULL);     
    }

    */
    //**************************************************
    
    /*    printf("Reads mapped with global variable at this moment %d\n", reads_mapped_global);
    printf("Reads mapped with local  variable at this moment %d\n", total_reads_unmapped);
    */
    //mapping_batch_free(mapping_batch_p);
    //list_item_free(sw_item_p);
    
    if (time_on) { timing_stop(SW_SERVER, sw_id, timing_p); }
    list_insert_item(sw_item_p, alignment_list_p); 
    //printf("RNA Server process batch finish!\n");
  } // end of while
  
  array_list_free(cals_list, NULL);

  // insert or free memory
  /*  if (write_batch_p != NULL) {
    if (write_batch_p->size > 0) {
      item_p = list_item_new(0, WRITE_ITEM, write_batch_p);
      list_insert_item(item_p, write_list_p);
    } else {
      write_batch_free(write_batch_p);
    }
    }*/

  for(k=0 ; k<depth ; k++) {
    free(sw_channels_p[k].ref_p);
  }

  
  free(sw_channels_p);
  sw_simd_input_free(sw_input_p);
  sw_simd_output_free(sw_output_p);
  sw_simd_context_free(context_p);

  free(allocate_cals_fusion_p);
  free(cals_per_chromosome);
  //free(mapping_reads_p);
  
  free(reference_fusion_p);
  free(reference);

  /*
  for (i = 0; i < allocated_mapping_reads; i++){
    array_list_free(allocate_mappings[i], NULL);
  }
  */
  free(allocate_mappings);

  list_decr_writers(alignment_list_p);
  /*
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
  */
  
  LOG_DEBUG_F("rna_server (%i): END\n", omp_get_thread_num());
  
  return;
  }
