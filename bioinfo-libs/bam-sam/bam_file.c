#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "bam_file.h"


//====================================================================================
//  bam_file.c
//
//  bam file handling methods
//====================================================================================

//-----------------------------------------------------
// bam_fopen
//-----------------------------------------------------

bam_file_t* bam_fopen(char* filename) {
	return bam_fopen(filename, NULL, (char*)"r");
} 

bam_file_t* bam_fopen(char* filename, bam_header_t* bam_header_p, char* mode) {
	
	bamFile bam_fd = bam_open(filename, mode);
	if (bam_fd == NULL) {
	      char log_message[200];
	      sprintf(log_message, "Error opening file '%.150s' in mode (%s) !!!!!\n", filename, mode);
	      LOG_FATAL(log_message);
	      return NULL;
	}
	
	bam_file_t* bam_file = (bam_file_t *) malloc(sizeof(bam_file_t));
	
	bam_file->filename = filename;
	bam_file->mode = mode;
	bam_file->bam_fd = bam_fd;
	bam_file->num_alignments = 0;
	//bam_file->num_alignments = bam_file_count(bam_file);
	
	// bam header is read once and stored in read mode
	// bam header is written in write mode if is not null

	if (mode[0] == 'r') {
	      bam_file->bam_header_p = bam_header_read(bam_fd);
	} else {
	      bam_file->bam_header_p = NULL;
	}
	
	
	// ---------------- P R I N T   H E A D E R ----------------------
	
// 	printf("bam_header_p->n_targets: %i\n", bam_file->bam_header_p->n_targets);
// 	
// 	int i;
// 	for (i=0; i < bam_file->bam_header_p->n_targets; i++) {
// 	  printf("bam_file->bam_header_p->target_name[%i]: %s\n", i, bam_file->bam_header_p->target_name[i]);
// 	  printf("bam_header_p->target_len[%i]: %u\n", i, bam_file->bam_header_p->target_len[i]);
// 	}
// 
// 	printf("bam_header_p->hash: %s\n", bam_file->bam_header_p->hash);
// 	printf("bam_header_p->rg2lib: %s\n", bam_file->bam_header_p->rg2lib);
// 	printf("bam_header_p->l_text: %i\n", bam_file->bam_header_p->l_text);
// 	printf("bam_header_p->text: %s\n", bam_file->bam_header_p->text);
	
// 	FILE* ftemp = fopen("/tmp/Human_NCBI37.hbam", "w");
// 	fprintf(ftemp, bam_file->bam_header_p->text);
// 	fclose(ftemp);

// 	char* text = (char*) calloc(16035, sizeof(char));
//   
// 	FILE* fd_header = fopen("./Human_NCBI37.hbam", "r");
// 
// 	if (fd_header != NULL) {
// 	  //printf("HEADER: %s\n", fgets(text, 16035, fd_header));
// 	  fread(text, 1, 16035, fd_header);
// 	  printf("BAM HEADER: \n%s\n", text);
// 	} else {
// 	  printf("HEADER CANNOT BE READ FROM FILE\n");
// 	}
	
// 	bamFile bam_for_header_fd = bam_open("/home/vrequena/appl/ngs-gpu/bioinfo-c/ngs-hpc-pipeline/bam-hpc-tools/Human_NCBI37.hbam", "w");
// 	bam_header_write(bam_for_header_fd, bam_file->bam_header_p);
// 	bam_close(bam_for_header_fd);

	// ---------------- P R I N T   H E A D E R ----------------------

	return bam_file;
}

//-----------------------------------------------------
// bam_fread_max_size
//-----------------------------------------------------

int bam_fread_max_size(bam_batch_t* batch_p, size_t batch_size, int base_quality, bam_file_t* bam_file_p) {
  int num_alignments = 0;
  int read_bytes = 0;
  
  long current_size = 0;
  
  bam1_t* alignment_p = bam_init1();
  while ( (current_size < batch_size) && ((read_bytes = bam_read1(bam_file_p->bam_fd, alignment_p)) > 0) ) {
    
    batch_p->alignments_p[num_alignments] = alignment_p;
    current_size += read_bytes;
    num_alignments++;
    if (batch_p->allocated_alignments<=num_alignments) {
      break;
    }
    
    //printf("num_aligments: %i, before bam_destroy1...\n", num_alignments);
    //bam_destroy1(batch_p->alignments_p[num_alignments]);
    //printf("after bam_destroy1...\n");
    
    alignment_p = bam_init1();
  }
  batch_p->num_alignments = num_alignments;
  
  return num_alignments;
}

//-----------------------------------------------------
// bam_fread_max_size_no_duplicates
//-----------------------------------------------------

int bam_fread_max_size_no_duplicates(bam_batch_t* batch_p, size_t batch_size, int base_quality, bam_file_t* bam_file_p, uint8_t* prev_seq, int* prev_seq_length, int* prev_seq_start_coordinate) {
  int num_alignments = 0;
  int read_bytes = 0;

  long current_size = 0;
  
  bam1_t* alignment_p = bam_init1();

    // --------------------- T E S T --------------------
//     bam_header_t* bam_header_p;
//     bamFile bam_test_file;
//     int first_pass = 0;
//     
//     if (FILE* file_aux = fopen("/tmp/bam_test.bam", "r"))  {
//       //printf("first_pass: %i\n", first_pass);
//     } else {
//       first_pass = 1;
//     }
    // --------------------- T E S T --------------------
  
  while ( (current_size < batch_size) && ((read_bytes = bam_read1(bam_file_p->bam_fd, alignment_p)) > 0) ) {

    // --------------------- T E S T --------------------

//     if ((num_alignments < 1) && (first_pass)) {
//       bam_header_p = bam_header_new(HUMAN, NCBI37);
//       bam_test_file = bam_open("/tmp/bam_test.bam", "w");
//       bam_header_write(bam_test_file, bam_header_p);
//     }
// 
//     if ((num_alignments < 100) && (first_pass)) {      
// 
//       alignment_t* alignment_record_p = alignment_new_by_bam(alignment_p, base_quality);      
//       alignment_print(alignment_record_p);
// 
//       bam1_t* bam_p = convert_to_bam(alignment_record_p, base_quality);
//       //alignment_p = convert_to_bam(alignment_record_p);
//       bam_print(bam_p,base_quality);
//       bam_write1(bam_test_file, bam_p);
// 
//       alignment_free(alignment_record_p);
//       bam_destroy1(bam_p);
//     }
// 
//     if ((num_alignments == 100) && (first_pass)) {
//       bam_close(bam_test_file);
//     }

//    printf("bam_file_p->filename: %s\n", bam_file_p->filename);
//    printf("HEADER VALIDATION: %i\n", bam_validate_header(bam_file_p));
    
    // --------------------- T E S T --------------------
    
    if (prev_seq != NULL) {
      if (bam_batch_compare_seq(prev_seq, *prev_seq_length, *prev_seq_start_coordinate, alignment_p->data, alignment_p->core.l_qseq, alignment_p->core.pos) != 0) continue;

      //if (alignment_p->core.flag & BAM_FDUP) continue;
    }
    
    prev_seq = alignment_p->data;
    prev_seq_length = &(alignment_p->core.l_qseq);
    prev_seq_start_coordinate = &(alignment_p->core.pos);
 
    batch_p->alignments_p[num_alignments] = alignment_p;
    current_size += read_bytes;
    num_alignments++;
    if (batch_p->allocated_alignments <= num_alignments) {
      break;
    }

    alignment_p = bam_init1();    
  }
  batch_p->num_alignments = num_alignments;
  
  //free bam alignment if not null and not empty
  if ((alignment_p != NULL) && (alignment_p->core.l_qseq == 0)) {
    free(alignment_p);
  }
  
  return num_alignments;
}

//-----------------------------------------------------
// bam_fread_max_size_by_chromosome
//-----------------------------------------------------
int bam_fread_max_size_by_chromosome(bam_batch_t* batch_p, size_t batch_size, int chromosome, bam_file_t* bam_file_p) {
  int num_alignments = 0;
  int read_bytes = 0;

  long current_size = 0;

  bam1_t* alignment_p = bam_init1();
  while ( (current_size < batch_size) && ((read_bytes = bam_read1(bam_file_p->bam_fd, alignment_p)) > 0) ) {
    
    if (alignment_p->core.tid == chromosome) {
      batch_p->alignments_p[num_alignments] = alignment_p;
      current_size += read_bytes;
      num_alignments++;
      if (batch_p->allocated_alignments<=num_alignments) {
	break;
      }
      
      alignment_p = bam_init1();
    }
  }

  batch_p->num_alignments = num_alignments;
  
  return num_alignments;
}

//-----------------------------------------------------
// bam_fwrite_header
//-----------------------------------------------------

void bam_fwrite_header(bam_header_t* bam_header_p, bam_file_t* bam_file_p) {
  bam_header_write(bam_file_p->bam_fd, bam_header_p);  
}

//-----------------------------------------------------
// bam_fwrite
//-----------------------------------------------------

int bam_fwrite(bam1_t* alignment_p, bam_file_t* bam_file_p) {
  return bam_write1(bam_file_p->bam_fd, alignment_p);
}

//-----------------------------------------------------
// bam_fwrite_array
//-----------------------------------------------------

int bam_fwrite_array(bam1_t** alignment_p, int length, bam_file_t* bam_file_p) {
  
  int current_size=0;
  
  for(int i=0; i<length; i++) {
    current_size += bam_write1(bam_file_p->bam_fd, alignment_p[i]);
  }
  
  return current_size;
}

//-----------------------------------------------------
// bam_fwrite_sorted_array
//-----------------------------------------------------

int bam_fwrite_sorted_array(bam1_t** alignment_p, int* indices, int length, bam_file_t* bam_file_p) {  
  int current_size=0;

  for(int i=0; i<length; i++) {
    current_size += bam_write1(bam_file_p->bam_fd, alignment_p[indices[i]]);
  }

  return current_size; 
}

#ifdef THRUST-GPU
int bam_fwrite_sorted_array(bam1_t** alignment_p, thrust::host_vector<int> indices, int length, bam_file_t* bam_file_p) {
  
  int current_size=0;
  //int indice;
  
  for(int i=0; i<length; i++) {
    current_size += bam_write1(bam_file_p->bam_fd, alignment_p[indices[i]]);
    //indice = indices[i];
    //printf("indice: %i\n", indice);
  }
  
  return current_size; 
}
#endif

//-----------------------------------------------------
// bam_fwrite_batch
//-----------------------------------------------------

int bam_fwrite_batch(bam_batch_t* batch_p, bam_file_t* bam_file_p) {
  int current_size = 0;
  
  for(int i=0; i<batch_p->num_alignments; i++) {
    current_size += bam_fwrite(batch_p->alignments_p[i], bam_file_p);
  }
  
  return current_size;
}

//-----------------------------------------------------
// alignment_fwrite
//-----------------------------------------------------

int alignment_fwrite(alignment_t* alignment_p, bam_file_t* bam_file_p) {
  int current_size = 0;
  bam1_t* bam_p;
  
  bam_p = convert_to_bam(alignment_p, 0);  
  current_size = bam_fwrite(bam_p, bam_file_p);
  
  return current_size;
}

//-----------------------------------------------------
// alignment_fwrite_array
//-----------------------------------------------------

int alignment_fwrite_array(alignment_t** alignment_p, int length, bam_file_t* bam_file_p) {
  int current_size = 0;
  bam1_t* bam_p;
  
  for(int i=0; i<length; i++) {
    bam_p = convert_to_bam(alignment_p[i], 0);  
    current_size += bam_write1(bam_file_p->bam_fd, bam_p);
  }
  
  return current_size;  
}

//-----------------------------------------------------
// alignment_fwrite_batch
//-----------------------------------------------------

int alignment_fwrite_batch(alignment_batch_t* batch_p, bam_file_t* bam_file_p) {
  int current_size = 0;
  bam1_t* bam_p;
  
  for(int i=0; i<batch_p->num_alignments; i++) {
    bam_p = convert_to_bam(batch_p->alignments_p[i], 0);
    current_size += bam_fwrite(bam_p, bam_file_p);
  }
  
  return current_size;  
}

//-----------------------------------------------------
// bam_validate_header
//-----------------------------------------------------

int bam_validate_header(bam_file_t* bam_file_p) {
  //bam_header_t* bam_header_p = bam_header_read(bam_file_p->bam_fd);
  bam_header_read(bam_file_p->bam_fd);
//   if (bam_header_p != NULL) {
//     return 1;
//   } else {
//     return 0;
//   }  
}

//-----------------------------------------------------
// bam_fcount
//-----------------------------------------------------

unsigned int bam_fcount(bam_file_t* bam_file) {
  
  bam1_t* bam_alignment_p = (bam1_t*) calloc(1, sizeof(bam1_t));

  bam_file->num_alignments = 0;  
  
  while (bam_read1(bam_file->bam_fd, bam_alignment_p) > 0) {
    bam_file->num_alignments++;
    free(bam_alignment_p->data);
    free(bam_alignment_p);
    bam_alignment_p = (bam1_t*) calloc(1, sizeof(bam1_t));
  }
  //printf("num_alignments: %i\n", bam_file->num_alignments);
  return bam_file->num_alignments;
}


//-----------------------------------------------------
// bam_fclose
//-----------------------------------------------------

void bam_fclose(bam_file_t* bam_file) {
  if (bam_file->bam_header_p != NULL) {
    bam_header_destroy(bam_file->bam_header_p);
  }
  
  bam_close(bam_file->bam_fd);
  free(bam_file);
}

//===================================================================================
//  bam batch functions
//===================================================================================

//----------------------------------------------------------------------------------
//  bam batch new
//----------------------------------------------------------------------------------

bam_batch_t* bam_batch_new(size_t batch_size, int type) {
    bam_batch_t* batch_p = (bam_batch_t*) calloc(1, sizeof(bam_batch_t));

    batch_p->allocated_alignments = (int) (1.4 * (batch_size / MEAN_COMPRESSED_ALIGNMENT_SIZE / 12));
    batch_p->alignments_p = (bam1_t**) calloc(batch_p->allocated_alignments, sizeof(bam1_t*));
    batch_p->type = type;
    batch_p->num_alignments = 0;

    //printf("allocating bam_batch: batch=%x, batch.alignments_p=%x\n", batch_p, batch_p->alignments_p);
    return batch_p;
}

//----------------------------------------------------------------------------------
//  bam batch free
//----------------------------------------------------------------------------------

void bam_batch_free(bam_batch_t* batch_p, int free_alignments) {

  if (batch_p != NULL) {

    if (batch_p->alignments_p != NULL) {
      if (free_alignments) {
	free_bam1(batch_p->alignments_p, batch_p->num_alignments);
	batch_p->alignments_p = NULL;
      }
    }    

    free(batch_p);
  }
}

//----------------------------------------------------------------------------------

void bam_batch_print(bam_batch_t* batch_p, FILE* fd) {
  for(int i=0; i<batch_p->num_alignments; i++) {
    print_bam1(batch_p->alignments_p[i], fd);
  }
}


//----------------------------------------------------------------------------------
//  bam1_t functions
//----------------------------------------------------------------------------------

void print_bam1(bam1_t* alignment_p, FILE* fd) {
  fprintf(fd, "bam1_t:%x\ttid=%i\tpos=%i\n", alignment_p, alignment_p->core.tid, alignment_p->core.pos);
}

//----------------------------------------------------------------------------------

void free_bam1(bam1_t** alignments_p, int num_alignments) {
  if (alignments_p != NULL) {

    for(int i=0; i<num_alignments; i++) {
      if (alignments_p[i] != NULL) {
	bam_destroy1(alignments_p[i]);
      }  
    }
    
    free(alignments_p);
  }
}

//----------------------------------------------------------------------------------


//-----------------------------------------------------
// bam_batch_compare_seq
//-----------------------------------------------------

int bam_batch_compare_seq(uint8_t* data1, int length_seq1, int start_seq1, uint8_t* data2, int length_seq2, int start_seq2) {
  
  if ((length_seq1 != length_seq2) || (start_seq1 != start_seq2)) return 0;
 
  for (int i=0; i < length_seq1; i++) {
      if (bam1_seqi(data1, i) != bam1_seqi(data2, i)) return 0;    
  }
    
  return 1; 
}







