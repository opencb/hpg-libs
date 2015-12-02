#include "genome.h"
 
#define NUCLEOTIDES_NUM 4

//---------------------------------------------------------------------------------
//                     GLOBAL VARIABLES

const char NUCLEOTIDES[] = {'A', 'C', 'G', 'T'};

const unsigned int TOTAL_CODES = NUCLEOTIDES_NUM * NUCLEOTIDES_NUM * NUCLEOTIDES_NUM * NUCLEOTIDES_NUM;

//------------------------------------------------------------------------------------

genome_t* genome_new(char* sequence_filename, 
		     char* directory, int mode) {
  const int MAXLINE = 4096;
  genome_t* genome_p = (genome_t*) calloc(1, sizeof(genome_t));
  
  // genome file
  //
  size_t dna_size;
  size_t j, i;
  char path[strlen(directory) + 512];
  unsigned int offset = 0;

  // read compressed genome file
  sprintf(path, "%s/%s", directory, sequence_filename);

  LOG_DEBUG("Loading Binary DNA File");
  genome_p->X = load_binary_dna(path, &dna_size);
  LOG_DEBUG("Load DNA File Done!");

  genome_p->code_table = load_array_codes_4_chunk(genome_p);

  if (mode == BWT_MODE) {
    // read index file
    sprintf(path, "%s/index", directory);
    unsigned int chromosome, num_chromosomes = 0;

    char line[MAXLINE];
    char value[1024];

    FILE *fd = fopen(path, "r");
    if(fd == NULL) { 
      LOG_FATAL_F("FILE: '%s' not found", path);
    }

    int ch;
    do {
      ch = fgetc(fd);
      if (ch == '\n')
	num_chromosomes++;
    } while (ch != EOF);

    fseek(fd, 0, SEEK_SET);

    genome_p->num_chromosomes = num_chromosomes;
    genome_p->chr_name = (char **) calloc(num_chromosomes, sizeof(char *));
    genome_p->chr_name_length = (size_t *) calloc(num_chromosomes, sizeof(size_t));
    genome_p->chr_size = (size_t *) calloc(num_chromosomes, sizeof(size_t));
    genome_p->chr_offset = (size_t *) calloc(num_chromosomes, sizeof(size_t));

    chromosome = 0;
    while (fgets(line, MAXLINE, fd)) {
      i = 0; j= 1;
      genome_p->chr_name[chromosome] = (char *) calloc(strlen(line) + 10, sizeof(char));
      while(line[j] != ' ' ) { 
	genome_p->chr_name[chromosome][i++] = line[j++];
      }
      genome_p->chr_name[chromosome][i] = '\0';
      //printf("--> (%d):: %s\n", strlen(genome_p->chr_name[num_chromosomes]), genome_p->chr_name[num_chromosomes]);
      genome_p->chr_name_length[chromosome] = strlen(genome_p->chr_name[chromosome]);
    
      j++;
      while(line[j] != ' ') { j++; }
    
      i=0;j++;
      //printf("START %i: %c\n", j, line[j]);

      while(line[j] != '\n') { value[i++] = line[j++]; }
      value[i] = '\0';
    
      sscanf(value, "%lu", &genome_p->chr_size[chromosome]);
      genome_p->chr_offset[chromosome] = offset;
      offset += (genome_p->chr_size[chromosome] + 1);
      //printf("%i\n", genome_p->chr_size[num_chromosomes]);
      chromosome++;
    }
    fclose(fd);
  } else {
    /*sprintf(path, "%s/params.txt", directory);
    char line[2048];
    int num_chroms;
    size_t genome_len;
    FILE *f_tab = fopen(directory, "r");

    // prefix
    fgets(line, 1024, f_tab);
    // k_value
    fgets(line, 1024, f_tab);
    // pre_length
    fgets(line, 1024, f_tab);
    // A_items
    fgets(line, 1024, f_tab);
    // IA_items
    fgets(line, 1024, f_tab);
    // num_suffixes
    fgets(line, 1024, f_tab);
    // genome_length
    fgets(line, 1024, f_tab);
    genome_len = atoi(line);
    // num_chroms
    fgets(line, 1024, f_tab);
    num_chroms = atoi(line);
    
    char chrom_name[1024];
    size_t chrom_len;

    for (int i = 0; i < num_chroms; i++) {
      fgets(line, 1024, f_tab);
      sscanf(line, "%s %lu\n", chrom_name, &chrom_len);
      //printf("chrom_name: %s, chrom_len: %lu\n", chrom_name, chrom_len);
      genome_p->chr_name[i] = strdup(chrom_name);

      genome_p->chr_size[i] = chrom_len;
      genome_p->chr_offset[i] = offset;

      offset += genome_p->chr_size[i];
    }

    genome_p->num_chromosomes  = num_chroms;
    */
  }

  genome_p->genome_length = offset;
  
  return genome_p;
  
}

//-----------------------------------------------------

void genome_free(genome_t* p) {
  if (p == NULL) {
    return;
  }

  if (p->code_table != NULL){
    for (unsigned int i = 0; i < TOTAL_CODES; i++) {
      free(p->code_table[i]);
    }
    free(p->code_table);
  }

  if (p->chr_name != NULL) {
    for (unsigned int i = 0; i < p->num_chromosomes; i++) {
      free(p->chr_name[i]);
    }
    free(p->chr_name);
  }

  if (p->chr_name_length) free(p->chr_name_length);
  if (p->chr_size) free(p->chr_size);
  if (p->chr_offset) free(p->chr_offset);
  if (p->X) free(p->X);
  
  free(p);
  
}

//------------------------------------------------------------------------------------

void genome_read_sequence(char* sequence, unsigned int strand, char* chromosome, 
			  unsigned long int* start_p, unsigned long int* end_p, 
			  genome_t* genome_p) {
  unsigned int i, chr = 0;

  for(i=0 ; i< genome_p->num_chromosomes ; i++) {
    if (strcmp(chromosome, (char *) genome_p->chr_name[i]) == 0) {
      chr = i;
      break;
    }
  }

  genome_read_sequence_by_chr_index(sequence, strand, chr, start_p, end_p, genome_p);
}

//------------------------------------------------------------------------------------

char* genome_get_chr_name(unsigned int chr, unsigned int* len, genome_t* genome_p) {
  *len = genome_p->chr_name_length[chr];
  return genome_p->chr_name[chr];
}

//------------------------------------------------------------------------------------

cp_hashtable *load_hashtable_codes() {
  cp_hashtable *t = cp_hashtable_create(400,
					cp_hash_istring,
					(cp_compare_fn)strcmp); 
  
  /*const unsigned char COMBINATORIAL = (NUCLEOTIDES_NUM * NUCLEOTIDES_NUM * NUCLEOTIDES_NUM) +  (NUCLEOTIDES_NUM * NUCLEOTIDES_NUM) +  NUCLEOTIDES_NUM;
  
  unsigned char *id_array = (unsigned char *)malloc(sizeof(unsigned char)*COMBINATORIAL); 
  */

  unsigned char id = 0;
  unsigned char *value;
  char combination[4];
  
  combination[3] = '\0';
  
  for(unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++){
    combination[0] = NUCLEOTIDES[nt_1];
    for(unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++){
      combination[1] = NUCLEOTIDES[nt_2];
      for(unsigned int nt_3 = 0; nt_3 < NUCLEOTIDES_NUM; nt_3++){
	combination[2] = NUCLEOTIDES[nt_3];
	value = (unsigned char *)malloc(sizeof(unsigned char));
	*value = id;
	cp_hashtable_put(t, strdup(combination), (void *)value);
	//cp_hashtable_put(t, strdup(combination), (void *)id);
	id++;
      }
    }
  }
  
  // printf("Table size %d\n", t->table_size);
  
  combination[2] = '\0';  
  for(unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++){
    combination[0] = NUCLEOTIDES[nt_1];
    for(unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++){
      combination[1] = NUCLEOTIDES[nt_2];
      value = (unsigned char *)malloc(sizeof(unsigned char));
      *value = id;
      cp_hashtable_put(t, strdup(combination), (void *)value);
      id++;	
    }
  }
  
  //printf("Table size %d\n", t->table_size);
  
  combination[1] = '\0';
  for(unsigned int nt = 0; nt < NUCLEOTIDES_NUM; nt++){
    combination[0] = NUCLEOTIDES[nt];
    value = (unsigned char *)malloc(sizeof(unsigned char));
    *value = id;
    cp_hashtable_put(t, strdup(combination), (void *)value);
    id++;	
  }


  return t;
}


cp_hashtable *load_hashtable_codes_4_chunk() {
  cp_hashtable *t = cp_hashtable_create(512,
					cp_hash_istring,
					(cp_compare_fn)strcmp); 
  
  /*const unsigned char COMBINATORIAL = (NUCLEOTIDES_NUM * NUCLEOTIDES_NUM * NUCLEOTIDES_NUM) +  (NUCLEOTIDES_NUM * NUCLEOTIDES_NUM) +  NUCLEOTIDES_NUM;
  
  unsigned char *id_array = (unsigned char *)malloc(sizeof(unsigned char)*COMBINATORIAL); 
  */

  unsigned char id = 0;
  unsigned char *value;
  char combination[5];
  
  combination[4] = '\0';
  
  for (unsigned int nt_0 = 0; nt_0 < NUCLEOTIDES_NUM; nt_0++) {
    combination[0] = NUCLEOTIDES[nt_0];
    for (unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++) {
      combination[1] = NUCLEOTIDES[nt_1];
      for (unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++) {
	combination[2] = NUCLEOTIDES[nt_2];
	for (unsigned int nt_3 = 0; nt_3 < NUCLEOTIDES_NUM; nt_3++) {
	  combination[3] = NUCLEOTIDES[nt_3];
	  value = (unsigned char *)malloc(sizeof(unsigned char));
	  *value = id;
	  cp_hashtable_put(t, strdup(combination), (void *)value);
	  //cp_hashtable_put(t, strdup(combination), (void *)id);
	  id++;
	}
      }
    }
  }

  // printf("Table size %d\n", t->table_size);
  /*
  combination[2] = '\0';  
  for(unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++){
    combination[0] = NUCLEOTIDES[nt_1];
    for(unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++){
      combination[1] = NUCLEOTIDES[nt_2];
      value = (unsigned char *)malloc(sizeof(unsigned char));
      *value = id;
      cp_hashtable_put(t, strdup(combination), (void *)value);
      id++;	
    }
  }
  */
  //printf("Table size %d\n", t->table_size);
  /*
  combination[1] = '\0';
  for(unsigned int nt = 0; nt < NUCLEOTIDES_NUM; nt++){
    combination[0] = NUCLEOTIDES[nt];
    value = (unsigned char *)malloc(sizeof(unsigned char));
    *value = id;
    cp_hashtable_put(t, strdup(combination), (void *)value);
    id++;	
  }
  */

  return t;

}


char** load_array_codes_4_chunk(genome_t *genome) {
  char **id_array = (char **)malloc(sizeof(char *)*(TOTAL_CODES));
  
  unsigned char id = 0;
  char combination[5];
  
  combination[4] = '\0';
  
  printf("TOTAL CODES = %i\n", TOTAL_CODES);
  //exit(-1);

  for (unsigned int nt_0 = 0; nt_0 < NUCLEOTIDES_NUM; nt_0++) {
    combination[0] = NUCLEOTIDES[nt_0];
    for (unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++) {
      combination[1] = NUCLEOTIDES[nt_1];
      for (unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++) {
	combination[2] = NUCLEOTIDES[nt_2];
	for (unsigned int nt_3 = 0; nt_3 < NUCLEOTIDES_NUM; nt_3++) {
	  combination[3] = NUCLEOTIDES[nt_3];
	  if (id >= TOTAL_CODES) {
	    printf("%i >= %i\n", id, TOTAL_CODES);
	    exit(-1);
	  }
	  id_array[id] = strdup(combination);
	  id++;
	}
      }
    }
  }
  
  genome->code_table = id_array;

  return id_array;

  /*  
  // printf("Table size %d\n", t->table_size);
  combination[2] = '\0';
  //combination[2] = '\0';  
  for(unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++){
    combination[0] = NUCLEOTIDES[nt_1];
    for(unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++){
      combination[1] = NUCLEOTIDES[nt_2];
      id_array[id] = strdup(combination);
      id++;	
    }
  }
  
  //printf("Table size %d\n", t->table_size);
  
  combination[1] = '\0';
  //id_array[id][1] = '\0';

  for(unsigned int nt = 0; nt < NUCLEOTIDES_NUM; nt++){
    combination[0] = NUCLEOTIDES[nt];
    id_array[id] = strdup(combination);
    id++;	
  }

  return id_array;
  */
}


void code_binary_file_generator(size_t chunk_size, char *dna_filename, char *dna_binary_filename,  cp_hashtable *t){
  
  if (chunk_size <= 0) { 
    chunk_size = 100000000; //100MB 
  }

  FILE *binary_fd, *fd;
  fd = fopen(dna_filename, "r");
  if (fd == NULL) {  LOG_FATAL_F("Error opening file %s", dna_filename); }

  binary_fd = fopen (dna_binary_filename, "wb");
  if (binary_fd == NULL) { LOG_FATAL_F("Error opening file %s", dna_binary_filename); }

  char *dna_chunk = (char *)malloc(sizeof(char)*chunk_size);
  
  size_t codes_allocate = chunk_size;
  unsigned char *code_values = (unsigned char *)malloc(sizeof(unsigned char)*codes_allocate);
  size_t code_pos = 0;

  size_t dna_len;
  char key[4];
  unsigned char max_chunk = 3;
  unsigned char actual_nt = 0;
  char *res;
  unsigned char value;
  unsigned char *value_ptr;

  fseek(fd, 0L, SEEK_END);
  size_t fd_size = ftell(fd);
  fseek(fd, 0L, SEEK_SET);

  char *genome_txt = (char *)malloc(sizeof(char)*fd_size);
  size_t actual_pos = 0;

  while (!feof(fd)) {
    res = fgets(dna_chunk, chunk_size, fd);
    if (dna_chunk[0] != '>') {
      dna_len = strlen(dna_chunk);
      //printf("Process (%i): %s", dna_len, dna_chunk);
      for (unsigned int c = 0; c < dna_len; c++) {
	if (dna_chunk[c] != '\n') {

	  if (dna_chunk[c] == 'n' || 
	      dna_chunk[c] == 'N') {
	    dna_chunk[c] = 'A';
	  }

	  //printf("Char(%i)[%i]: %c\n", c, actual_nt, dna_chunk[c]);
	  if (dna_chunk[c] == 'a' || 
	      dna_chunk[c] == 'c' || 
	      dna_chunk[c] == 'g' || 
	      dna_chunk[c] == 't' ) {
	    //printf("Convert %c in %c\n", dna_chunk[c], dna_chunk[c] - 32);
	    dna_chunk[c] = dna_chunk[c] - 32;
	  }
	} else {
	  dna_chunk[c] = '\0';
	}
      }
      size_t len = strlen(dna_chunk);
      memcpy(&genome_txt[actual_pos], dna_chunk, len);
      actual_pos += len;
    }
  }

  for (size_t c = 0; c < actual_pos; c++) {
    key[actual_nt++] = genome_txt[c];
    if (actual_nt ==  max_chunk) {
      key[actual_nt] = '\0';
      //printf("Store: %s\n", key);
      //value = cp_hashtable_get(t, key);
      //printf("1.%p\n", value_ptr);
      
      value_ptr = cp_hashtable_get(t, key);
      //printf("Value : %i\n", cp_hashtable_get(t, key));
      if (value_ptr == NULL) { exit(-1); }
      //printf("2.%p\n", value_ptr);
      value = *value_ptr;
      
      code_values[code_pos++] = value;
      //printf("Stored code %d == %s : %d\n", value, key, code_pos - 1);
      if (code_pos >= codes_allocate) {
	//printf("Write Ids in file...\n");
	fwrite(code_values, sizeof(unsigned char), code_pos, binary_fd);
	code_pos = 0;
      }
      actual_nt = 0;
    }
  }

  printf("end\n");
    
  if(actual_nt > 0) {
    key[actual_nt] = '\0';
    //printf("Store: %s\n", key);
    value_ptr = cp_hashtable_get(t, key);
    value = *value_ptr;
    //value = (unsigned char)cp_hashtable_get(t, key);
    //value = *value_ptr;
    code_values[code_pos++] = value;	
  }

  if (code_pos >= 0) {
    fwrite(code_values, sizeof(unsigned char), code_pos, binary_fd);
    code_pos = 0;
  }
      
  fclose(fd);
  fclose(binary_fd);
}


void code_binary_file_generator_4_chunk(size_t chunk_size, char *dna_filename, char *dna_binary_filename,  cp_hashtable *t){
  
  FILE *binary_fd, *fd;
  fd = fopen(dna_filename, "r");
  if (fd == NULL) {  LOG_FATAL_F("Error opening file %s", dna_filename); }

  binary_fd = fopen (dna_binary_filename, "wb");
  if (binary_fd == NULL) { LOG_FATAL_F("Error opening file %s", dna_binary_filename); }

  char *dna_chunk = (char *)malloc(sizeof(char)*chunk_size);
  
  size_t code_pos = 0;

  size_t dna_len;
  char key[5];
  unsigned char max_chunk = 4;
  unsigned char actual_nt = 0;
  char *res;
  unsigned char value;
  unsigned char *value_ptr;

  fseek(fd, 0L, SEEK_END);
  size_t fd_size = ftell(fd);
  fseek(fd, 0L, SEEK_SET);

  char *genome_txt = (char *)malloc(sizeof(char)*fd_size);
  size_t actual_pos = 0;

  while (!feof(fd)) {
    res = fgets(dna_chunk, chunk_size, fd);
    if (dna_chunk[0] != '>') {
      dna_len = strlen(dna_chunk);
      for (unsigned int c = 0; c < dna_len; c++) {
	if (dna_chunk[c] != '\n') {
	  if (dna_chunk[c] == 'n' || 
	      dna_chunk[c] == 'N') {
	    dna_chunk[c] = 'A';
	  }
	  if (dna_chunk[c] == 'a' || 
	      dna_chunk[c] == 'c' || 
	      dna_chunk[c] == 'g' || 
	      dna_chunk[c] == 't' ) {
	    dna_chunk[c] = dna_chunk[c] - 32;
	  }
	} else {
	  dna_chunk[c] = '\0';
	}
      }
      size_t len = strlen(dna_chunk);
      memcpy(&genome_txt[actual_pos], dna_chunk, len);
      actual_pos += len;
    }
  }

  printf("Reading file end!!\n");
  

  size_t codes_allocate = actual_pos;//chunk_size;
  unsigned char *code_values = (unsigned char *)malloc(sizeof(unsigned char)*codes_allocate);
  
  for (size_t c = 0; c < actual_pos; c++) {
    key[actual_nt++] = genome_txt[c];
    if (actual_nt ==  max_chunk) {
      key[actual_nt] = '\0';
      value_ptr = cp_hashtable_get(t, key);      
      if (value_ptr == NULL) { printf("Error encoding: %s\n", key); exit(-1); }
      value = *value_ptr;
      code_values[code_pos++] = value;
      actual_nt = 0;
    }
  }
  
  fwrite(code_values, sizeof(unsigned char), code_pos, binary_fd);

  fclose(fd);
  fclose(binary_fd);

  /*
  //Test Process
  //===============================================================================//
  size_t genome_len = 3095693982;
  FILE *fd_S = fopen("/data/hpg/sa_18/hs.73.fa.S", "r");
  char *S = (char *) malloc(genome_len);  
  fread(S, sizeof(char), genome_len, fd_S);
  
  char ref[1024];
  char ref_ok[1024];
  char ref_ok_2[1024];
  
  genome_t *genome = genome_new(dna_binary_filename, "/home/hmartinez/appl.SA/hpg-aligner", SA_MODE);

  size_t start = 500000000;
  size_t end = 500000100;

  strncpy(ref_ok, &genome_txt[start], end - start);
  strncpy(ref_ok_2, &S[start], end - start);
  genome_read_sequence_sa(ref, &start, &end, genome);

  printf("Total Size %lu: %s\n", actual_pos, ref);
  printf("Total Size %lu: %s\n", actual_pos, ref_ok);
  printf("Total Size %lu: %s\n", actual_pos, ref_ok_2);

  start = 3095690000;
  end = 3095690100;

  strncpy(ref_ok, &genome_txt[start], end - start);
  strncpy(ref_ok_2, &S[start], end - start);
  genome_read_sequence_sa(ref, &start, &end, genome);
  ref[start - end] = '\0';
  
  printf("Total Size %lu: %s\n", actual_pos, ref);
  printf("Total Size %lu: %s\n", actual_pos, ref_ok);
  printf("Total Size %lu: %s\n", actual_pos, ref_ok_2);

  start = 1095693000;
  end = 1095693100;

  strncpy(ref_ok, &genome_txt[start], end - start);
  strncpy(ref_ok_2, &S[start], end - start);
  genome_read_sequence_sa(ref, &start, &end, genome);
  ref[start - end] = '\0';
  
  printf("Total Size %lu: %s\n", actual_pos, ref);
  printf("Total Size %lu: %s\n", actual_pos, ref_ok);
  printf("Total Size %lu: %s\n", actual_pos, ref_ok_2);

  start = 20000;
  end = 20100;

  strncpy(ref_ok, &genome_txt[start], end - start);
  strncpy(ref_ok_2, &S[start], end - start);
  genome_read_sequence_sa(ref, &start, &end, genome);
  ref[start - end] = '\0';
  
  printf("Total Size %lu: %s\n", actual_pos, ref);
  printf("Total Size %lu: %s\n", actual_pos, ref_ok);
  printf("Total Size %lu: %s\n", actual_pos, ref_ok_2);

  //===============================================================================//
  */
}


unsigned char *load_binary_dna(char *dna_binary_filename, size_t *size){
  FILE *binary_fd = fopen (dna_binary_filename, "rb");
  if (!binary_fd) {
    LOG_FATAL_F("Error to opening '%s' file\n", dna_binary_filename);
  }

  struct stat st;                                                                                                                                                 
  stat(dna_binary_filename, &st);     
  *size = st.st_size;                                                                                                               
  LOG_DEBUG_F("Size File %lu\n", *size);

  unsigned char *dna_encoding = (unsigned char *)malloc(sizeof(unsigned char)*(*size));

  int res = fread(dna_encoding, sizeof(unsigned char), *size, binary_fd);
  if (!res) {
    LOG_FATAL_F("Error, '%s' file is empty!\n", dna_binary_filename);
  }

  fclose(binary_fd);

  return dna_encoding;
}

void get_genome_sequence(char *sequence, unsigned int chromosome, unsigned int strand, size_t start, size_t end, char **array_codes, char *dna_encoding){
  size_t group_start = start/3;
  size_t group_end   = end/3;

  size_t nucleotide_start = start%3;
  size_t nucleotide_end   = end%3;

  //size_t total_nt = end - start;
  //size_t total_groups = total_nt/3;
  //size_t rest_groups = total_nt%3;

  unsigned int seq_pos = 0;
  


  
  if (start >= end)
    return;

  for(unsigned int i = nucleotide_start; i < 3; i++){
    sequence[seq_pos++] = array_codes[(unsigned char)dna_encoding[group_start]][i];
  }

  group_start++;
  
  while (group_start < group_end) {
    sequence[seq_pos++] = array_codes[(unsigned char)dna_encoding[group_start]][0];
    sequence[seq_pos++] = array_codes[(unsigned char)dna_encoding[group_start]][1];
    sequence[seq_pos++] = array_codes[(unsigned char)dna_encoding[group_start]][2];
    group_start++;
  }

  if (group_start <= group_end) {
    for(unsigned int i = 0; i <= nucleotide_end; i++){
      sequence[seq_pos++] = array_codes[(unsigned char)dna_encoding[group_start]][i];
    }
  }

  sequence[seq_pos] = '\0';
}


void genome_read_sequence_by_chr_index(char* sequence, unsigned int strand, 
				       unsigned int chr, size_t *start_p,
				       size_t *end_p, genome_t* genome_p) {

  size_t s, e;
  
  //if (*start_p < 1) (*start_p) = 1;
  //if (*start_p > genome_p->chr_size[chr]) (*start_p) = genome_p->chr_size[chr];
  //if (*end_p < 1) (*end_p) = 1;
  //if (*end_p > genome_p->chr_size[chr]) (*end_p) = genome_p->chr_size[chr];

  s = (*start_p) + genome_p->chr_offset[chr] - 1;
  e = (*end_p) + genome_p->chr_offset[chr] - 1;

  size_t group_start = s/4;
  size_t group_end   = e/4;

  size_t nucleotide_start = s%4;
  size_t nucleotide_end   = e%4;

  unsigned int seq_pos = 0;
  char *str_tmp;
  
  //printf("Get sequence 1\n");
  str_tmp = genome_p->code_table[genome_p->X[group_start]];
  for(unsigned int i = nucleotide_start; i < 4; i++){
    sequence[seq_pos++] = str_tmp[i];
  }
  group_start++;

  //printf("Get sequence 2\n");
  while (group_start < group_end) {
    str_tmp = genome_p->code_table[genome_p->X[group_start]];
    sequence[seq_pos++] = str_tmp[0];
    sequence[seq_pos++] = str_tmp[1];
    sequence[seq_pos++] = str_tmp[2];
    sequence[seq_pos++] = str_tmp[3];
    group_start++;
  }

  //printf("Get sequence 3\n");
  if (group_start <= group_end) {  
    str_tmp = genome_p->code_table[genome_p->X[group_start]];
    for(unsigned int i = 0; i <= nucleotide_end; i++){
      sequence[seq_pos++] = str_tmp[i];
    }
  }


  sequence[seq_pos] = '\0';


  /*
  size_t group_start = s/3;
  size_t group_end   = e/3;

  size_t nucleotide_start = s%3;
  size_t nucleotide_end   = e%3;

  unsigned int seq_pos = 0;
  
  //printf("Genome start %lu and genome end %lu = %lu\n", s, e, e - s);
  for(unsigned int i = nucleotide_start; i < 3; i++){
    sequence[seq_pos++] = genome_p->code_table[genome_p->X[group_start]][i];
  }

  group_start++;
  
  while (group_start < group_end) {
    sequence[seq_pos++] = genome_p->code_table[genome_p->X[group_start]][0];
    sequence[seq_pos++] = genome_p->code_table[genome_p->X[group_start]][1];
    sequence[seq_pos++] = genome_p->code_table[genome_p->X[group_start]][2];
    group_start++;
  }
  
  if (group_start <= group_end) {
    for(unsigned int i = 0; i <= nucleotide_end; i++){
      sequence[seq_pos++] = genome_p->code_table[genome_p->X[group_start]][i];
    }
  }

  sequence[seq_pos] = '\0';
  */

  

}

void genome_read_sequence_sa(char* sequence, size_t *start_p,
			     size_t *end_p, genome_t* genome_p) {
  
  //printf("Call genome %lu - %lu < %lu ??\n", *start_p, *end_p, genome_p->genome_length);
  
  size_t s, e;
  
  s = (*start_p);
  e = (*end_p);
  
  size_t group_start = s/4;
  size_t group_end   = e/4;

  size_t nucleotide_start = s%4;
  size_t nucleotide_end   = e%4;

  unsigned int seq_pos = 0;
  char *str_tmp;
  
  //printf("Get sequence 1\n");
  str_tmp = genome_p->code_table[genome_p->X[group_start]];
  for(unsigned int i = nucleotide_start; i < 4; i++){
    sequence[seq_pos++] = str_tmp[i];
  }
  group_start++;

  //printf("Get sequence 2\n");
  while (group_start < group_end) {
    str_tmp = genome_p->code_table[genome_p->X[group_start]];
    sequence[seq_pos++] = str_tmp[0];
    sequence[seq_pos++] = str_tmp[1];
    sequence[seq_pos++] = str_tmp[2];
    sequence[seq_pos++] = str_tmp[3];
    group_start++;
  }

  //printf("Get sequence 3\n");
  if (group_start <= group_end) {  
    str_tmp = genome_p->code_table[genome_p->X[group_start]];
    for(unsigned int i = 0; i <= nucleotide_end; i++){
      sequence[seq_pos++] = str_tmp[i];
    }
  }


  sequence[seq_pos] = '\0';

}

void genome_read_sequence_sa_backup(char* sequence, size_t *start_p,
				    size_t *end_p, genome_t* genome_p) {

  //printf("Call genome %lu - %lu < %lu ??\n", *start_p, *end_p, genome_p->genome_length);

  size_t s, e;
  
  s = (*start_p);
  e = (*end_p);
  
  size_t group_start = s/3;
  size_t group_end   = e/3;

  size_t nucleotide_start = s%3;
  size_t nucleotide_end   = e%3;

  unsigned int seq_pos = 0;
  char *str_tmp;

  str_tmp = genome_p->code_table[genome_p->X[group_start]];
  for(unsigned int i = nucleotide_start; i < 3; i++){
    sequence[seq_pos++] = str_tmp[i];
  }
  group_start++;

  while (group_start < group_end) {
    str_tmp = genome_p->code_table[genome_p->X[group_start]];
    sequence[seq_pos++] = str_tmp[0];
    sequence[seq_pos++] = str_tmp[1];
    sequence[seq_pos++] = str_tmp[2];
    group_start++;
  }

  if (group_start <= group_end) {  
    str_tmp = genome_p->code_table[genome_p->X[group_start]];
    for(unsigned int i = 0; i <= nucleotide_end; i++){
      sequence[seq_pos++] = str_tmp[i];
    }
  }


  sequence[seq_pos] = '\0';

}


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------                                                                                  
void generate_codes(char *dna_binary_filename, char *dna_filename){
  LOG_DEBUG("Loading hashtable Codes ...\n");
  cp_hashtable *t = load_hashtable_codes();
  LOG_DEBUG("Loading done!\n");

  LOG_DEBUG("Genrate Binary Genome File...\n");
  code_binary_file_generator(100000000, dna_filename,
                             dna_binary_filename, t);
  LOG_DEBUG("Generate done! Happy usage!\n");
}

void generate_codes_4_chunk(char *dna_binary_filename, char *dna_filename){
  LOG_DEBUG("Loading hashtable Codes ...\n");
  cp_hashtable *t = load_hashtable_codes_4_chunk();
  LOG_DEBUG("Loading done!\n");

  LOG_DEBUG("Genrate Binary Genome File...\n");
  code_binary_file_generator_4_chunk(100000000, dna_filename,
				     dna_binary_filename, t);
  LOG_DEBUG("Generate done! Happy usage!\n");
}

//------------------------------------------------------------------------------------
