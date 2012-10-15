#include <math.h>
#include <sys/stat.h>

#include "genome.h"

#define NUCLEOTIDES_NUM 5

//---------------------------------------------------------------------------------
//                     GLOBAL VARIABLES

const char NUCLEOTIDES[] = {'A', 'C', 'G', 'N', 'T'};
const unsigned char TOTAL_CODES = pow(NUCLEOTIDES_NUM, 3) + pow(NUCLEOTIDES_NUM, 2) + NUCLEOTIDES_NUM; 

//------------------------------------------------------------------------------------
/*
genome_t* genome_new(char* sequence_filename, char* chromosome_filename) {
  const int MAXLINE = 1024;
  genome_t* genome_p = (genome_t*) calloc(1, sizeof(genome_t));
  // genome file
  //
  size_t dna_size;

  printf("Loading Binary File\n");
  genome_p->X = load_binary_dna(sequence_filename, &dna_size);
  printf("Load done!\n");
  genome_p->code_table = load_array_codes();
  
  // read index file
  unsigned int num_chromosomes = 0;
  unsigned int offset = 0;

  char* p;
  char line[MAXLINE];
  FILE *fd = fopen("/home/hmartinez/genome/human64/chromosome_index.txt", "r");
  while (fgets(line, MAXLINE, fd) ) {
    p = strrchr(line, '\n'); *p = '\0';
    p = strrchr(line, '\t'); *p = '\0';
    //printf("%s\n", line);
    strcpy((char*) genome_p->chr_name[num_chromosomes], line);
    genome_p->chr_name_length[num_chromosomes] = strlen(genome_p->chr_name[num_chromosomes]);
    sscanf(p+1, "%i", &genome_p->chr_size[num_chromosomes]);
    genome_p->chr_offset[num_chromosomes] = offset;

    offset += (genome_p->chr_size[num_chromosomes] + 1);
    num_chromosomes++;
  }
  fclose(fd);
  //printf("In genome %d chromosome and %d\n", num_chromosomes, offset);
  genome_p->num_chromosomes = num_chromosomes;

  return genome_p;
}*/

genome_t* genome_new(char* sequence_filename, char* directory) {
  const int MAXLINE = 1024;
  genome_t* genome_p = (genome_t*) calloc(1, sizeof(genome_t));
  // genome file
  //
  size_t dna_size;
  size_t j, i;
  char path[strlen(directory) + 512];
  
  strcpy(path, directory);
  strcat(path, "/index");

  printf("Loading Binary File\n");
  genome_p->X = load_binary_dna(sequence_filename, &dna_size);
  printf("Load done!\n");
  genome_p->code_table = load_array_codes();

  // read index file
  //
  unsigned int num_chromosomes = 0;
  unsigned int offset = 0;
  char* p;
  char line[MAXLINE];
  char value[1024];
  FILE *fd = fopen(path, "r");
  if(fd == NULL) { printf("FILE: '%s' not found\n", path);exit(-1); }
  //printf("%s\n", path);
  //FILE *fd = fopen("/home/hmartinez/BenchMarks/HomoSapiens_BWT_Index/chromosome_index.txt", "r");
  while (fgets(line, MAXLINE, fd) ) {
    //printf("%s\n", line);
    i = 0; j= 1;
    while(line[j] != ' ' ){genome_p->chr_name[num_chromosomes][i++] = line[j++];}
    genome_p->chr_name[num_chromosomes][i] = '\0';
    //printf("--> (%d):: %s\n", strlen(genome_p->chr_name[num_chromosomes]), genome_p->chr_name[num_chromosomes]);
    genome_p->chr_name_length[num_chromosomes] = strlen(genome_p->chr_name[num_chromosomes]);
    
    j++;
    while(line[j] != ' '){j++;}
    
    i=0;j++;
    //printf("START %i: %c\n", j, line[j]);

    while(line[j] != '\n'){value[i++] = line[j++];}
    value[i] = '\0';
    
    sscanf(value, "%i", &genome_p->chr_size[num_chromosomes]);
    genome_p->chr_offset[num_chromosomes] = offset;
    offset += (genome_p->chr_size[num_chromosomes] + 1);
    //printf("%i\n", genome_p->chr_size[num_chromosomes]);
    num_chromosomes++;
  }

  fclose(fd);
  
  //printf("In genome %d chromosome and %d\n", num_chromosomes, offset);

  genome_p->num_chromosomes = num_chromosomes;

  return genome_p;
}

//-----------------------------------------------------

void genome_free(genome_t* genome_p) {
  if (genome_p == NULL) {
    return;
  }

  if (genome_p->code_table != NULL){
    for (unsigned int i = 0; i < TOTAL_CODES; i++){
      free(genome_p->code_table[i]);
    }
    free(genome_p->code_table);
  }

  if (genome_p->X != NULL) free(genome_p->X);

  free(genome_p);
}

//------------------------------------------------------------------------------------

void genome_read_sequence(char* sequence, unsigned int strand, char* chromosome, 
			  unsigned long int* start_p, unsigned long int* end_p, 
			  genome_t* genome_p) {
  unsigned int i, chr;

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

cp_hashtable *load_hasthable_codes() {
  cp_hashtable *t = cp_hashtable_create(400,
					cp_hash_istring,
					strcasecmp); 
  
  /*const unsigned char COMBINATORIAL = (NUCLEOTIDES_NUM * NUCLEOTIDES_NUM * NUCLEOTIDES_NUM) +  (NUCLEOTIDES_NUM * NUCLEOTIDES_NUM) +  NUCLEOTIDES_NUM;
  
  unsigned char *id_array = (unsigned char *)malloc(sizeof(unsigned char)*COMBINATORIAL); 
  */
  unsigned char id = 0;
  char combination[4];
  
  combination[3] = '\0';
  
  for(unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++){
    combination[0] = NUCLEOTIDES[nt_1];
    for(unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++){
      combination[1] = NUCLEOTIDES[nt_2];
      for(unsigned int nt_3 = 0; nt_3 < NUCLEOTIDES_NUM; nt_3++){
	combination[2] = NUCLEOTIDES[nt_3];
	cp_hashtable_put(t, strdup(combination), id);
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
      cp_hashtable_put(t, strdup(combination), id);
      id++;	
    }
  }
  
  //printf("Table size %d\n", t->table_size);
  
  combination[1] = '\0';
  for(unsigned int nt = 0; nt < NUCLEOTIDES_NUM; nt++){
    combination[0] = NUCLEOTIDES[nt];
    cp_hashtable_put(t, strdup(combination), id);
    id++;	
  }

  return t;

}


char** load_array_codes() {
  char **id_array = (char **)malloc(sizeof(char *)*TOTAL_CODES);
  
  unsigned char id = 0;
  char combination[4];
  
  combination[3] = '\0';
  
  for(unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++){
    combination[0] = NUCLEOTIDES[nt_1];
    for(unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++){
      combination[1] = NUCLEOTIDES[nt_2];
      for(unsigned int nt_3 = 0; nt_3 < NUCLEOTIDES_NUM; nt_3++){
	combination[2] = NUCLEOTIDES[nt_3];
	id_array[id] = strdup(combination);
	//printf("id_array=%s - combination=%s\n", id_array[id], combination);
	id++;
      }
    }
  }
  
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
}


void code_binary_file_generator(size_t chunk_size, char *dna_filename, char *dna_binary_filename,  cp_hashtable *t){
  
  if (chunk_size <= 0) { 
    chunk_size = 100000000; //100MB 
  }

  FILE *binary_fd, *fd;
  fd = fopen(dna_filename, "r");
  if (fd == NULL) {  printf("Error opening file %s\n", dna_filename); exit(-1); }

  binary_fd = fopen (dna_binary_filename, "wb");
  if (binary_fd == NULL) { printf("Error opening file %s\n", dna_binary_filename); exit(-1); }

  char *dna_chunk = (char *)malloc(sizeof(char)*chunk_size);
  
  size_t codes_allocate = chunk_size;
  unsigned char *code_values = (unsigned char *)malloc(sizeof(unsigned char)*codes_allocate);
  size_t code_pos = 0;

  size_t dna_len;
  char key[4];
  unsigned char max_chunk = 3;
  unsigned char actual_nt = 0;
  
  unsigned char value;
  size_t nt = 0;
  
  /*
  while (!feof(fd)) {
    fgets(dna_chunk, chunk_size, fd);
    dna_len = strlen(dna_chunk);
    while(nt < dna_len) {
      while((actual_nt < key_chunk) && (dna_chunk[nt] != '\n') && (dna_chunk[nt] != '\0') && (nt < dna_len)){
	key[actual_nt] = dna_chunk[nt];
	nt++;
	actual_nt++;
      }

      key[actual_nt] = '\0';
      actual_nt = 0;
      value = (unsigned char)cp_hashtable_get(t, key);
      code_values[code_pos] = value; 
      code_pos++;
   
      if (code_pos >= codes_allocate) {
	fwrite(code_values, sizeof(unsigned char), code_pos, binary_fd);
	code_pos = 0;
      }
      //printf("Key:%s - Value:%i\n", key, value);
    }
    }*/

  while (!feof(fd)) {
    fgets(dna_chunk, chunk_size, fd);
    dna_len = strlen(dna_chunk);
    for (unsigned int c = 0; c < dna_len; c++) {
      key[actual_nt++] = dna_chunk[c];
      if (actual_nt ==  max_chunk){
	key[actual_nt] = '\0';
	value = (unsigned char)cp_hashtable_get(t, key);
	code_values[code_pos++] = value;
	//printf("Stored code %d == %s : %d\n", value, key, code_pos - 1);
	if (code_pos >= codes_allocate) {
	  printf("Write Ids...\n");
	  fwrite(code_values, sizeof(unsigned char), code_pos, binary_fd);
	  code_pos = 0;
	}
   	actual_nt = 0;
      }
    }
  }

  if(actual_nt > 0){
    key[actual_nt] = '\0';
    value = (unsigned char)cp_hashtable_get(t, key);
    code_values[code_pos++] = value;	
  }

  if (code_pos >= 0) {
    fwrite(code_values, sizeof(unsigned char), code_pos, binary_fd);
    code_pos = 0;
  }
      
  fclose(fd);
  fclose(binary_fd);
}

unsigned char *load_binary_dna(char *dna_binary_filename, size_t *size){
  FILE *binary_fd = fopen (dna_binary_filename, "rb");
  
  struct stat st;                                                                                                                                                     
  stat(dna_binary_filename, &st);     
  *size = st.st_size;                                                                                                                                           
  printf("size file=%i\n", *size);

  unsigned char *dna_encoding = (unsigned char *)malloc(sizeof(unsigned char)*(*size));
  fread(dna_encoding, sizeof(unsigned char), *size, binary_fd);
  
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
  
  char nt_code[4];
  unsigned char id;
  
  if (start >= end)
    return;

  for(unsigned int i = nucleotide_start; i < 3; i++){
    sequence[seq_pos++] = array_codes[dna_encoding[group_start]][i];
  }

  group_start++;
  
  while (group_start < group_end) {
    sequence[seq_pos++] = array_codes[dna_encoding[group_start]][0];
    sequence[seq_pos++] = array_codes[dna_encoding[group_start]][1];
    sequence[seq_pos++] = array_codes[dna_encoding[group_start]][2];
    group_start++;
  }

  if (group_start <= group_end) {
    for(unsigned int i = 0; i <= nucleotide_end; i++){
      sequence[seq_pos++] = array_codes[dna_encoding[group_start]][i];
    }
  }

  sequence[seq_pos] = '\0';
}


void genome_read_sequence_by_chr_index(char* sequence, unsigned int strand, 
				       unsigned int chr, size_t *start_p,
				       size_t *end_p, genome_t* genome_p) {

  size_t s, e;
  
  if (*start_p < 1) (*start_p) = 1;
  if (*start_p > genome_p->chr_size[chr]) (*start_p) = genome_p->chr_size[chr];
  if (*end_p < 1) (*end_p) = 1;
  if (*end_p > genome_p->chr_size[chr]) (*end_p) = genome_p->chr_size[chr];

  s = (*start_p) + genome_p->chr_offset[chr] - 1;
  e = (*end_p) + genome_p->chr_offset[chr] - 1;
  
  size_t group_start = s/3;
  size_t group_end   = e/3;

  size_t nucleotide_start = s%3;
  size_t nucleotide_end   = e%3;

  unsigned int seq_pos = 0;
  
  char nt_code[4];
  unsigned char id;
  
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

}

/*
>>>>>>> 19eebf394c74de43c96416a96971dce85ac38d21
void genome_read_sequence_by_chr_index(char* sequence, unsigned int strand, unsigned int chr, unsigned long int* start_p, unsigned long int* end_p, genome_t* genome_p) {

  unsigned long int i, j, s, e;
  
  if (*start_p < 1) (*start_p) = 1;
  if (*start_p > genome_p->chr_size[chr]) (*start_p) = genome_p->chr_size[chr];
  if (*end_p < 1) (*end_p) = 1;
  if (*end_p > genome_p->chr_size[chr]) (*end_p) = genome_p->chr_size[chr];

  s = (*start_p) + genome_p->chr_offset[chr] - 1;
  e = (*end_p) + genome_p->chr_offset[chr] - 1;
  j = 0;

  if (strand == 0) {
    for(i=s ; i<=e ; i++) {
      sequence[j++] = genome_p->X[i];
    }
  } else {
    char c;
    for(i=e ; i>=s ; i--) {
      c = genome_p->X[i];
      if      (c == 'A' || c == 'a') { sequence[j++] = 'T';  } 
      else if (c == 'C' || c == 'c') { sequence[j++] = 'G';  } 
      else if (c == 'G' || c == 'g') { sequence[j++] = 'C';  }
      else if (c == 'T' || c == 't') { sequence[j++] = 'A';  } 
      else         	             { sequence[j++] = 'N'; }
    }    
  }
  
  sequence[j]='\0';

}

<<<<<<< HEAD
//------------------------------------------------------------------------------------

char* genome_get_chr_name(unsigned int chr, unsigned int* len, genome_t* genome_p) {
  *len = genome_p->chr_name_length[chr];
  return genome_p->chr_name[chr];
}

//------------------------------------------------------------------------------------
=======
*/

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------                                                                                  
void generate_codes(char *dna_binary_filename, char *dna_filename){
  printf("Loading hashtable Codes ...\n");
  cp_hashtable *t = load_hasthable_codes();
  printf("Loading done!\n");

  printf("Genrate Binary Genome File...\n");
  code_binary_file_generator(100000000, dna_filename,
                             dna_binary_filename, t);
  printf("Generate done! Happy usage!\n");
}
/*
void generate_codes(){
  char *dna_binary_filename = "dna_compression.bin";
  char *dna_filename = "/home/hmartinez/BenchMarks/HomoSapiens_BWT_Index/all.oneline";//"pruebecilla.fa";
  size_t dna_size;
  printf("Loading hashtable...\n");
  //  cp_hashtable *t = load_hasthable_codes();
  printf("Loading done!\n");
  printf("Genrate Binary Genome File...\n");
  //code_binary_file_generator(100000000, dna_filename, dna_binary_filename, t);
  printf("Generate done\n");
  printf("Loading Binary File\n");
  unsigned char *dna_encoding = load_binary_dna(dna_binary_filename, &dna_size);
  printf("Load done!\n");

  char **array_codes = load_array_codes();
  
  size_t start = 10010;
  size_t end = 10050;

  char *sequence = (char *)malloc(sizeof(char)*(end - start + 1));

  get_genome_sequence(sequence, 0, 1, start, end, array_codes, dna_encoding);
  printf("(%i-%i) = %s\n", start, end, sequence);

  //for (int i = 0; i < dna_size; i++) { printf("%s", array_codes[dna_encoding[i]] ); }

  /*
  for (unsigned int i = 0; i < dna_size; i++) {
    printf("%i-", dna_encoding[i]);
    }*/
  
  /* unsigned char **values = (unsigned char **)cp_hashtable_get_values(t);

  for(int i = 0; i < 155; i++){
    printf("%i\n", values[i]);
  }
  exit(-1);
  
 
  }*/
//------------------------------------------------------------------------------------
