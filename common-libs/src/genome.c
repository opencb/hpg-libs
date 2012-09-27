
#include "genome.h"

//-----------------------------------------------------

genome_t* genome_new(char* sequence_filename, char* chromosome_filename) {
  const int MAXLINE = 1024;
  genome_t* genome_p = (genome_t*) calloc(1, sizeof(genome_t));
  // genome file
  //
  FILE* fd = fopen(sequence_filename, "r");
  fseek(fd, 0L, SEEK_END);
  genome_p->nX = ftell(fd);
  fseek(fd, 0L, SEEK_SET);
  genome_p->X = (char*) calloc(1, genome_p->nX);
  fread(genome_p->X, 1, genome_p->nX, fd);
  fclose(fd);

  // read index file
  //
  unsigned int num_chromosomes = 0;
  unsigned int offset = 0;

  char* p;
  char line[MAXLINE];
  fd = fopen(chromosome_filename, "r");
  while (fgets(line, MAXLINE, fd) ) {
    p = strrchr(line, '\n'); *p = '\0';
    p = strrchr(line, '\t'); *p = '\0';
    //    printf("%s\n", line);
    strcpy((char*) genome_p->chr_name[num_chromosomes], line);
    genome_p->chr_name_length[num_chromosomes] = strlen(genome_p->chr_name[num_chromosomes]);
    sscanf(p+1, "%i", &genome_p->chr_size[num_chromosomes]);
    genome_p->chr_offset[num_chromosomes] = offset;

    offset += (genome_p->chr_size[num_chromosomes] + 1);
    num_chromosomes++;
  }
  fclose(fd);
  //  printf("In genome %d chromosome and %d\n", num_chromosomes, offset);
  genome_p->num_chromosomes = num_chromosomes;

  return genome_p;
}

//-----------------------------------------------------

void genome_free(genome_t* genome_p) {
  if (genome_p == NULL) {
    return;
  }

  if (genome_p->X != NULL) free(genome_p->X);
  free(genome_p);
}

//-----------------------------------------------------

void genome_read_sequence(char* sequence, unsigned int strand, char* chromosome, unsigned long int* start_p, unsigned long int* end_p, genome_t* genome_p) {
  unsigned int i, chr;

  for(i=0 ; i< genome_p->num_chromosomes ; i++) {
    if (strcmp(chromosome, (char *) genome_p->chr_name[i]) == 0) {
      chr = i;
      break;
    }
  }

  genome_read_sequence_by_chr_index(sequence, strand, chr, start_p, end_p, genome_p);
}

//-----------------------------------------------------

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

//------------------------------------------------------------------------------------

char* genome_get_chr_name(unsigned int chr, unsigned int* len, genome_t* genome_p) {
  *len = genome_p->chr_name_length[chr];
  return genome_p->chr_name[chr];
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
