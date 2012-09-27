
#include "BW_search.h"
#include "buffers.h"
#include "fastq_ex_batch.h"
#include "bw_server.h"
#include "commons.h"
#include <omp.h>
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

/*
extern unsigned int *h_C, *h_O, *h_I;
extern unsigned int nC, nO, sizO, nI, mI, mO;

extern unsigned int *S;
extern unsigned int nS;
*/



void bwt_server_omp(bw_server_input_omp_t* input_p, coi_t* coi_data, coi_t* coi_data_rev)
{
	  
	  list_item_t *item_p = NULL;
	  int read;
	  int batchMaps=0;
	  int length=0;
	  
	  list_item_t* kl_item;
	  list_item_t* kl_item_rev;
	  
	  list_t* kl_list_p = input_p->kl_list_p;
	  kl_t *kl_p;
	  
	  int bwt_id = omp_get_thread_num();
	  
	  if(input_p->mode==NORMAL_MODE)
	  {
	      fastq_ex_batch_t *batch;
	      printf("bwt_search_omp:START\n");
	      while ( (item_p=list_remove_item(input_p->read_list_p)) != NULL ) {
		 // printf(" BW extract batch!\n");
		  batch=(fastq_ex_batch_t *)item_p->data_p;
		  
		  kl_p = kl_new(NORMAL_MODE, (void*) batch);
		  kl_p->k_p = (unsigned int*) malloc((unsigned int) batch->num_reads * 2 * sizeof(unsigned int));
		  kl_p->l_p = (unsigned int*) malloc((unsigned int) batch->num_reads * 2 * sizeof(unsigned int));
		  

		  #pragma omp parallel for private(length) shared(kl_p) schedule(dynamic, MAX(1, batch->num_reads/( omp_get_num_threads()*10 ))) 
		  for(read=0; read< batch->num_reads;  read++)
		  {
		    length=batch->data_indices[read+1]-batch->data_indices[read]-1;
		    
		    BWExactSearchBackwardCompNew(&(batch->seq[batch->data_indices[read]]),    0,  length-1, 0,     coi_data->mO-2,     coi_data->_C,     coi_data->_O,     coi_data->sizO,     coi_data->_I,     coi_data->mI, read, kl_p);
		    BWExactSearchBackwardCompRevNew(&(batch->seq[batch->data_indices[read]]), 0,  length-1, 0, coi_data_rev->mO-2, coi_data_rev->_C, coi_data_rev->_O, coi_data_rev->sizO, coi_data_rev->_I, coi_data_rev->mI, read+batch->num_reads, kl_p);
		    
		  }
		  
		  kl_item=list_item_new(batch->source_id, KL_ITEM, kl_p);
		  list_insert_item(kl_item, kl_list_p);
		  
		  list_item_free(item_p);
		//printf(" BW extract batch Finish!\n"); 
	      }
	      list_decr_writers(kl_list_p);
	      printf("bwt_search_omp:FINISH\n");
	    // printf("Total maps %lld. Total reads %lld => %d\%\n", nMaps, nReads, (100*nMaps)/nReads);
	  }else
	  {
	      split_batch_t *split_batch_p;
	      int split_size;
	      int splits_per_read;
	      int num_splits;
	      int sp;
	      char search[1024], split[1024];
	      int len;
	      int totsplits=0;
	      printf("bwt_split_search_omp:START\n");
	      while ( (item_p=list_remove_item(input_p->read_list_p)) != NULL ) {
		  //printf(" BW split extract batch!\n");
		  split_batch_p=(split_batch_t *)item_p->data_p;
		  split_size = split_batch_p->split_size;
		  num_splits=split_batch_p->num_splits;
		  
		  printf("Total Reads = %d\n", split_batch_p->num_reads);
		  
		  kl_p = kl_new(SPLIT_MODE, (void*) split_batch_p);
		  kl_p->k_p = (unsigned int*) malloc((unsigned int) num_splits * 2 * sizeof(unsigned int));
		  kl_p->l_p = (unsigned int*) malloc((unsigned int) num_splits * 2 * sizeof(unsigned int));
		  
		  #pragma omp parallel for private(length, sp) shared(kl_p) schedule(dynamic, MAX(1, split_batch_p->num_reads/( omp_get_num_threads()*10 ))) 
		  for(read=0; read< split_batch_p->num_reads;  read++)
		  {
		      splits_per_read=split_batch_p->splits_per_read;
		      totsplits+=splits_per_read;
		      
		      len = split_batch_p->read_indices_p[read+1] - split_batch_p->read_indices_p[read] - 1 ;
		      decodeBases(search, &(split_batch_p->read_p[split_batch_p->read_indices_p[read]]), len);
		      
		     // printf("BWT Search Split read=%d, len=%d, search=%s\n", read,  len, search);
		      for (sp=0 ; sp<splits_per_read ; sp++) {
			    
			    //printf("Split search pos=%d, indice=%d\n", sp + (read*splits_per_read), split_batch_p->split_indices_p[sp + (read*splits_per_read)] );
			    BWExactSearchBackwardCompNew(&(split_batch_p->read_p[split_batch_p->split_indices_p[sp + (read*splits_per_read)]]),    0,  split_size-1, 0,     coi_data->mO-2,     coi_data->_C,     coi_data->_O,     coi_data->sizO,     coi_data->_I,     coi_data->mI, sp+(read*splits_per_read) , kl_p);
			    BWExactSearchBackwardCompRevNew(&(split_batch_p->read_p[split_batch_p->split_indices_p[sp + (read*splits_per_read)]]), 0,  split_size-1, 0, coi_data_rev->mO-2, coi_data_rev->_C, coi_data_rev->_O, coi_data_rev->sizO, coi_data_rev->_I, coi_data_rev->mI, sp+(read*splits_per_read)+(num_splits) , kl_p);
			    
			    //decodeBases(split, &(split_batch_p->read_p[split_batch_p->split_indices_p[sp + (read*splits_per_read)]]), split_size);
			    //printf("BWT Split Search read %s\n", split);
			   // printf("Search data read=%d +pos=%d, -pos=%d, +(k=%d, l=%d)--(k=%d, l=%d)\n", read, sp+(read*splits_per_read), sp+(read*splits_per_read)+(num_splits), kl_p->k_p[sp+(read*splits_per_read)], kl_p->l_p[sp+(read*splits_per_read)], kl_p->k_p[sp+(read*splits_per_read)+(num_splits)], kl_p->l_p[sp+(read*splits_per_read)+(num_splits)]);
			
		      }
		  }
		 // printf(" BW split extract batch Finish!\n"); 
		  kl_item=list_item_new(item_p->id, KL_ITEM, kl_p);
		  list_insert_item(kl_item, kl_list_p);  
		  
		  list_item_free(item_p);
	    }
	    list_decr_writers(kl_list_p);
	    printf("bw_split_search_omp:FINISH\n");
	}
	
	// printf("Total maps %lld. Total reads %lld => %d\%\n", nMaps, nReads, (100*nMaps)/nReads);
}



void server_read_kl(list_t* kl_list_p)
{
  list_item_t* item_p;
  fastq_ex_batch_t* batch;
  char *search;
  kl_t *kl_p;
  int read;
  int strand;
  int mapsPlus=0, mapsMinus=0;
  int len;
  while ( (item_p=list_remove_item(kl_list_p)) != NULL ) {
 
      mapsPlus=0; mapsMinus=0;
      kl_p=(kl_t *)item_p->data_p;
      batch=(fastq_ex_batch_t*)kl_p->batch_p;
      
      
      for(read=0; read< batch->num_reads;  read++)
      {
	len = batch->data_indices[read+1] - batch->data_indices[read] - 1;
	decodeBases(search, &(batch->seq[batch->data_indices[read]]), len);
	printf("Split read: %s\n", search);
	for(strand=0; strand<2; strand++)
	  if(!( kl_p->k_p[read+(batch->num_reads*strand)] > kl_p->l_p[read+(batch->num_reads*strand)]))
	  {
	    if(strand==0)
	      mapsPlus++;
	    else
	      mapsMinus++;
	    //printf("Read %d: strand=%d k=%d, l=%d\n", read, strand, kl_p->k_p[read], kl_p->l_p[read]);
	  }/*else
	  {
	    printf("Read %d: No mapped\n", read);
	  }*/
      }
      kl_free(kl_p);
     // printf("Process Batch. Match strand + %d and Match strand - %d\n", mapsPlus, mapsMinus);
  }
  
}


void server_read_split(list_t* split_list_p)
{
  list_item_t* item_p;
  split_batch_t* split_batch_p;
  int split_size, splits_per_read;
  char search[1024], split[1024];
  int read_len;
  int read;
  int strand;
  int len;
  int sp;
  while ( (item_p=list_remove_item(split_list_p)) != NULL ) {
      split_batch_p=(split_batch_t *)item_p->data_p;
      printf("Extract split read...\n");
      
      split_size = split_batch_p->split_size;
      splits_per_read = split_batch_p->splits_per_read;
      
      for(read=0; read< split_batch_p->num_reads;  read++)
      {
	read_len = split_batch_p->read_indices_p[read+1] - split_batch_p->read_indices_p[read];
	decodeBases(search, &(split_batch_p->read_p[split_batch_p->read_indices_p[read]]), read_len);
	
	printf("Read No mapped %s:\n split in %d and split size %d\n", search, splits_per_read, split_size);
	  for (sp=0 ; sp<splits_per_read ; sp++) {
	      decodeBases(split, &(split_batch_p->read_p[split_batch_p->split_indices_p[sp + (read*splits_per_read)]]), split_size);
	      printf("Split %s\n", split);
	  }
	
      }
      split_batch_free(split_batch_p);
     // printf("Process Batch. Match strand + %d and Match strand - %d\n", mapsPlus, mapsMinus);
  }
  
}

void server_read_kl_split(list_t* kl_list_p)
{
  list_item_t* item_p;
  split_batch_t* split_batch_p;
  char search[1024], split[1024];
  kl_t *kl_p;
  int read;
  int len;
  int num_reads, num_splits, split_size, splits_per_read, read_len;
  int const STRAND=2;
  int strand, sp;
  int k_aux;
  int l_aux;

  while ( (item_p=list_remove_item(kl_list_p)) != NULL ) {
      printf("Extract one element to list. Total %d\n", kl_list_p->length);
      kl_p = (kl_t*) item_p->data_p;
      split_batch_p = (split_batch_t*) kl_p->batch_p;
      num_reads = split_batch_p->num_reads;
      num_splits = split_batch_p->num_splits;
  
      split_size = split_batch_p->split_size;
      splits_per_read = split_batch_p->splits_per_read;
      for (read=0; read < num_reads; read++) {
      
	read_len = split_batch_p->read_indices_p[read+1] - split_batch_p->read_indices_p[read];
	decodeBases(search, &(split_batch_p->read_p[split_batch_p->read_indices_p[read]]), read_len);
	
	for(strand=0 ; strand<STRAND ; strand++) {
	  printf("Read %d strand %d: %s\n", read, strand, search);

	  for (sp=0 ; sp<splits_per_read ; sp++) {
	    decodeBases(split, &(split_batch_p->read_p[split_batch_p->split_indices_p[sp + (read*splits_per_read)]]), split_size);
	    k_aux = kl_p->k_p[sp + (read*splits_per_read) + (num_splits * strand)];
	    l_aux = kl_p->l_p[sp + (read*splits_per_read) + (num_splits * strand)];
	    printf("  Split %d: %s\n", sp, split);
	    if(k_aux > l_aux)
	      printf("     K=%d, l=%d \n", k_aux, l_aux);
	    else
	      printf("     No Map\n");
	  }
	}
   }
   kl_free(kl_p);
     // printf("Process Batch. Match strand + %d and Match strand - %d\n", mapsPlus, mapsMinus);
  }
  
}


