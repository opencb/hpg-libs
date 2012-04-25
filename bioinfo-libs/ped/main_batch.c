#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include <log.h>

#include "ped_batch.h"
#include "ped_file_structure.h"
#include "ped_file.h"
#include "ped_read.h"
#include "ped_reader.h"
#include "ped_write.h"

#include "list.h"

int main (int argc, char *argv[])
{
    size_t max_batches = 20;
    size_t batch_size = 2000;
    list_t *read_list = (list_t*) malloc (sizeof(list_t));
    list_init("batches", 1, max_batches, read_list);

    int ret_code;
    double start, stop, total;
    char *filename = (char*) malloc ((strlen(argv[1])+1) * sizeof(char));
    strncat(filename, argv[1], strlen(argv[1]));
    ped_file_t* file;
    
    init_log_custom(1, 1, NULL);
    
#pragma omp parallel sections private(start, stop, total) lastprivate(file)
{
    #pragma omp section
    {
        LOG_DEBUG_F("Thread %d reads the PED file\n", omp_get_thread_num());
        // Reading
        start = omp_get_wtime();
        
        file = ped_open(filename);
        ret_code = ped_read_batches(read_list, batch_size, file);
        
        stop = omp_get_wtime();
        total = (stop - start);
        
        if (ret_code) {
            LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code);
        }
        LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
        LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        
        // Writing to a new file
        if (argc == 3) 
        {
            start = omp_get_wtime();
        
            ret_code = ped_write(file, argv[2]);
            
            stop = omp_get_wtime();
            total = (stop - start);
            
            if (ret_code) {
                LOG_ERROR_F("[%dW] Error code = %d\n", omp_get_thread_num(), ret_code);
            }
            LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        }
        
        list_decr_writers(read_list);
        
        ped_close(file, 0);
    }
    #pragma omp section
    {
        printf("1st log debug\n");
        LOG_DEBUG_F("OMP num threads = %d\n", omp_get_num_threads());
        LOG_DEBUG_F("Thread %d prints info\n", omp_get_thread_num());
        printf("after 1st log debug\n");
        
        start = omp_get_wtime();
        
        int i = 0;
        list_item_t* item = NULL;
        FILE *out = fopen("result.ped", "w");
        while ( (item = list_remove_item(read_list)) != NULL ) {
            if (i % 200 == 0) 
            {
                int debug = 1;
                LOG_DEBUG_F("Batch %d reached by thread %d - %zu/%zu records \n", i, omp_get_thread_num(), 
                    ((ped_batch_t*) item->data_p)->length, ((ped_batch_t*) item->data_p)->max_length);
            }
            
//             ped_write_to_file(file, out);
//             ped_batch_print(stdout, item->data_p);
            write_ped_batch(item->data_p, out);
            ped_batch_free(item->data_p);
            list_item_free(item);
            i++;
        }
        fclose(out);
        
        stop = omp_get_wtime();
        total = (stop - start);
        
        LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
        LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
    }
}

    free(read_list);

    return 0;
}
