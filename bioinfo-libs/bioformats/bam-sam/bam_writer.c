
#include "bam_writer.h"

/* ******************************************************
 *      	Function implementations    		*
 * *****************************************************/

bam_writer_t* bam_writer_new(char* filename, alignments_list_t* alignments_list_p, bam_header_t* bam_header_p, int mode, int chromosome) {
    bam_writer_t* writer_p = (bam_writer_t*) calloc(1, sizeof(bam_writer_t));

    // open the output file....
    writer_p->bam_file_p = bam_fopen_mode(filename, bam_header_p, "w");
    writer_p->mode = mode;
    writer_p->chromosome = chromosome;

    writer_p->alive = 1;
    pthread_mutex_init(&(writer_p->alive_lock), NULL);

    writer_p->alignments_list_p = alignments_list_p;
    writer_p->bam_write_alignment_count = (int*) calloc(num_of_chromosomes, sizeof(int));

    bam_fwrite_header(bam_header_p, writer_p->bam_file_p);

    return writer_p;
}

void bam_writer_free(bam_writer_t* writer_p, int all) {
    if (writer_p == NULL) {
        return;
    }

    // close the output file and exiting....
    if (writer_p->bam_file_p != NULL) {
        bam_fclose(writer_p->bam_file_p);
        writer_p->bam_file_p = NULL;
    }
    
    if (all) {
        if (writer_p->alignments_list_p != NULL) {
            free(writer_p->alignments_list_p);
            writer_p->alignments_list_p = NULL;
        }
    }
    
    if (writer_p->bam_write_alignment_count != NULL) {
        free(writer_p->bam_write_alignment_count);
        writer_p->bam_write_alignment_count = NULL;
    }
    
    free(writer_p);
    writer_p = NULL;
}

void bam_writer_start(bam_writer_t* writer_p) {
    // create and launch pthread....
    LOG_DEBUG("Launching bam write server thread....\n");

    if (writer_p->mode == CHROMOSOME_MODE) {
        pthread_create(&(writer_p->thread), NULL, bam_writer_chromosome_thread_function, (void*) writer_p);
    } else if (writer_p->mode == SEQUENTIAL_MODE) {
        pthread_create(&(writer_p->thread), NULL, bam_writer_sequential_thread_function, (void*) writer_p);
    }
}

unsigned int bam_writer_join(bam_writer_t* writer_p) {
    void* r;
    pthread_join(writer_p->thread, &r);
    return (uintptr_t) r;
}

void bam_writer_set_alive(bam_writer_t* writer_p, int alive) {
    pthread_mutex_lock(&(writer_p->alive_lock));
    writer_p->alive = alive;
    pthread_mutex_unlock(&(writer_p->alive_lock));
}

int bam_writer_get_alive(bam_writer_t* writer_p) {
    int alive;
    pthread_mutex_lock(&(writer_p->alive_lock));
    alive = writer_p->alive;
    pthread_mutex_unlock(&(writer_p->alive_lock));
    return alive;
}

/* **************************************************************
 *      	Thread functions implementations    		*
 * *************************************************************/

void* bam_writer_chromosome_thread_function(void* param_p) {
    bam_writer_t* writer_p = (bam_writer_t*) param_p;
    bam1_t* alignment_p;

    bam_writer_set_alive(writer_p, 1);

    char log_message[200];
    sprintf(log_message, "Thread-WRITE: START, for file '%.150s' and chromosome %i\n", writer_p->bam_file_p->filename, writer_p->chromosome);
    LOG_DEBUG(log_message);

    int current_chromosome = writer_p->chromosome;
    int num_alignments = 0;
    int write_bytes = 0, total_write_bytes = 0;

    chrom_alignments_t* chrom_alignments_p;

    while ((chrom_alignments_p = alignments_list_get_chrom_alignment(current_chromosome, writer_p->alignments_list_p)) == NULL) {
        sched_yield();
        usleep(10000);
    }

    while (1) {
        if ((chrom_alignments_is_complete(chrom_alignments_p)) && (alignment_p = chrom_alignments_get_alignment(chrom_alignments_p, num_alignments)) != NULL) {
            if (time_flag) {
                start_timer(t1_write);
            }
            write_bytes = bam_fwrite(alignment_p, writer_p->bam_file_p);
            total_write_bytes += write_bytes;
            num_alignments++;
            if (time_flag) {
                stop_timer(t1_write, t2_write, write_time);
            }

            if ((chrom_alignments_is_complete(chrom_alignments_p)) && (num_alignments == chrom_alignments_p->alignment_count)) {
                break;
            } else if (!chrom_alignments_is_complete(chrom_alignments_p)) {
                sched_yield();
                usleep(10000);
            }
        } else {
            if ((chrom_alignments_is_complete(chrom_alignments_p)) && (num_alignments == chrom_alignments_p->alignment_count)) {
                break;
            } else {
                sched_yield();
                usleep(10000);
            }
        }
    }

    bam_fclose(writer_p->bam_file_p);
    writer_p->bam_file_p = NULL;

    bam_writer_set_alive(writer_p, 0);

    sprintf(log_message, "Thread-WRITE: END for chromosome %i \n", current_chromosome);
    LOG_DEBUG(log_message);

    pthread_exit((void*)num_alignments);
}

void* bam_writer_sequential_thread_function(void* param_p) {
    struct timespec ts;
    ts.tv_sec = 0;
    ts.tv_nsec = 1000000;

    bam_writer_t* writer_p = (bam_writer_t*) param_p;

    bam_writer_set_alive(writer_p, 1);

    char log_message[200];
    sprintf(log_message, "Thread-WRITE: START, for file '%.150s' and chromosome %i\n", writer_p->bam_file_p->filename, writer_p->chromosome);
    LOG_DEBUG(log_message);

    int current_chromosome = 0;
    int num_alignments = 0;
    int total_write_bytes = 0;

    chrom_alignments_t* chrom_alignments_p;

    while ((chrom_alignments_p = alignments_list_get_chrom_alignment(current_chromosome, writer_p->alignments_list_p)) == NULL) {
        sched_yield();
        usleep(10000);
    }

    while (1) {
        if ((current_chromosome < num_of_chromosomes) && (chrom_alignments_is_complete(chrom_alignments_p))) {

            if (time_flag) {
                start_timer(t1_write);
            }

            total_write_bytes = bam_fwrite_sorted_array(chrom_alignments_p->bam_alignments_p, chrom_alignments_p->indices_p, chrom_alignments_p->alignment_count, writer_p->bam_file_p);
            num_alignments += chrom_alignments_p->alignment_count;

            if (time_flag) {
                stop_timer(t1_write, t2_write, write_time);
            }

            current_chromosome++;

            while (((chrom_alignments_p = alignments_list_get_chrom_alignment(current_chromosome, writer_p->alignments_list_p)) == NULL) && (current_chromosome < num_of_chromosomes)) {
                sched_yield();
                usleep(10000);
            }
        } else if (current_chromosome == num_of_chromosomes) {  
            break;
        } else {
            sched_yield();
            usleep(10000);
        }
    }

    bam_fclose(writer_p->bam_file_p);
    writer_p->bam_file_p = NULL;

    bam_writer_set_alive(writer_p, 0);

    sprintf(log_message, "Thread-WRITE: END with total alignments: %i \n", num_alignments);
    LOG_DEBUG(log_message);

    pthread_exit((void*)num_alignments);
}
