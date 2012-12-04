#include "bam_data_batch.h"

/* ******************************************************
 *      	Private functions 	    		*
 * *****************************************************/

void bam_data_batch_init_cigar_(bam_data_batch_t* bam_data_batch_p, bam_batch_t* bam_batch_p, int cigar_operations);
void bam_data_batch_init_id_seq_(bam_data_batch_t* bam_data_batch_p, bam_batch_t* bam_batch_p, int id_seq_length);

/* **************************************************************
 *      	Public functions implementations	  	*
 * *************************************************************/

bam_data_batch_t* bam_data_batch_new(int num_alignments) {
    bam_data_batch_t* bam_data_batch_p = (bam_data_batch_t*) calloc(1, sizeof(bam_data_batch_t));
    bam_data_batch_p->num_alignments = num_alignments;
    bam_data_batch_p->core_data_p = (bam_data_core_t*) calloc((num_alignments + 1), sizeof(bam_data_core_t));
    bam_data_batch_p->cigar_data_p = NULL;
    bam_data_batch_p->id_seq_data_p = NULL;

    return bam_data_batch_p;
}

bam_data_batch_t* bam_data_batch_init(bam_data_batch_t* bam_data_batch_p, bam_batch_t* bam_batch_p) {
    int cigar_operations = 0, id_seq_length = 0, last_alignments_position = 0, first_chromosome_position = 0;
    short int num_chromosomes_in_batch = 0;

    short int chromosomes[num_of_chromosomes];
    int start_positions[num_of_chromosomes];
    int end_positions[num_of_chromosomes];

    bam_data_core_t* core_data_p;

    for (int i = 0; i < bam_batch_p->num_alignments; i++) {
        core_data_p = &(bam_data_batch_p->core_data_p[i]);

        core_data_p->strand = bam1_strand(bam_batch_p->alignments_p[i]);
        core_data_p->map_quality = bam_batch_p->alignments_p[i]->core.qual;
        core_data_p->chromosome = bam_batch_p->alignments_p[i]->core.tid;
        core_data_p->mate_chromosome = bam_batch_p->alignments_p[i]->core.mtid;
        core_data_p->isize = abs(bam_batch_p->alignments_p[i]->core.isize);
        core_data_p->start_coordinate = bam_batch_p->alignments_p[i]->core.pos;
        core_data_p->alignment_length = bam_cigar2qlen(&(bam_batch_p->alignments_p[i]->core), bam1_cigar(bam_batch_p->alignments_p[i]));

        cigar_operations += bam_batch_p->alignments_p[i]->core.n_cigar;
        id_seq_length += bam_batch_p->alignments_p[i]->core.l_qname;

        if (bam_batch_p->alignments_p[i]->core.flag & BAM_FREAD1) {
            core_data_p->paired_end = PAIRED_END1;
        } else if (bam_batch_p->alignments_p[i]->core.flag & BAM_FREAD2) {
            core_data_p->paired_end = PAIRED_END2;
        } else {
            core_data_p->paired_end = NO_PAIRED_END;
        }

        if ((last_alignments_position != 0) && (bam_batch_p->alignments_p[i]->core.tid == bam_batch_p->alignments_p[i-1]->core.tid)) {
            last_alignments_position = max(last_alignments_position, bam_batch_p->alignments_p[i]->core.pos + bam_batch_p->alignments_p[i]->core.l_qseq);
        } else {
            first_chromosome_position = bam_batch_p->alignments_p[i]->core.pos;

            chromosomes[num_chromosomes_in_batch] = bam_batch_p->alignments_p[i]->core.tid;

            start_positions[num_chromosomes_in_batch] = first_chromosome_position;
            if (last_alignments_position != 0) {
                end_positions[num_chromosomes_in_batch - 1] = last_alignments_position;
            }

            num_chromosomes_in_batch++;
            last_alignments_position = bam_batch_p->alignments_p[i]->core.pos + bam_batch_p->alignments_p[i]->core.l_qseq;
        }
    }

    // count for the last chromosome
    end_positions[num_chromosomes_in_batch - 1] = last_alignments_position;

    bam_data_batch_p->num_cigar_operations = cigar_operations;
    bam_data_batch_p->id_seq_length = id_seq_length;
    bam_data_batch_p->last_alignments_position = last_alignments_position;
    bam_data_batch_p->num_chromosomes_in_batch = num_chromosomes_in_batch;

    bam_data_batch_p->chromosomes = (short int*) calloc(num_chromosomes_in_batch, sizeof(short int));
    bam_data_batch_p->start_positions = (int*) calloc(num_chromosomes_in_batch, sizeof(int));
    bam_data_batch_p->end_positions = (int*) calloc(num_chromosomes_in_batch, sizeof(int));

    memcpy(bam_data_batch_p->chromosomes, chromosomes, num_chromosomes_in_batch * sizeof(short int));
    memcpy(bam_data_batch_p->start_positions, start_positions, num_chromosomes_in_batch * sizeof(int));
    memcpy(bam_data_batch_p->end_positions, end_positions, num_chromosomes_in_batch * sizeof(int));

    bam_data_batch_init_cigar_(bam_data_batch_p, bam_batch_p, cigar_operations);
    bam_data_batch_init_id_seq_(bam_data_batch_p, bam_batch_p, id_seq_length);

    return bam_data_batch_p;
}

bam_data_batch_t* bam_data_batch_by_chromosome_init(bam_data_batch_t* bam_data_batch_p, bam_batch_t* bam_batch_p, int start_alignment) {
    int cigar_operations = 0, id_seq_length = 0, bam_data_num_alignments = 0;

    for (int i = start_alignment; i < bam_batch_p->num_alignments; i++) {
        if ((i > start_alignment) && (bam_data_batch_p->core_data_p[i-1].chromosome != bam_batch_p->alignments_p[i]->core.tid)) {
	    break;
	}

        bam_data_batch_p->core_data_p[i].strand = bam1_strand(bam_batch_p->alignments_p[i]);
        bam_data_batch_p->core_data_p[i].map_quality = bam_batch_p->alignments_p[i]->core.qual;
        bam_data_batch_p->core_data_p[i].chromosome = bam_batch_p->alignments_p[i]->core.tid;
        bam_data_batch_p->core_data_p[i].start_coordinate = bam_batch_p->alignments_p[i]->core.pos;
        bam_data_batch_p->core_data_p[i].alignment_length = bam_batch_p->alignments_p[i]->core.l_qseq;
        cigar_operations += bam_batch_p->alignments_p[i]->core.n_cigar;
        id_seq_length += bam_batch_p->alignments_p[i]->core.l_qname;

        if (bam_batch_p->alignments_p[i]->core.flag & BAM_FREAD1) {
            bam_data_batch_p->core_data_p[i].paired_end = PAIRED_END1;
        } else if (bam_batch_p->alignments_p[i]->core.flag & BAM_FREAD2) {
            bam_data_batch_p->core_data_p[i].paired_end = PAIRED_END2;
        } else {
            bam_data_batch_p->core_data_p[i].paired_end = NO_PAIRED_END;
        }

        bam_data_num_alignments++;
    }

    bam_data_batch_p->num_alignments = bam_data_num_alignments;
    bam_data_batch_p->num_cigar_operations = cigar_operations;
    bam_data_batch_p->id_seq_length = id_seq_length;

    bam_data_batch_init_cigar_(bam_data_batch_p, bam_batch_p, cigar_operations);
    bam_data_batch_init_id_seq_(bam_data_batch_p, bam_batch_p, id_seq_length);

    return bam_data_batch_p;
}

void bam_data_batch_free(bam_data_batch_t* bam_data_batch_p) {
    if (bam_data_batch_p->core_data_p != NULL) {
        free(bam_data_batch_p->core_data_p);
    }
    if (bam_data_batch_p->cigar_data_p != NULL) {
        free(bam_data_batch_p->cigar_data_p);
    }
    if (bam_data_batch_p->id_seq_data_p != NULL) {
        free(bam_data_batch_p->id_seq_data_p);
    }

    if (bam_data_batch_p->chromosomes != NULL) {
        free(bam_data_batch_p->chromosomes);
    }
    if (bam_data_batch_p->start_positions != NULL) {
        free(bam_data_batch_p->start_positions);
    }
    if (bam_data_batch_p->end_positions != NULL) {
        free(bam_data_batch_p->end_positions);
    }

    free(bam_data_batch_p);
}

bam_data_core_t* bam_data_batch_get_core(bam_data_batch_t* bam_data_batch_p) {
    return bam_data_batch_p->core_data_p;
}

uint32_t* bam_data_batch_get_cigar_data(bam_data_batch_t* bam_data_batch_p) {
    return bam_data_batch_p->cigar_data_p;
}

char* bam_data_batch_get_id_seq_data(bam_data_batch_t* bam_data_batch_p) {
    return bam_data_batch_p->id_seq_data_p;
}

/* **************************************************************
 *      	Private functions implementations	  	*
 * *************************************************************/

void bam_data_batch_init_cigar_(bam_data_batch_t* bam_data_batch_p, bam_batch_t* bam_batch_p, int cigar_operations) {
    int index = 0;
    uint32_t* cigar_alignment_p;

    bam_data_batch_p->cigar_data_p = (uint32_t*) calloc(cigar_operations, sizeof(uint32_t));
    bam_data_batch_p->core_data_p[0].cigar_index = index;

    for (int i = 0; i < bam_batch_p->num_alignments; i++) {
        cigar_alignment_p = bam1_cigar(bam_batch_p->alignments_p[i]);

        if (cigar_alignment_p != NULL) {
            memcpy((void *) &(bam_data_batch_p->cigar_data_p[index]), (void *) cigar_alignment_p, bam_batch_p->alignments_p[i]->core.n_cigar * sizeof(uint32_t));
        } else {
            LOG_DEBUG("cigar is NULL!!!\n");
        }

        index += bam_batch_p->alignments_p[i]->core.n_cigar;
        bam_data_batch_p->core_data_p[i+1].cigar_index = index;
    }

    bam_data_batch_p->core_data_p[bam_batch_p->num_alignments].cigar_index = cigar_operations;
}

void bam_data_batch_init_id_seq_(bam_data_batch_t* bam_data_batch_p, bam_batch_t* bam_batch_p, int id_seq_length) {
    int index = 0;
    char* id_seq_alignment_p;

    bam_data_batch_p->id_seq_data_p = (char*) calloc(id_seq_length, sizeof(char));
    bam_data_batch_p->core_data_p[0].id_seq_index = index;

    for (int i = 0; i < bam_batch_p->num_alignments; i++) {
        id_seq_alignment_p = bam1_qname(bam_batch_p->alignments_p[i]);

        if (id_seq_alignment_p != NULL) {
            memcpy((void *) &(bam_data_batch_p->id_seq_data_p[index]), (void *) id_seq_alignment_p, bam_batch_p->alignments_p[i]->core.l_qname * sizeof(char));
        } else {
	    LOG_DEBUG("id seq is NULL!!!\n");
        }

        index += bam_batch_p->alignments_p[i]->core.l_qname;
        bam_data_batch_p->core_data_p[i+1].id_seq_index = index;
    }

    bam_data_batch_p->core_data_p[bam_batch_p->num_alignments].id_seq_index = id_seq_length;
}
