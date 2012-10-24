#include <stdlib.h>

#include "fastq_file.h"

/* ******************************************************
 *    		Function implementations  		*
 * ******************************************************/


fastq_file_t *fastq_fopen(char *filename) {
	return fastq_fopen_mode(filename, (char*)"r");
}

fastq_file_t *fastq_fopen_mode(char *filename, char *mode) {
	FILE *fd = fopen(filename, mode);
	char log_message[50];

	if (fd == NULL) {
		//        sprintf(log_message, "Error opening file: %s, mode (%s) !!!!!\n", filename, mode);
		LOG_FATAL_F("Error opening file: %s, mode (%s)\n", filename, mode);
		return NULL;
	}

	fastq_file_t* fq_file = (fastq_file_t*) malloc(sizeof(fastq_file_t));

	fq_file->filename = filename;
	fq_file->mode = mode;
	fq_file->fd = fd;

	return fq_file;
}



size_t fastq_fread_se(array_list_t *reads, size_t num_reads, fastq_file_t *fq_file) {
	size_t count = 0;
	char header1[MAX_READ_ID_LENGTH];
	char sequence[MAX_READ_SEQUENCE_LENGTH];
	char header2[MAX_READ_ID_LENGTH];
	char qualities[MAX_READ_SEQUENCE_LENGTH];
	int header_length, sequence_length, quality_length;
	fastq_read_t *read;

	while (count < num_reads && fgets(header1, MAX_READ_ID_LENGTH, fq_file->fd) != NULL) {
		fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
		fgets(header2, MAX_READ_ID_LENGTH, fq_file->fd);
		fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);

		header_length = strlen(header1);
		sequence_length = strlen(sequence);
		quality_length = strlen(qualities);

		// '\n' char is removed, but '\0' is left
		chomp_at(header1, header_length - 1);
		chomp_at(sequence, sequence_length - 1);
		chomp_at(qualities, quality_length - 1);

		read = fastq_read_new(header1, sequence, qualities);
		array_list_insert(read, reads);

		count++;
	}

	return count;
}

size_t fastq_fread_bytes_se(array_list_t *reads, size_t bytes, fastq_file_t *fq_file) {
	size_t accumulated_size = 0;
	char header1[MAX_READ_ID_LENGTH];
	char sequence[MAX_READ_SEQUENCE_LENGTH];
	char header2[MAX_READ_ID_LENGTH];
	char qualities[MAX_READ_SEQUENCE_LENGTH];
	int header_length, sequence_length, quality_length;
	fastq_read_t *read;

	while (accumulated_size < bytes && fgets(header1, MAX_READ_ID_LENGTH, fq_file->fd) != NULL) {
		fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
		fgets(header2, MAX_READ_ID_LENGTH, fq_file->fd);
		fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);

		header_length = strlen(header1);
		sequence_length = strlen(sequence);
		quality_length = strlen(qualities);

		// '\n' char is removed, but '\0' is left
		chomp_at(header1, header_length - 1);
		chomp_at(sequence, sequence_length - 1);
		chomp_at(qualities, quality_length - 1);

		read = fastq_read_new(header1, sequence, qualities);
		array_list_insert(read, reads);

		accumulated_size += header_length + sequence_length + quality_length;
	}

	return accumulated_size;
}



/*
 * OLD, BUT IN USE!!
 */


int fastq_fread(fastq_read_t *read, fastq_file_t *fq_file) {
	return fastq_fread_num_reads(read, 1, fq_file);
}

int fastq_fread_num_reads(fastq_read_t *buffer_fq_reads, int num_reads, fastq_file_t *fq_file) {
	int count = 0;
	char header1[MAX_READ_ID_LENGTH];
	char sequence[MAX_READ_SEQUENCE_LENGTH];
	char header2[MAX_READ_ID_LENGTH];
	char qualities[MAX_READ_SEQUENCE_LENGTH];
	int header_length, sequence_length, quality_length;

	while (count < num_reads && fgets(header1, MAX_READ_ID_LENGTH, fq_file->fd) != NULL) {
		fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
		fgets(header2, MAX_READ_ID_LENGTH, fq_file->fd);
		fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);

		header_length = strlen(header1);
		sequence_length = strlen(sequence);
		quality_length = strlen(qualities);

		chomp_at(header1, header_length - 1);
		chomp_at(sequence, sequence_length - 1);
		chomp_at(qualities, quality_length - 1);

		buffer_fq_reads[count].id = (char*)malloc(sizeof(char) * header_length);
		buffer_fq_reads[count].sequence = (char*)malloc(sizeof(char) * sequence_length);
		buffer_fq_reads[count].quality = (char*)malloc(sizeof(char) * quality_length);

		strcpy(buffer_fq_reads[count].id, header1);
		strcpy(buffer_fq_reads[count].sequence, sequence);
		strcpy(buffer_fq_reads[count].quality, qualities);

		count++;
	}

	return count;
}

int fastq_fread_max_size(fastq_read_t *buffer_fq_reads, unsigned long max_size, fastq_file_t *fq_file) {
	int count = 0;
	unsigned long accumulated_size = 0;
	char header1[MAX_READ_ID_LENGTH];
	char sequence[MAX_READ_SEQUENCE_LENGTH];
	char header2[MAX_READ_ID_LENGTH];
	char qualities[MAX_READ_SEQUENCE_LENGTH];
	int header_length, sequence_length, quality_length;

	while (accumulated_size <= max_size && fgets(header1, MAX_READ_ID_LENGTH, fq_file->fd) != NULL) {
		fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
		fgets(header2, MAX_READ_ID_LENGTH, fq_file->fd);
		fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);

		header_length = strlen(header1);
		sequence_length = strlen(sequence);
		quality_length = strlen(qualities);

		chomp_at(header1, header_length - 1);
		chomp_at(sequence, sequence_length - 1);
		chomp_at(qualities, quality_length - 1);

		buffer_fq_reads[count].id = (char*)malloc(sizeof(char) * header_length);
		buffer_fq_reads[count].sequence = (char*)malloc(sizeof(char) * sequence_length);
		buffer_fq_reads[count].quality = (char*)malloc(sizeof(char) * quality_length);

		strcpy(buffer_fq_reads[count].id, header1);
		strcpy(buffer_fq_reads[count].sequence, sequence);
		strcpy(buffer_fq_reads[count].quality, qualities);
		accumulated_size += header_length + sequence_length + quality_length;

		count++;
	}

	return count;
}

int fastq_fread_batch_max_size(fastq_batch_t *buffer_fq_read_batch, unsigned long max_size, fastq_file_t *fq_file) {
	unsigned long accumulated_size = 0;

	char header1[MAX_READ_ID_LENGTH];
	char sequence[MAX_READ_SEQUENCE_LENGTH];
	char header2[MAX_READ_ID_LENGTH];
	char qualities[MAX_READ_SEQUENCE_LENGTH];
	int header_length, sequence_length, quality_length;

	int count = 0;
	buffer_fq_read_batch->header_indices[count] = 0;
	buffer_fq_read_batch->data_indices[count] = 0;

	while (accumulated_size <= (max_size - 1024) && fgets(header1, MAX_READ_ID_LENGTH, fq_file->fd) != NULL) {
		fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
		fgets(header2, MAX_READ_ID_LENGTH, fq_file->fd);
		fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);

		header_length = strlen(header1);
		sequence_length = strlen(sequence);
		quality_length = strlen(qualities);

		if (sequence_length == quality_length) {
			// remove '\n' character, now length includes '\0' character
			chomp(header1);
			chomp(sequence);
			chomp(qualities);

			count++;

			strcpy(&(buffer_fq_read_batch->header[buffer_fq_read_batch->header_indices[count-1]]), header1);
			strcpy(&(buffer_fq_read_batch->seq[buffer_fq_read_batch->data_indices[count-1]]), sequence);
			strcpy(&(buffer_fq_read_batch->quality[buffer_fq_read_batch->data_indices[count-1]]), qualities);

			if (count*sizeof(int) >= buffer_fq_read_batch->data_indices_size) {

				// maybe realloc function can be used here
				int size = (count + 100) * sizeof(int);

				// copying data indices
				int* p = (int*) malloc(size);
				memset((void *) p, 0, size);
				memcpy((void*) p, (void*) buffer_fq_read_batch->data_indices, count * sizeof(int));

				free(buffer_fq_read_batch->data_indices);

				buffer_fq_read_batch->data_indices = p;

				// copying header indices
				p = (int*) malloc(size);
				memset((void *) p, 0, size);
				memcpy((void*) p, (void*) buffer_fq_read_batch->header_indices, count * sizeof(int));

				free(buffer_fq_read_batch->header_indices);

				buffer_fq_read_batch->header_indices = p;
				buffer_fq_read_batch->data_indices_size = size;
			}

			buffer_fq_read_batch->data_indices[count] = buffer_fq_read_batch->data_indices[count-1] + sequence_length;
			buffer_fq_read_batch->header_indices[count] = buffer_fq_read_batch->header_indices[count-1] + header_length;

			accumulated_size += sequence_length + quality_length;
		} else {
			LOG_DEBUG("Read has different length in sequence and quality");
		}

	}

	buffer_fq_read_batch->num_reads = count;

	return buffer_fq_read_batch->num_reads;
}

int fastq_fread_index_positions(fastq_read_t* buffer_reads, int *index_positions, fastq_file_t *fq_file) {
	const int max_length = 512;

	int count = 0;
	int index = 0;
	char header[max_length];
	char sequence[max_length];
	char plus[max_length];
	char quality[max_length];

	while (index_positions != NULL && fgets(header, max_length, fq_file->fd) != NULL) {
		fgets(sequence, max_length, fq_file->fd);
		fgets(plus, max_length, fq_file->fd);
		fgets(quality, max_length, fq_file->fd);

		if (count == index_positions[index]) {
			strcpy(buffer_reads[index].id, trim(header));
			strcpy(buffer_reads[index].sequence, trim(sequence));
			strcpy(buffer_reads[index].quality, trim(quality));

			index++;
		}

		count++;
	}

	return count;
}

int fastq_fwrite(fastq_read_t* buffer_reads, int num_writes, fastq_file_t *fq_file) {
	int count = 0;

	while (count < num_writes) {
		fprintf(fq_file->fd, "%s\n", buffer_reads->id);
		fprintf(fq_file->fd, "%s\n", buffer_reads->sequence);
		fprintf(fq_file->fd, "+\n");
		fprintf(fq_file->fd, "%s\n", buffer_reads->quality);

		buffer_reads++;
		count++;
	}

	return count;
}

unsigned int fastq_fcount(fastq_file_t *fq_file) {
	return fq_file->num_reads;
}

//-----------------------------------------------------
// fastq_remove_Ns
//-----------------------------------------------------
/*
void fastq_remove_Ns(fastq_read_t* buffer_reads, qc_read_t* qc_read, int max_N_per_read) {
    int count = 0;
    int index = 0;

    while (buffer_reads[count].id != NULL) {
        if (qc_read[count].counters[N] <= max_N_per_read) {
            buffer_reads[index] = buffer_reads[count];
            index++;
        }

        count++;
    }

    while (index <= count) {
        index++;
    }
}
 */
//-----------------------------------------------------
// fastq_fclose
//-----------------------------------------------------

void fastq_fclose(fastq_file_t* fq_file) {
	fclose(fq_file->fd);
	free(fq_file);
}

size_t consume_input(int c, char **data, size_t max_len, int position_in_data) {
	assert(data);

	(*data)[position_in_data] = c;
	// Text too long to be stored in 'data', realloc
	if (position_in_data == max_len - 1) {
		char *aux = realloc(*data, max_len + 10000);
		if (aux) {
			*data = aux;
			return max_len + 10000;
		} else {
			LOG_FATAL("Could not allocate enough memory for reading input VCF file\n");
		}
	}
	return max_len;
}
