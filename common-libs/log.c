#include <stdio.h>
#include <time.h>
#include <string.h>

#include "log.h"
#include "string_utils.h"

int log_level;
int verbose;
char *log_filename;

void set_log_verbose(int v) {
	verbose = v;
}

void set_log_level(int level) {
	log_level = level;
}

void set_log_filename(char *filename) {
	log_filename = (char*) calloc(1024, sizeof(char));
	strcpy(log_filename, filename);
}

void print_log_message(int level, char *log_level_word, char *filename, int num_line, const char *func, char *msg) {
	if((level >= log_level) && (verbose || strcmp(log_filename, "") != 0)) {
		time_t rawtime;
		time(&rawtime);
		char* str_time = ctime(&rawtime);
		chomp(str_time);

		// if 'verbose' logs are printed in stdout
		if(verbose) {
			printf("%s\t%s\t%s [%i] in %s(): %s\n", str_time, log_level_word,
			filename, num_line, func, msg);
		}

		// if 'log_file' has been set up then logs are printed
		// logs are ALWAYS printed in log_file independently of 'verbose'		
		if(log_filename != NULL) {
			FILE *log_file = fopen(log_filename, "a");
			fprintf(log_file, "%s\t%s\t%s [%i] in %s(): %s\n", str_time, log_level_word, filename, num_line, func, msg);
			fclose(log_file);
		}

	}
}
