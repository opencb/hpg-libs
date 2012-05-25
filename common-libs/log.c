#include "log.h"

/* **********************************************
 *  		Global variables		*
 * *********************************************/

int log_level;
int verbose;
char *log_filename;


void init_log() {
    set_log_verbose(0);
    set_log_level(LOG_INFO_LEVEL);
    set_log_filename("output.log");
}

void init_log_custom(int level, int verbose, char* log_filename) {
    set_log_level(level);
    set_log_verbose(verbose);
    set_log_filename(log_filename);
}

void set_log_verbose(int v) {
    verbose = v;
}

void set_log_level(int level) {
    log_level = level;
}

void set_log_filename(char *filename) {
    if (filename != NULL && strlen(filename) > 0) {
        log_filename = (char*) calloc (strlen(filename) + 1, sizeof(char));
        strncat(log_filename, filename, strlen(filename));
    }
}

void print_log_message(int level, char *log_level_word, char *filename, int num_line, const char *func, char *msg) {
    if((level >= log_level) && (verbose || log_filename != NULL)) {
        time_t rawtime;
        time(&rawtime);
        char* str_time = ctime(&rawtime);
        chomp(str_time);

        // if 'verbose' logs are printed in stderr
        if(verbose) {
            fprintf(stderr, "%s\t%s\t%s [%i] in %s(): %s\n", str_time, log_level_word, filename, num_line, func, msg);
        }

        // if 'log_file' has been set up then logs are printed
        // logs are ALWAYS printed in log_file independently of 'verbose' parameter
        if(log_filename != NULL) {
            FILE *log_file = fopen(log_filename, "a");
            fprintf(log_file, "%s\t%s\t%s [%i] in %s(): %s\n", str_time, log_level_word, filename, num_line, func, msg);
            fclose(log_file);
        }

    }
}

void print_log_message_with_format(int level, char *log_level_word, char *filename, int num_line, const char *func, char *msg_format, ...) {
    if((level >= log_level) && (verbose || log_filename != NULL)) {
        time_t rawtime;
        time(&rawtime);
        char* str_time = ctime(&rawtime);
        chomp(str_time);

        va_list args;
        va_start(args, msg_format);
        
        // if 'verbose' logs are printed in stdout
        if(verbose) {
            fprintf(stderr, "%s\t%s\t%s [%i] in %s(): ", str_time, log_level_word, filename, num_line, func);
            if (args == NULL) printf("OMFGexplosion! -- v\n");
            vfprintf(stderr, msg_format, args);
        }

        // if 'log_file' has been set up then logs are printed
        // logs are ALWAYS printed in log_file independently of 'verbose'       
        if(log_filename != NULL) {
            va_list argsf;
            memcpy(argsf, args, sizeof(args));
            va_start(argsf, msg_format);
            
            FILE *log_file = fopen(log_filename, "a");
            fprintf(log_file, "%s\t%s\t%s [%i] in %s(): ", str_time, log_level_word, filename, num_line, func);
            vfprintf(log_file, msg_format, argsf);
            fclose(log_file);
            
            va_end(argsf);
        }
        
        va_end(args);
    }
}
