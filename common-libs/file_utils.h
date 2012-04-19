#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "string_utils.h"


#define MAX_LENGTH_CONFIG_LINE 		256
#define MAX_FILENAME_LENGTH			256
#define MAX_FULL_PATH_LENGTH		2048


typedef struct launch_option {
    char option_name[MAX_LENGTH_CONFIG_LINE+2];
    char option_value[MAX_LENGTH_CONFIG_LINE];
} launch_option_t;


void *mmap_file(size_t *len, const char *filename);


char* fgets_no_ln(char *s, int n, FILE *f);

int exists(const char *path);

int is_file(const char *path);

int is_directory(const char *path);

unsigned long count_lines(const char *filename);


int copy(const char *dest, const char *ori);

int move(const char *dest, const char *ori);


int touch(const char *path);

int create_directory(const char *path);

int delete_directory(const char *path);

int delete_files_by_extension(const char *dir_path, const char *extension);


unsigned long size(const char *path);

int create_directory(const char *path);


char** parse_conf_file(char *filename);

int parse_conf_file2(char **argvs, char *filename);

char* get_filename_from_path(char* path, char* filename_p);

#endif	/*  FILE_UTILS_H   */
