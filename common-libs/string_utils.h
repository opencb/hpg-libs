
#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

int equals(const char *str1, const char *str2);

int equals_ignore_case(const char *str1, const char *str2);

int is_numeric(const char *str);


int contains(const char *str, const char *search);

int starts_with(const char *str, const char *search);

int ends_with(const char *str, const char *search);


char* to_lower_case(char *str);

char* to_upper_case(char *str);


char chop(char *str);

char chop_at(char *str, int position);

char chomp(char *str);

char chomp_at(char *str, int position);


char* remove_char(char *str, char c);

char* remove_char_at(char *str, int position);

char* remove_str(char *str, const char *search_str);

char* remove_start(char *str, int lentgh);

char* remove_end(char *str, int length);


char* str_replace(char *st, const char *orig, const char *repl, int max_line_length);

int array_concat(char **dest, int orig1_length, const char **orig1, int orig2_length, const char **orig2);


char* trim(char *str);

char* ltrim2(char *str);

char* rtrim2(char *str);

char* strip(char *str);

char* lstrip(char *str);

char* rstrip(char *str);


char** split(const char *str);

char** splitn(const char *str, int limit);

/**
 * Case-insensitive string comparison. Inspired in non-standard function:
 * http://www.opensource.apple.com/source/gpatch/gpatch-2/patch/strcasecmp.c
 * 
 * @param s1
 * 	First string to compare
 * @param s2
 * 	Second string to compare
 * 
 * @return
 * 	Less, equals or greater than zero if s1 is lexicographically less 
 * 	than, equal to or greater than s2.
 */
int strcasecmp(const char *s1, const char *s2);
//int strcasecmp(const char *s1, const char *s2);
//char strcasecmp(const char *s1, const char *s2);

//-----------------------------------------------------------
// functions to encode and decode nucleotied sequences
//-----------------------------------------------------------

enum bases{DD=-1,AA=0,CC=1,GG=2,TT=3}; //Lexicographic order
static const char alph_rep[] ={'A','C','G','T'};


void initTable();
char* encodeBases(char *dest, char* src, unsigned int length);
char* decodeBases(char *dest, char* src, unsigned int length);

#endif	/*    STRING_UTILS_H	*/
