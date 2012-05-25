#include "string_utils.h"

int equals(const char *str1, const char *str2) {
    return strcmp(str1, str2) == 0;
}

int equals_ignore_case(const char *str1, const char *str2) {
    int str1_len = strlen(str1);
    int str2_len = strlen(str2);
    if (str1_len == str2_len) {
        int i = 0;
        for (i = 0; i < str1_len; i++) {
            if (str1[i] != str2[i]) {
                if (str1[i] >= 65 && str1[i] <= 90) {
                    if (str1[i] + 32 != str2[i]) {
                        return 0;
                    }
                } else {
                    if (str1[i] >= 97 && str1[i] <= 122) {
                        if (str1[i] - 32 != str2[i]) {
                            return 0;
                        }
                    } else {
                        // they are different and they are no letters
                        return 0;
                    }
                }
            }
        }
        // they have same size and chars ignoring case
        return 1;
    } else {
        // they have different sizes
        return 0;
    }
}

int is_numeric(const char *str) {
    int str_len = strlen(str);
    
    for (int i = 0; i < str_len; i++) {
        if (!isdigit(str[i])) {
            return 0;
        }
    }
    
    return 1;
}


int starts_with(const char *str, const char *search) {
    int str_len = strlen(str);
    int search_len = strlen(search);

    if (search_len > str_len) {
        return 0;
    }

    return strncmp(str, search, search_len) == 0;
}

int ends_with(const char *str, const char *search) {
    int str_len = strlen(str);
    int search_len = strlen(search);

    if (search_len > str_len) {
        return 0;
    }

    return strncmp(str + (str_len - search_len), search, search_len) == 0;
}

char* to_lower_case(char *str) {
    int str_len = strlen(str);

    for (int i = 0; i < str_len; i++) {
        if (str[i] >= 65 && str[i] <= 90) {
            str[i] += 32;
        }
    }
    return str;
}

char* to_upper_case(char *str) {
    int str_len = strlen(str);
    
    for (int i = 0; i < str_len; i++) {
        if (str[i] >= 97 && str[i] <= 122) {
            str[i] -= 32;
        }
    }
    return str;
}

char chop(char *str) {
    return chop_at(str, strlen(str) - 1);
}

char chop_at(char *str, int position) {
    char last_char = str[position];
    str[position] = '\0';
    return last_char;
}

char chomp(char *str) {
    return chomp_at(str, strlen(str) - 1);
}

char chomp_at(char *str, int position) {
    if (str[position] == '\n' || str[position] == '\r') {
        char new_line = str[position];
        str[position] = '\0';
        return new_line;
    } else {
        return '\0';
    }
}

char* remove_char(char *str, char c) {
    int str_len = strlen(str);
    int count = 0;
    
    for (int i = 0; i < str_len; i++) {
        if (str[i] == c) {
            count++;
        }
    }
    char *str_aux = (char*)malloc((str_len - count + 1) * sizeof(char));
    
    for (int i = 0, j = 0; i < str_len + 1; i++) {
        if (str[i] != c) {
            str_aux[j] = str[i];
            j++;
        }
    }
    
    strcpy(str, str_aux);
    free(str_aux);
    return str;
}

char* remove_char_at(char *str, int position) {
    int str_len = strlen(str);    
    char *str_aux = (char*)malloc((str_len) * sizeof(char));
    
    for (int i = 0, j = 0; i < str_len + 1; i++) {
        if (i != position) {
            str_aux[j] = str[i];
            j++;
        }
    }
    
    strcpy(str, str_aux);
    free(str_aux);
    return str;
}

char* remove_str(char *str, const char *search_str) {
    int str_length = strlen(str);
    int search_length = strlen(search_str);
    int num_repeats = 0;

    if (str_length > search_length) {
        int search_indices[10];
		int j;

        for (int i = 0; i < str_length - search_length + 1; i++) {
            j = 0;
            if (str[i] == search_str[j]) {
                for (j = 1; j < search_length; j++) {
                    if (str[i+j] != search_str[j]) {
                        break;
                    }
                }
                if (j == search_length) {
                    search_indices[num_repeats] = i;
                    num_repeats++;
                }
            }
        }
    }
    
    return str;
}

char* remove_start(char *str, int num_chars) {
    if (!str) {
        return NULL;
    }
    
    int i = 0;
    while (str[i] != '\0') {
        str[i] = str[num_chars];
        i++;
        num_chars++;
    }
    
    str[i] = '\0';
    return str;
}

char* remove_end(char *str, int num_chars) {
    int length = strlen(str);
    
    if (length > num_chars) {
        str[length - num_chars] = '\0';
    } else {
        str[0] = '\0';
    }
    
    return str;
}

char* str_replace(char *str, const char *orig, const char *repl, int max_line_length) {
    char buffer[max_line_length];
    char *ch;

    if (!(ch = strstr(str, orig))) {
        return str;
    }

    strncpy(buffer, str, ch - str);
    buffer[ch - str] = 0;
    sprintf(buffer + (ch - str), "%s%s", repl, ch + strlen(orig));

    return str_replace(buffer, (const char*)orig, (const char*)repl, max_line_length);
}

int array_concat(char **dest, int orig1_length, const char **orig1, int orig2_length, const char **orig2) {
    int i, j;
    for (i = 0; i < orig1_length; i++) {
        dest[i] = (char*)malloc(strlen(orig1[i]) * sizeof(char) + 1);

        strcpy(dest[i], orig1[i]);
    }

    for (i = orig1_length, j = 0; j < orig2_length; i++, j++) {
        dest[i] = (char *)malloc(strlen(orig2[j]) * sizeof(char) + 1);
        strcpy(dest[i], orig2[j]);
    }

    return i;
}

char* trim(char* str) {
    return ltrim2(rtrim2(str));
}

char* ltrim2(char *str) {
    if (!str) {
        return NULL;
    }
    
    int i = 0;
    int j = 0;
    
    while (isspace(str[j])) {
        j++;
    }
    
    while (str[j] != '\0') {
        str[i] = str[j];
        i++;
        j++;
    }
    
    str[i] = '\0';
    return str;
}

char* rtrim2(char *str) {
    int i = strlen(str) - 1;
    
    while (i >= 0 && isspace(str[i])) {
        --i;
    }
    
    str[i + 1] = '\0';
    return str;
}

char* strip(char *str) {
    return lstrip(rstrip(str));
}

char* lstrip(char *str) {
    if (!str) {
        return NULL;
    }
    
    int i = 0;
    int j = 0;
    
    while (str[j] == ' ') {
        j++;
    }
    
    while (str[j] != '\0') {
        str[i] = str[j];
        i++;
        j++;
    }
    
    str[i] = '\0';
    return str;
}

char* rstrip(char* str) {
    int i = strlen(str) - 1;
    
    while (i >= 0 && str[i] == ' ') {
        --i;
    }
    
    str[i+1] = '\0';
    return str;
}

char* ltrim(char* string, int num_chars) {
    int length = strlen(string) - num_chars;
    char* cut_string = (char*) malloc(length * sizeof(char));

    strncpy(cut_string, string + num_chars, length);
    string = cut_string;

    return string;
}

char* rtrim(char* string, int num_chars) {
    int index;
    index = strlen(string) - num_chars;

    string[index] = '\0';

    return string;
}


char** split(char* str, const char *delimiters, int *num_substrings) {
    return splitn(str, delimiters, INT_MAX, num_substrings);
}

char** splitn(char* str, const char *delimiters, int limit, int *num_substrings) {
    int i = 0;
    int max_substrings = 16;
    char **split_text = (char**) malloc (max_substrings * sizeof(char*));
    char **tmp_realloc;
    
    char *token, *token_dest;
    char *save_strtok;
    
    while ((token = strtok_r(str, delimiters, &save_strtok)) && i < limit) {
        token_dest = (char*) calloc (strlen(token)+1, sizeof(char));
        strcat(token_dest, token);
        split_text[i] = token_dest;
        
        str = NULL;
        i++;
        if (i == max_substrings) {
            printf("realloc'ing\n");
            tmp_realloc = realloc(split_text, (max_substrings + 16) * sizeof(char*));
            if (tmp_realloc) {
                max_substrings += 16;
                split_text = tmp_realloc;
            } else {
                // Impossibility of reallocation avoids the function to end successfully
                return NULL;
            }
            printf("realloc'd\n");
        }
    }
    
    // Reallocate memory for a little optimization
    tmp_realloc = realloc(split_text, (i + 1) * sizeof(char*));
    if (tmp_realloc) {
        split_text = tmp_realloc;
    }
    
    *num_substrings = i;
    return split_text;
}


int strcasecmp(const char *s1, const char *s2)
{
    const char *p1 = (const char *) s1;
    const char *p2 = (const char *) s2;
    int c1, c2;

    do
    {
        c1 = isupper(*p1) ? tolower(*p1) : *p1;
        c2 = isupper(*p2) ? tolower(*p2) : *p2;

        if (c1 == '\0' || c2 == '\0') {
            break;
        }

        ++p1;
        ++p2;
    } while (c1 == c2);

    return c1 - c2;
}

int table[128] = {};

void initTable() {
    table['a'] = AA;
    table['A'] = AA;
    table['c'] = CC;
    table['C'] = CC;
    table['t'] = TT;
    table['T'] = TT;
    table['g'] = GG;
    table['G'] = GG;
    table['n'] = AA;
    table['N'] = AA;
}

char* encodeBases(char* dest, char* src, unsigned int length) {    
    for (unsigned int i = 0; i < length; i++) {
        dest[i] = table[(int)src[i]];
    }
    
    dest[length] = '\0';
    return dest;
}

char* decodeBases(char* dest, char* src, unsigned int length) {
    for (unsigned int i = 0; i < length; i++) {
        dest[i] = alph_rep[(int)src[i]];
    }
    
    dest[length] = '\0';
    return dest;
}
