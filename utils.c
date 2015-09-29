#include "stdio.h"
#include "stdlib.h"
#include "string.h"

void missing_value_error(char *arg) {
    fprintf(stderr, "Missing value for argument: %s\n", arg);
    exit(1);
}

void invalid_value_error(char *arg) {
    fprintf(stderr, "Invalid value for argument: %s\n", arg);
    exit(1);
}

// used by qsort
int double_comp(const void *a, const void *b) {
    if (*((const double *) a) < *((const double *) b))
        return -1;
    return *((const double *) a) > *((const double *) b);
}

// input must be sorted
double median(double *arr, int len) {
    if (len % 2 == 1)
        return arr[(len-1)/2];
    else
        return (arr[len/2] + arr[len/2 - 1]) / 2.0;
}

FILE *open_file(char *path) {
    FILE *file = fopen(path, "r");
    if (file == NULL) {
        fprintf(stderr, "Cannot read file: %s\n", path);
        exit(1);
    }
    return file;
}

int count_lines_in_file(FILE *file) {
    int n_lines = 0;
    int tmp_char;
    while ((tmp_char = getc(file)) != EOF) {
        if (tmp_char == '\n') n_lines++;
    } rewind(file);
    return n_lines;
}

void read_sample_name(char *dest, char *src) {
    char *slash_pos;
    while ((slash_pos = strchr(src, '/')) != NULL)
        src = slash_pos + 1;
    char *first_dot = strchr(src, '.');
    if (first_dot == NULL) {
        strcpy(dest, src);
    } else {
        int len = first_dot - src;
        strncpy(dest, src, len);
        dest[len] = '\0';
    }
}
