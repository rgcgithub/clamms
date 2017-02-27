#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#define CHR_X 23
#define CHR_Y 24
#define CHR_M 25

#define X_ASCII 88
#define Y_ASCII 89
#define M_ASCII 77
#define T_ASCII 84
#define DIGIT_ASCII_OFFSET 48

// read chromosome id into unsigned char
// and skip the next character after
// (':' in GATK coverage input, '\t' in samtools depth and windows input)
int read_chr(FILE *input, unsigned char *chr) {
    int tmp = getc(input);
    if (tmp == EOF) return 0;

    if (tmp == X_ASCII) { *chr = CHR_X; getc(input); return 1; }
    if (tmp == Y_ASCII) { *chr = CHR_Y; getc(input); return 1; }
    if (tmp == M_ASCII) {
        *chr = CHR_M;
        tmp = getc(input);
        if (tmp == T_ASCII) getc(input);
        return 1;
    }
    if (tmp < DIGIT_ASCII_OFFSET || tmp > (DIGIT_ASCII_OFFSET+9)) {
        *chr = 0;
        tmp = getc(input);
        return 2; // Illegal chromosome character. Likely the header.
    }

    tmp -= DIGIT_ASCII_OFFSET;
    if (tmp < 3) {
        int tmp2 = getc(input) - DIGIT_ASCII_OFFSET;
        if (tmp2 >= 0 && tmp2 <= 9) {
            *chr = 10*tmp + tmp2;
            getc(input);
        } else {
            *chr = tmp;
        }
    } else {
        *chr = tmp;
        getc(input);
    }

    return 1;
}

void skip_to_end_of_line(FILE *input) {
    int tmp;
    while ((tmp = getc(input)) != '\n') {
        if (tmp == EOF) {
            fputs("ERROR: malformed input\n", stderr);
            exit(1);
        }
    }
}

void print_chr(unsigned char chr,
               unsigned char *window_chr,
               int *window_start,
               int *window_end,
               double *window_cov,
               int *chr_start_idx) {
    int i;
    char chr_str[64];
         if (chr == CHR_X) { chr_str[0] = 'X'; chr_str[1] = '\0'; }
    else if (chr == CHR_Y) { chr_str[0] = 'Y'; chr_str[1] = '\0'; }
    else if (chr == CHR_M) {
        chr_str[0] = 'M'; chr_str[1] = 'T'; chr_str[2] = '\0';
    } else sprintf(chr_str, "%hhu", chr);

    i = chr_start_idx[chr];
    while (window_chr[i] == chr) {
        printf("%s\t%d\t%d\t%.6g\n",
               chr_str, window_start[i], window_end[i], window_cov[i]);
        i++;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s sample.gatk_readDepth_1x_q30.out windows.bed\n\n", argv[0]);
        fputs("Computes depth of coverage for intervals specified in windows.bed.\n", stderr);
        fputs("Also changes chr sort order from GATK's 1,2,3.. to BED's 1,10,11,..\n\n", stderr);
        return 1;
    }

    FILE *base_cov = fopen(argv[1], "r");
    if (base_cov == NULL) {
        fprintf(stderr, "Cannot read coverage file: %s\n", argv[1]);
        return 1;
    }


    FILE *windows = fopen(argv[2], "r");
    if (windows == NULL) {
        fprintf(stderr, "Cannot read windows file: %s\n", argv[2]);
        return 1;
    }

    char *line = NULL;
    char *pos;
    size_t line_len;
    ssize_t bytes_read;

    int i, j;
    int tmp_char;
    int n_windows = 0;

    while ((tmp_char = getc(windows)) != EOF) {
        if (tmp_char == '\n') n_windows++;
    } rewind(windows);

    int chr_start_idx[26];
    for (i = 0; i < 26; i++)
        chr_start_idx[i] = -1;

    unsigned char last_chr = 0;
    unsigned char *window_chr = (unsigned char *) malloc(n_windows * sizeof(unsigned char));
    int *window_start  = (int *) malloc(n_windows * sizeof(int));
    int *window_end    = (int *) malloc(n_windows * sizeof(int));
    double *window_cov = (double *) malloc(n_windows * sizeof(double));
    memset(window_cov, 0.0, n_windows*sizeof(double));

    i = 0;
    while (read_chr(windows, window_chr+i)) {
        fscanf(windows, "%d\t%d", window_start+i, window_end+i);
        skip_to_end_of_line(windows);
        if (window_chr[i] != last_chr) {
            last_chr = window_chr[i];
            chr_start_idx[window_chr[i]] = i;
        }
        i++;
    }

    unsigned char chr;
    last_chr = 0;
    int locus, read_depth;
    int cur_window, tot_cov, n_bases;
    n_bases = 0;
    // Check for header, advance line if it exists
    if (read_chr(base_cov, &chr) < 2)
        rewind(base_cov);
    else
        skip_to_end_of_line(base_cov);
    while (read_chr(base_cov, &chr)) {
        fscanf(base_cov, "%d\t%d", &locus, &read_depth);
        skip_to_end_of_line(base_cov);
        
        if (chr != last_chr) {
            if (n_bases > 0 && window_cov[cur_window] == 0.0)
                window_cov[cur_window] = (double) tot_cov / (double) n_bases;
            last_chr = chr;
            cur_window = chr_start_idx[chr];
            tot_cov = 0; n_bases = 0;
        } else if (chr != window_chr[cur_window]) {
            continue;
        }

        if (locus > window_end[cur_window]) {
            if (n_bases == 0)
                window_cov[cur_window] = 0.;
            else
                window_cov[cur_window] = (double) tot_cov / (double) n_bases;
            cur_window++;
            while (locus > window_end[cur_window] && chr == window_chr[cur_window])
                cur_window++;
            tot_cov = 0; n_bases = 0;
        }

        if (locus < window_start[cur_window]+1) continue;

        tot_cov += read_depth;
        n_bases++;
    }
    if (n_bases > 0 && window_cov[cur_window] == 0.0)
        window_cov[cur_window] = (double) tot_cov / (double) n_bases;

    print_chr(1, window_chr, window_start, window_end, window_cov, chr_start_idx);
    for (i = 10; i <= 19; i++)
        print_chr(i, window_chr, window_start, window_end, window_cov, chr_start_idx);
    print_chr(2, window_chr, window_start, window_end, window_cov, chr_start_idx);
    for (i = 20; i <= 22; i++)
        print_chr(i, window_chr, window_start, window_end, window_cov, chr_start_idx);
    for (i = 3; i <= 9; i++)
        print_chr(i, window_chr, window_start, window_end, window_cov, chr_start_idx);
    print_chr(CHR_M, window_chr, window_start, window_end, window_cov, chr_start_idx);
    print_chr(CHR_X, window_chr, window_start, window_end, window_cov, chr_start_idx);
    print_chr(CHR_Y, window_chr, window_start, window_end, window_cov, chr_start_idx);

    free(window_chr);
    free(window_start);
    free(window_end);
    free(window_cov);

    return 0;
}
