#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "utils.h"

typedef struct {
    double min_gc;        // %
    double max_gc;        // %
    int gc_window_size;   // # bases in the window used to calculate GC %
    int n_gc_bins;
    double min_mappability;
} Options;

typedef struct { 
    int n_bins;
    int min_gc;    // # bases
    int max_gc;    // # bases
    int gc_range;  // max_gc - min_gc        
    int bin_size;  // # distinct values per bin
    int bin_size_odd; // boolean
    int left_edge_bin_size;  // bins at the left and right edges
    int right_edge_bin_size; // may be larger than other bins
                             // if (gc_range+1) % n_bins != 0
    double left_edge_center;
    double right_edge_center;
} GCBinParams;

Options parse_args(int argc, char *argv[], int arg_start) {
    Options options;
    options.min_gc = 0.3;
    options.max_gc = 0.7;
    options.gc_window_size = 200;
    options.min_mappability = 0.75;

    int i;
    int n_gc_bins_set = 0;
    for (i = arg_start; i < argc; i += 2) {
        if (strcmp(argv[i], "--min_gc") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.min_gc = strtod(argv[i+1], NULL);
            if (options.min_gc <= 0.0 || options.min_gc > 1.0)
                invalid_value_error(argv[i]);
        } else if (strcmp(argv[i], "--max_gc") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.max_gc = strtod(argv[i+1], NULL);
            if (options.max_gc <= 0.0 || options.max_gc > 1.0)
                 invalid_value_error(argv[i]);
        } else if (strcmp(argv[i], "--n_gc_bins") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.n_gc_bins = (int) strtol(argv[i+1], NULL, 10);
            if (options.n_gc_bins < 10) invalid_value_error(argv[i]);
            n_gc_bins_set = 1;
        } else if (strcmp(argv[i], "--min_mappability") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.min_mappability = strtod(argv[i+1], NULL);
            if (options.min_mappability <= 0.0 || options.min_mappability > 1.0)
                invalid_value_error(argv[i]);
        } else {
            fprintf(stderr, "Unrecognized argument: %s\n", argv[i]);
            fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
            exit(1);
        }
    }

    if (!n_gc_bins_set)
        options.n_gc_bins = (int) floor((1.0 + options.gc_window_size *
                                        (options.max_gc - options.min_gc)) / 3.0);

    return options;
}

GCBinParams calc_gc_bin_params(Options *options) {
    GCBinParams gc_bin_params;
    gc_bin_params.n_bins = options->n_gc_bins;
    gc_bin_params.min_gc = (int)  ceil(options->min_gc * options->gc_window_size);
    gc_bin_params.max_gc = (int) floor(options->max_gc * options->gc_window_size);
    gc_bin_params.gc_range = gc_bin_params.max_gc - gc_bin_params.min_gc;
    gc_bin_params.bin_size = (int) floor((gc_bin_params.gc_range + 1) / gc_bin_params.n_bins);
    gc_bin_params.bin_size_odd = (gc_bin_params.bin_size % 2 == 1 ? 1 : 0);
    gc_bin_params.left_edge_bin_size = gc_bin_params.bin_size;
    gc_bin_params.right_edge_bin_size = gc_bin_params.bin_size;
    int remainder = (gc_bin_params.gc_range + 1) % options->n_gc_bins;
    if (remainder % 2 == 0) {
        gc_bin_params.left_edge_bin_size  += remainder / 2;
        gc_bin_params.right_edge_bin_size += remainder / 2;
    } else {
        gc_bin_params.left_edge_bin_size  += (remainder + 1) / 2;
        gc_bin_params.right_edge_bin_size += (remainder - 1) / 2;
    }

    gc_bin_params.left_edge_center = (double) (gc_bin_params.left_edge_bin_size - 1) / 2.0;
    gc_bin_params.right_edge_center =    (double) gc_bin_params.gc_range
                                      - ((double) gc_bin_params.right_edge_bin_size - 1) / 2.0;

    return gc_bin_params;
}

int get_gc_bin(int gc, GCBinParams *params) {
    int x = gc - params->min_gc;
    if (x < params->left_edge_bin_size) {
        return 0;
    } else if (x > params->gc_range - params->right_edge_bin_size) {
        return params->n_bins - 1;
    } else {
        x -= params->left_edge_bin_size;
        return 1 + x / params->bin_size;
    };
}

// estimate median(coverage | gc)
//     from median(coverage | gc bin)
double get_normalizing_factor(int gc, double *gc_meds, GCBinParams *params) {
    int x = gc - params->min_gc;
    if (x < params->left_edge_bin_size + params->bin_size / 2) {
        double rise  = gc_meds[1] - gc_meds[0];
        double  run  =   (double) params->left_edge_bin_size / 2.0
                       + (double) params->bin_size / 2.0;
        double slope = rise / run;
        return gc_meds[0] + slope * ((double) x - params->left_edge_center);
    } else if (x > params->gc_range - (params->right_edge_bin_size + params->bin_size / 2)) {
        double rise  = gc_meds[params->n_bins-2] - gc_meds[params->n_bins-1];
        double  run  =   (double) params->right_edge_bin_size / 2.0
                       + (double) params->bin_size / 2.0;
        double slope = rise / run;
        return gc_meds[params->n_bins-1] + slope * (params->right_edge_center - (double) x);
    } else {
        x -= params->left_edge_bin_size;
        int bin = 1 + x / params->bin_size;
        int remainder = x % params->bin_size;
        if (params->bin_size_odd && remainder == params->bin_size / 2) {
            return gc_meds[bin];
        } else {
            double t = ((double) remainder + (double) (remainder+1)) / (2.0 * params->bin_size);
            if (t < 0.5) {
                bin--;
                t += 0.5;
            } else {
                t -= 0.5;
            }
            return (1.0 - t) * gc_meds[bin] + t * gc_meds[bin+1];
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s coverage.bed windows.bed [OPTIONS] >normalized.coverage.bed\n\n", argv[0]);
        fputs("Normalizes a sample's coverage track relative to it's overall median depth\n", stderr);
        fputs("and also corrects coverage biases due to sequence GC content.\n\n", stderr);
        fputs("  --min_gc             Windows with GC fraction less than this are filtered.\n", stderr);
        fputs("                       Default = 0.3\n", stderr);
        fputs("  --max_gc             Windows with GC fraction greater than this are filtered.\n", stderr);
        fputs("                       Default = 0.7\n", stderr);
        fputs("  --n_gc_bins          # of bins to use when estimating the gc bias curve.\n", stderr);
        fputs("                       Default = floor((1+200*(max_gc-min_gc))/3).\n", stderr);
        fputs("  --min_mappability    Windows with mean mappability less than this are filtered.\n", stderr);
        fputs("                       Default = 0.75\n", stderr);
        fputs("                       Mappability scores are taken from windows.bed (column 8).\n\n", stderr);
        return 1;
    }

    FILE *coverage = open_file(argv[1]);
    FILE *windows  = open_file(argv[2]);

    Options options = parse_args(argc, argv, 3);

    char *line = NULL;
    char *pos, *last_pos;
    size_t line_len;
    ssize_t bytes_read;

    int i, j;
    int max_mm_cn;
    int n_windows = count_lines_in_file(windows);

    // read window info
    char   **window_coords      = (char  **) malloc(n_windows * sizeof(char *));
    int     *window_gc          = (int    *) malloc(n_windows * sizeof(int));
    double  *window_mappability = (double *) malloc(n_windows * sizeof(double));  
    int     *window_blacklisted = (int    *) malloc(n_windows * sizeof(int));
    for (i = 0; i < n_windows; i++)
        window_blacklisted[i] = 0;

    i = 0;
    while ((bytes_read = getline(&line, &line_len, windows)) != -1) {
        pos = strchr(line, '\t');
        pos = strchr(pos+1, '\t');
        pos = strchr(pos+1, '\t');
        int coord_str_len = (int) (pos - line);
        window_coords[i] = (char *) malloc((coord_str_len+1) * sizeof(char));
        strncpy(window_coords[i], line, coord_str_len);
        window_coords[i][coord_str_len] = '\0';

        pos = strchr(pos+1, '\t');
        window_gc[i] = (int) strtol(pos, &pos, 10);
        pos = strchr(pos+1, '\t');
        window_mappability[i] = strtod(pos, &pos);
        if (*pos == '\t') {
            sscanf(++pos, "%d", &max_mm_cn);
            if (max_mm_cn < 0)
                window_blacklisted[i] = 1;
        }
        i++;
    }
    
    // initialize gc bins
    GCBinParams gc_bin_params = calc_gc_bin_params(&options);
    int max_bin_size = n_windows / 5; // don't want to implement resizeable array (or C++)
                                      // so I just allocate a big chunk of memory and assume it'll be enough
                                      // it's only 1 sample being processed, so memory usage shouldn't be an issue
    double **gc_bins = (double **) malloc(gc_bin_params.n_bins * sizeof(double *));
    int *gc_bin_n_elems  = (int *) malloc(gc_bin_params.n_bins * sizeof(int));
    for (i = 0; i < gc_bin_params.n_bins; i++) {
        gc_bins[i] = (double *) malloc(max_bin_size * sizeof(double));
        gc_bin_n_elems[i] = 0;
    }

    // read coverage values and put them into gc bins
    i = 0; // window # (autosome only)
    j = 0; // line #
    double *cov = (double *) malloc(n_windows * sizeof(double));
    while ((bytes_read = getline(&line, &line_len, coverage)) != -1) {
         j++;

         char chr = *line;
         pos = strchr(line, '\t');
         pos = strchr(pos+1, '\t');
         pos = strchr(pos+1, '\t');
         int coord_str_len = (int) (pos - line);
         if (strncmp(window_coords[i], line, coord_str_len) != 0) {
             fprintf(stderr, "Coordinates              [ %*.*s ] at line %d of coverage file\n"
                             "do not match coordinates [ %s ] at line %d of windows file\n",
                     coord_str_len, coord_str_len, line, j, window_coords[i], j);
             break;
         }
         
        cov[i] = strtod(pos, &pos);
        if (window_gc[i] < gc_bin_params.min_gc ||
            window_gc[i] > gc_bin_params.max_gc ||
            window_mappability[i] < options.min_mappability) {
            window_blacklisted[i] = 1;
        } else if (!(window_blacklisted[i] || chr == 'X' || chr == 'Y')) {
            int bin = get_gc_bin(window_gc[i], &gc_bin_params);
            gc_bins[bin][gc_bin_n_elems[bin]++] = cov[i];
        }
        
        i++;
    }

    // compute median(coverage | sample, gc_bin)
    double *gc_meds = (double *) malloc(gc_bin_params.n_bins * sizeof(double));
    for (i = 0; i < options.n_gc_bins; i++) {
        qsort(gc_bins[i], gc_bin_n_elems[i], sizeof(double), double_comp);
        gc_meds[i] = median(gc_bins[i], gc_bin_n_elems[i]);
    }

// for debugging the gc-bias curve estimation
/*
    for (i = 61; i <= 141; i++) {
        fprintf(stderr, "%d\t%.2f\t%.2f\n", i, gc_meds[get_gc_bin(i, &gc_bin_params)], get_normalizing_factor(i, gc_meds, &gc_bin_params));
    }
*/
    // output normalized coverage values
    for (i = 0; i < n_windows; i++) {
        printf("%s\t%.6g\n",
               window_coords[i],
               cov[i] / get_normalizing_factor(window_gc[i], gc_meds, &gc_bin_params));
    }

    // cleanup

    if (line) free (line);
    for (i = 0; i < n_windows; i++)
        free(window_coords[i]);
    for (i = 0; i < options.n_gc_bins; i++)
        free(gc_bins[i]);
    free(window_coords);
    free(window_gc);
    free(window_mappability);
    free(gc_bins);
    free(gc_bin_n_elems);
    free(gc_meds);
    free(cov);

    fclose(coverage);
    fclose(windows);
}

