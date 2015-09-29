#ifndef UTILS_H
#define UTILS_H

#define MAX_CN                 6

#define SIGMA_RATIO_CN1 0.707107 // sigma_haploid = this * sigma_diploid
                                 // poisson distribution variance = mean
                                 // so var(hap) = 1/2 var(dip)
#define SIGMA_RATIO_CN3 1.224745
#define SIGMA_RATIO_CN4 1.414214
#define SIGMA_RATIO_CN5 1.581139
#define SIGMA_RATIO_CN6 1.732051

#define HOM_DEL_THRESHOLD  -0.98 // if the exponential distribution model for hom del coverage
                                 // is fit by the EM algorithm to have a very small mean,
                                 // it gets replaced by a uniform distribution from -1 to this

void missing_value_error(char *arg);
void invalid_value_error(char *arg);
int double_comp(const void *a, const void *b);
double median(double *arr, int len);
FILE* open_file(char *path);
int count_lines_in_file(FILE *file);
void read_sample_name(char *dest, char *src);

#endif
