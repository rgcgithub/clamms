#ifndef HMM_H
#define HMM_H

#define N_STATES  3 // DEL, NORM, DUP
                    // hom/het distinctions for DEL and DUP are not part of the HMM
                    // but are rather estimated seperately for the CNV calls that the HMM makes
#define N_CHROM  25 // 1-22, X, Y, M
#define CHR_X    23
#define CHR_Y    24
#define CHR_M    25

// array indexes for HMM states
#define DEL  0
#define NORM 1
#define DUP  2

// possible values for NORM state
// NORM = diploid for autosome and female chrX
// NORM = haploid for male chrX/Y and both sexes chrM
// NORM = not present for female chrY
#define NOT_PRESENT 0
#define HAPLOID     1
#define DIPLOID     2

// directions to run Viterbi algorithm
#define FORWARD   1
#define BACKWARD -1

// variables for defining GC-content weights
#define GC_BUFFER 0.01

int get_next_window(int window, int n_windows,
                    unsigned char *window_chr, char *max_cn);
int get_prev_window(int window, int n_windows,
                    unsigned char *window_chr, char *max_cn);
double transition_prob(int from, int to, double p, double f);
double homozygous_del_log_likelihood(double cov, int dist_is_point, double lambda);
double gaussian_log_likelihood(double cov, int cn, double sigma_dip);
int expected_copy_number(char sex, unsigned char chr);
void read_model_data(FILE *input,
                     int n_windows,
                     unsigned char *window_chr,
                     int *window_start,
                     int *window_end,
                     char *max_cn,
                     unsigned char *hom_del_flag,
                     double *window_gc,
                     double *lambda,
                     double *mu_dip,
                     double *sigma_dip,
                     double *model_conf);
void read_coverage_data(FILE *input,
                        int n_windows,
                        unsigned char *window_chr,
                        int *window_start,
                        int *window_end,
                        double *cov,
                        double *mu_dip);
void calc_base_model_conf(int n_windows,
                          double gc_min,
                          double gc_max,
                          unsigned char *window_chr,
                          int *window_start,
                          int *window_end,
                          char *max_cn,
                          double *window_gc,
                          double *model_conf);
void calc_sample_specific_model_conf(int n_windows,
                                     char sex,
                                     unsigned char *window_chr,
                                     char *max_cn,
                                     double *cov,
                                     unsigned char *hom_del_flag,
                                     double *lambda,
                                     double *sigma_dip,
                                     double *model_conf);
void calc_cn_emission_logp(int n_windows,
                           char sex,
                           unsigned char *window_chr,
                           char *max_cn,
                           double *cov,
                           unsigned char *hom_del_flag,
                           double *lambda,
                           double *sigma_dip,
                           double **cn_emission_logp);
void calc_hmm_state_emission_logp(int n_windows,
                                  char sex,
                                  unsigned char *window_chr,
                                  char *max_cn,
                                  double **cn_emission_logp,
                                  double **hmm_state_emission_logp);
unsigned char *viterbi(int n_windows,
                       int direction,
                       unsigned char *window_chr,
                       int *window_start,
                       int *window_end,
                       char *max_cn,
                       double *model_conf,
                       double **hmm_state_emission_logp,
                       double cnv_rate,
                       double mean_cnv_length);
void mask_sequence(int n_windows, char *max_cn,
                   unsigned char *seq1, unsigned char *seq2);
void forward_backward(int n_windows,
                      unsigned char *window_chr,
                      int *window_start,
                      int *window_end,
                      char *max_cn,
                      double *model_conf,
                      double **hmm_state_emission_logp,
                      double cnv_rate,
                      double mean_cnv_length,
                      double **forward_scaled_prob,
                      double **backward_scaled_prob);

#endif
