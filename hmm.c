#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "hmm.h"
#include "utils.h"
#include "ltqnorm.c"

int get_next_window(int window, int n_windows,
                    unsigned char *window_chr,
                    char *max_cn) {
    int next = window + 1;
    while (max_cn[next] < 0 && next < n_windows-1)
        next++;
    if (max_cn[next] < 0 || window_chr[next] != window_chr[window])
        return -1;
    else
        return next;
}

int get_prev_window(int window, int n_windows,
                    unsigned char *window_chr,
                    char *max_cn) {
    int prev = window - 1;
    while (max_cn[prev] < 0 && prev > 0)
        prev--;
    if (max_cn[prev] < 0 || window_chr[prev] != window_chr[window])
        return -1;
    else
        return prev;
}

double transition_prob(int from, int to, double p, double f) {
    if (from == NORM) {
        if (to == NORM)         return 1.0 - 2.0*p;
        else                    return p;
    } else if (from == DEL) {
        if (to == DEL)          return f + (1.0 - f) * p;
        else if (to == NORM)    return (1.0 - f) * (1.0 - 2.0*p);
        else                    return (1.0 - f) * p;
    } else if (from == DUP) {
        if (to == DUP)          return f + (1.0 - f) * p;
        else if (to == NORM)    return (1.0 - f) * (1.0 - 2.0*p);
        else                    return (1.0 - f) * p;
    } else {
        fprintf(stderr, "Invalid call to transition_prob: from %d to %d\n", from, to);
        exit(1);
    }
}

// cov should be already centered+scaled
double homozygous_del_log_likelihood(double cov, int dist_is_point, double lambda) {
    if (dist_is_point) {
        if (cov < -0.98)    return  3.912023; // log(50)
        else if (cov < 0.0) return -13.81551; // log(1e-6)
        else                return -100.0;
    } else {
        double log_lambda = log(lambda);
        return log_lambda - lambda * (cov + 1.0);
    }
}

// cov should be already centered+scaled
double gaussian_log_likelihood(double cov, int cn, double sigma_dip) {
    double var = sigma_dip * sigma_dip;
    double const_term = -0.5 * log(M_2_PI * var);
    double scale_term = 1.0 / (2.0 * var);
    double dev = cov + 1.0 - 0.5 * cn;
    return const_term - scale_term * dev * dev;
}

int expected_copy_number(char sex, unsigned char chr) {
    if (sex == 'F' && chr == CHR_Y)
        return NOT_PRESENT;
    else if ((sex == 'M' && (chr == CHR_X || chr == CHR_Y)) || chr == CHR_M)
        return HAPLOID;
    else
        return DIPLOID;
}

double gc_confidence_term(double gc, double gc_min, double gc_max) {
    if (gc < gc_min || gc > gc_max) return 0.0;
    double constant;
    if (gc <= 0.5) {
        constant = 1.0/( 0.5 - gc_min + GC_BUFFER );
    }
    else if (gc > 0.5) {
        constant = 1.0/( gc_max - 0.5 + GC_BUFFER );
    }
    double x = constant * fabs(gc - 0.5);
    double x_2  = x    * x;
    double x_4  = x_2  * x_2;
    double x_8  = x_4  * x_4;
    double x_16 = x_8  * x_8;
    double x_18 = x_16 * x_2;
    double y    = 1.0  - x_18;
    double y_2  = y    * y;
    double y_4  = y_2  * y_2;
    double y_8  = y_4  * y_4;
    double y_16 = y_8  * y_8;
    double y_18 = y_16 * y_2;
    if (y_18 < 0.0) return 0.0;
    else return y_18;
}

double cov_confidence_term(int expected_cn,
                           char max_cn,
                           double coverage,
                           unsigned char hom_del_flag,
                           double lambda,
                           double sigma_dip) {
    if (hom_del_flag && coverage < HOM_DEL_THRESHOLD)
        return 1.0;

    int i;
    int ml_nonzero_cn = 1;
    double ml_cov_sigma  = fabs(coverage + 0.5) / (sigma_dip * SIGMA_RATIO_CN1);
    for (i = 2; i <= max_cn; i++) {
        double model_sigma = sigma_dip;
             if (i == 3) model_sigma *= SIGMA_RATIO_CN3;
        else if (i == 4) model_sigma *= SIGMA_RATIO_CN4;
        else if (i == 5) model_sigma *= SIGMA_RATIO_CN5;
        else if (i == 6) model_sigma *= SIGMA_RATIO_CN6;
        double cov_sigma = fabs(coverage - (-1.0 + 0.5*i)) / model_sigma;
        if (cov_sigma < ml_cov_sigma) {
            ml_cov_sigma  = cov_sigma;
            ml_nonzero_cn = i;
        }
    }

    if (ml_nonzero_cn == 1 && !hom_del_flag) {
        double expo_tail_prob = exp(-lambda * (coverage + 1.0));
        double gaus_tail_prob = 0.5 * erfc(ml_cov_sigma / sqrt(2.0));
        if (expo_tail_prob > gaus_tail_prob) {
            if (expo_tail_prob > 0.5)
                return 1.0;
            else
                ml_cov_sigma = -ltqnorm(expo_tail_prob); // gaussian inverse cdf
        }
    }

    double max_sigma = 3.598858; // -> y_18 @ 3-sigma = 0.5
    if (ml_nonzero_cn > expected_cn)
        max_sigma = 2.9990483;   // -> y_18 @ 2.5-sigma = 0.5

    if (ml_cov_sigma < 0.0)
        return 1.0; // shouldn't happen, but checking just in case there's floating point error or something
    else if (ml_cov_sigma > max_sigma)
        return 0.0;

    double x = ml_cov_sigma / max_sigma;
    double x_2  = x    * x;
    double x_4  = x_2  * x_2;
    double x_8  = x_4  * x_4;
    double x_16 = x_8  * x_8;
    double x_18 = x_16 * x_2;
    double y    = 1.0  - x_18;
    double y_2  = y    * y;
    double y_4  = y_2  * y_2;
    double y_8  = y_4  * y_4;
    double y_16 = y_8  * y_8;
    double y_18 = y_16 * y_2;
    if (y_18 < 0.0) return 0.0;
    else return y_18;
}

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
                     double *model_conf) {
    int i = 0;
    char *line = NULL;
    char *pos;
    size_t line_len;
    ssize_t bytes_read;

    while ((bytes_read = getline(&line, &line_len, input)) != -1) {
             if (*line == 'X') window_chr[i] = CHR_X;
        else if (*line == 'Y') window_chr[i] = CHR_Y;
        else if (*line == 'M') window_chr[i] = CHR_M;
        else sscanf(line, "%hhu", window_chr + i);

        pos = line;
        pos = strchr(pos+1, '\t');
        sscanf(pos+1, "%u", window_start + i);
        pos = strchr(pos+1, '\t');
        sscanf(pos+1, "%u", window_end + i);
        pos = strchr(pos+1, '\t');
        pos++;
        sscanf(pos, "%hhd", max_cn + i);
        //if (max_cn[i] == 2) max_cn[i] == 3; // we limit common del regions to CN 2
                                            // for purposes of fitting model parameters
                                            // but we allow the possibility of a rare CN 3
                                            // when making the final CNV calls
        window_gc[i] = strtod(pos+2, &pos);
        pos = strchr(pos+1, '\t');
        pos++;
        sscanf(pos, "%hhu", hom_del_flag + i);
        lambda[i] = strtod(pos+1, &pos);
        mu_dip[i] = strtod(pos+1, &pos);
        sigma_dip[i] = strtod(pos+1, &pos);

        i++;
    }
}

void read_coverage_data(FILE *input,
                        int n_windows,
                        unsigned char *window_chr,
                        int *window_start,
                        int *window_end,
                        double *cov,
                        double *mu_dip) {
    int i;
    char *line = NULL;
    char *pos;
    size_t line_len;
    ssize_t bytes_read;

    for (i = 0; i < n_windows; i++) {
        bytes_read = getline(&line, &line_len, input);
        if (bytes_read == -1) {
            fputs("ERROR: coverage file and models file must be for exactly the same set of genomic windows\n.", stderr);
            exit(1);
        }

        unsigned char chr;
        int start, end;

             if (*line == 'X') chr = CHR_X;
        else if (*line == 'Y') chr = CHR_Y;
        else if (*line == 'M') chr = CHR_M;
        else sscanf(line, "%hhu", &chr);

        pos = line;
        pos = strchr(pos+1, '\t');
        sscanf(pos+1, "%u", &start);
        pos = strchr(pos+1, '\t');
        sscanf(pos+1, "%u", &end);

        if (chr != window_chr[i] || start != window_start[i] || end != window_end[i]) {
            fputs("ERROR: coverage file and models file must be for exactly the same set of genomic windows\n.", stderr);
            exit(1);
        }

        pos = strchr(pos+1, '\t');
        cov[i] = (strtod(pos+1, NULL) - mu_dip[i]) / mu_dip[i];
    }
}

void calc_base_model_conf(int n_windows,
                          double gc_min,
                          double gc_max,
                          unsigned char *window_chr,
                          int *window_start,
                          int *window_end,
                          char *max_cn,
                          double *window_gc,
                          double *model_conf) {
    int i;
    for (i = 0; i < n_windows; i++) {
        if (max_cn[i] < 0) continue;
        int prev_window = get_prev_window(i, n_windows, window_chr, max_cn);
        int next_window = get_next_window(i, n_windows, window_chr, max_cn);
        int prev_window_end   = -1;
        int next_window_start = -1;
        if (prev_window != -1) prev_window_end   = window_end[prev_window];
        if (next_window != -1) next_window_start = window_start[next_window];
        
        model_conf[i] = gc_confidence_term(window_gc[i], gc_min, gc_max);
        if (window_start[i] == prev_window_end || window_end[i] == next_window_start)
            model_conf[i] = fmin(model_conf[i], 2.0/3.0);
    }
}

// coverage should already be centered+scaled around mu_dip
void calc_sample_specific_model_conf(int n_windows,
                                     char sex,
                                     unsigned char *window_chr,
                                     char *max_cn,
                                     double *cov,
                                     unsigned char *hom_del_flag,
                                     double *lambda,
                                     double *sigma_dip,
                                     double *model_conf) {
    int i;
    for (i = 0; i < n_windows; i++) {
        if (max_cn[i] < 0) continue;
        double cov_conf = cov_confidence_term(
            expected_copy_number(sex, window_chr[i]),
            max_cn[i], cov[i], hom_del_flag[i], lambda[i], sigma_dip[i]);
        model_conf[i] = fmin(model_conf[i], cov_conf);
    }
}

void calc_cn_emission_logp(int n_windows,
                           char sex,
                           unsigned char *window_chr,
                           char *max_cn,
                           double *cov,
                           unsigned char *hom_del_flag,
                           double *lambda,
                           double *sigma_dip,
                           double **cn_emission_logp) {
    int i, j;
    for (i = 0; i < n_windows; i++) {
        if (max_cn[i] < 0) continue;
        int norm_cn = expected_copy_number(sex, window_chr[i]);
        if (norm_cn == HAPLOID) max_cn[i] = 2;

        // first just compute the likelihoods

        cn_emission_logp[i][0] = homozygous_del_log_likelihood(cov[i], hom_del_flag[i], lambda[i]);
        for (j = 1; j <= MAX_CN; j++)
            cn_emission_logp[i][j] = gaussian_log_likelihood(cov[i], j, sigma_dip[i]);

        // now use bayes theorem to go from likelihoods to actual probabilities
        // we use a uniform prior, since our prior beliefs about CNV probabilities
        // are already encoded in the transition matrix of the HMM
        double evidence = 0.0;
        for (j = 0; j <= max_cn[i]; j++)
            evidence += exp(cn_emission_logp[i][j]);
        evidence = log(evidence);
        for (j = 0; j <= max_cn[i]; j++)
            cn_emission_logp[i][j] -= evidence;
    }
}

void calc_hmm_state_emission_logp(int n_windows,
                                  char sex,
                                  unsigned char *window_chr,
                                  char *max_cn,
                                  double **cn_emission_logp,
                                  double **hmm_state_emission_logp) {
    int i, j;
    for (i = 0; i < n_windows; i++) {
        int norm_cn = expected_copy_number(sex, window_chr[i]);
        if (norm_cn == NOT_PRESENT) {
            hmm_state_emission_logp[i][DEL]  = -100.0;
            hmm_state_emission_logp[i][NORM] = 0.0;
            hmm_state_emission_logp[i][DUP]  = -100.0;
        } else if (norm_cn == HAPLOID) {
            hmm_state_emission_logp[i][DEL]  = cn_emission_logp[i][0];
            hmm_state_emission_logp[i][NORM] = cn_emission_logp[i][1];
            hmm_state_emission_logp[i][DUP]  = cn_emission_logp[i][2];
        } else {
            hmm_state_emission_logp[i][DEL]  = log(exp(cn_emission_logp[i][0]) +
                                                   exp(cn_emission_logp[i][1]));
            hmm_state_emission_logp[i][NORM] = cn_emission_logp[i][2];
            double dup_tmp = 0.0;
            for (j = 3; j <= max_cn[i]; j++)
                dup_tmp += exp(cn_emission_logp[i][j]);
            hmm_state_emission_logp[i][DUP] = log(dup_tmp);
        }
    }
}

unsigned char *viterbi(int n_windows,
                       int direction,
                       unsigned char *window_chr,
                       int *window_start,
                       int *window_end,
                       char *max_cn,
                       double *model_conf,
                       double **hmm_state_emission_logp,
                       double cnv_rate,
                       double mean_cnv_length) {
    int i, j;
    unsigned char **ml_prev_state  = (unsigned char **) malloc(n_windows * sizeof(unsigned char *));
    unsigned char  *ml_final_state = (unsigned char  *) malloc(N_CHROM   * sizeof(unsigned char));
    for (i = 0; i < n_windows; i++)
        ml_prev_state[i] = (unsigned char *) malloc(N_STATES * sizeof(unsigned char));
    unsigned char  *ml_state_seq   = (unsigned char  *) malloc(n_windows * sizeof(unsigned char *));

    int start_window, end_window, delta;
    if (direction == FORWARD) {
        start_window = 0;
        end_window = n_windows;
        delta = 1;
    } else {
        start_window = n_windows-1;
        end_window = -1;
        delta = -1;
    }

    unsigned char last_chr = 0;
    int last_coord;

    double v_prev[N_STATES];  // viterbi state scores for previous window
    double v_cur[N_STATES];   // viterbi state scores for current window

    unsigned char tmp_state;
    double tmp_logp;
    for (i = start_window; i != end_window; i += delta) {
        if (max_cn[i] < 0) continue;

        // if at a new chromosome, restart the algorithm
        if (window_chr[i] != last_chr) {
            if (last_chr != 0) {
                tmp_state = DEL; tmp_logp = v_prev[DEL];
                if (v_prev[NORM] > tmp_logp) { tmp_state = NORM; tmp_logp = v_prev[NORM]; }
                if (v_prev[DUP]  > tmp_logp) { tmp_state = DUP;  tmp_logp = v_prev[DUP];  }
                ml_final_state[last_chr-1] = tmp_state;
            }
            last_chr = window_chr[i];
            if (direction == FORWARD)
                last_coord = window_start[i];
            else
                last_coord = window_end[i];
            v_prev[DEL]  = log(cnv_rate);
            v_prev[NORM] = log(1.0 - 2.0*cnv_rate);
            v_prev[DUP]  = log(cnv_rate);
        }

        v_cur[DEL]  = hmm_state_emission_logp[i][DEL];
        v_cur[NORM] = hmm_state_emission_logp[i][NORM];
        v_cur[DUP]  = hmm_state_emission_logp[i][DUP];

        // if we're not confident in the validity of the model at this window
        // don't give as much weight to its emission likelihoods
        // note that this doesn't bias in favor of the DIP state:
        // the effect of the window's predictions on all states is discounted
        v_cur[DEL]  *= model_conf[i];
        v_cur[NORM] *= model_conf[i];
        v_cur[DUP]  *= model_conf[i];

        // attenuation factor for state transitions
        // the probability of a CNV in the last window being extended to this one
        // is proportional to this factor
        double attenuation;
        if (direction == FORWARD)
            attenuation = exp(-((double)(window_start[i]-last_coord)) / mean_cnv_length);
        else
            attenuation = exp(-((double)(last_coord-window_end[i])) / mean_cnv_length);

        // find most likely previous state for DEL
        double del_del  = v_prev[DEL]  + log(transition_prob(DEL,  DEL, cnv_rate, attenuation));
        double norm_del = v_prev[NORM] + log(transition_prob(NORM, DEL, cnv_rate, attenuation));
        double dup_del  = v_prev[DUP]  + log(transition_prob(DUP,  DEL, cnv_rate, attenuation));
        tmp_state = DEL; tmp_logp = del_del;
        if (norm_del > tmp_logp) { tmp_state = NORM; tmp_logp = norm_del; }
        if (dup_del  > tmp_logp) { tmp_state = DUP;  tmp_logp = dup_del;  }
        ml_prev_state[i][DEL] = tmp_state;
        v_cur[DEL] += tmp_logp;

        // find most likely previous state for NORM
        double del_norm  = v_prev[DEL]  + log(transition_prob(DEL,  NORM, cnv_rate, attenuation));
        double norm_norm = v_prev[NORM] + log(transition_prob(NORM, NORM, cnv_rate, attenuation));
        double dup_norm  = v_prev[DUP]  + log(transition_prob(DUP,  NORM, cnv_rate, attenuation));
        tmp_state = DEL; tmp_logp = del_norm;
        if (norm_norm > tmp_logp) { tmp_state = NORM; tmp_logp = norm_norm; }
        if (dup_norm  > tmp_logp) { tmp_state = DUP;  tmp_logp = dup_norm;  }
        ml_prev_state[i][NORM] = tmp_state;
        v_cur[NORM] += tmp_logp;

        // find most likely previous state for DUP
        double del_dup  = v_prev[DEL]  + log(transition_prob(DEL,  DUP, cnv_rate, attenuation));
        double norm_dup = v_prev[NORM] + log(transition_prob(NORM, DUP, cnv_rate, attenuation));
        double dup_dup  = v_prev[DUP]  + log(transition_prob(DUP,  DUP, cnv_rate, attenuation));
        tmp_state = DEL; tmp_logp = del_dup;
        if (norm_dup > tmp_logp) { tmp_state = NORM; tmp_logp = norm_dup; }
        if (dup_dup  > tmp_logp) { tmp_state = DUP;  tmp_logp = dup_dup;  }
        ml_prev_state[i][DUP] = tmp_state;
        v_cur[DUP] += tmp_logp;

        // update lagging statistics
        if (direction == FORWARD)
            last_coord = window_start[i];
        else
            last_coord = window_end[i];
        v_prev[DEL]  = v_cur[DEL];
        v_prev[NORM] = v_cur[NORM];
        v_prev[DUP]  = v_cur[DUP];
    }

    tmp_state = DEL; tmp_logp = v_prev[DEL];
    if (v_prev[NORM] > tmp_logp) { tmp_state = NORM; tmp_logp = v_prev[NORM]; }
    if (v_prev[DUP]  > tmp_logp) { tmp_state = DUP;  tmp_logp = v_prev[DUP];  }
    ml_final_state[last_chr-1] = tmp_state;

    // backtrack through the DAG to find the maximum likelihood state sequence
    if (direction == FORWARD) {
        start_window = n_windows-1;
        end_window = -1;
        delta = -1;
    } else {
        start_window = 0;
        end_window = n_windows;
        delta = 1;
    }

    last_chr = 0;
    unsigned char last_state;
    int lookbehind = -delta;
    for (i = start_window; i != end_window; i += delta) {
        if (max_cn[i] < 0) { lookbehind -= delta; continue; };
        if (window_chr[i] != last_chr) {
            last_chr = window_chr[i];
            ml_state_seq[i] = ml_final_state[window_chr[i]-1];
        } else {
            ml_state_seq[i] = ml_prev_state[i+lookbehind][last_state];
        }

        last_state = ml_state_seq[i];
        lookbehind = -delta;
    }

    for (i = 0; i < n_windows; i++)
        free(ml_prev_state[i]);
    free(ml_prev_state);
    free(ml_final_state);
    return ml_state_seq;
}

void mask_sequence(int n_windows, char *max_cn,
                   unsigned char *seq1, unsigned char *seq2) {
    int i;
    for (i = 0; i < n_windows; i++) {
        if (max_cn[i] < 0) continue;
        if (seq1[i] != seq2[i]) seq1[i] = NORM;
    }
}

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
                      double **backward_scaled_prob) {
    int i, j;
    double attenuation;
    // emission likelihoods for the window (raw, not log)
    double E_del, E_norm, E_dup;
    // the forward/backward scaled prob for the last non-filtered window
    double prev_del, prev_norm, prev_dup;
    double scale_factor;

    unsigned char last_chr;
    int last_window;
    int last_coord;

    // compute forward posteriors
    last_chr = 0;
    for (i = 0; i < n_windows; i++) {
        if (max_cn[i] < 0) continue;

        E_del  = exp(hmm_state_emission_logp[i][DEL]  * model_conf[i]);
        E_norm = exp(hmm_state_emission_logp[i][NORM] * model_conf[i]);
        E_dup  = exp(hmm_state_emission_logp[i][DUP]  * model_conf[i]);

        if (window_chr[i] != last_chr) {
            last_chr = window_chr[i];
            forward_scaled_prob[i][DEL]  = E_del  * cnv_rate;
            forward_scaled_prob[i][NORM] = E_norm * (1.0 - 2.0*cnv_rate);
            forward_scaled_prob[i][DUP]  = E_dup  * cnv_rate;
        } else {
            attenuation = exp(-((double)(window_start[i]-last_coord)) / mean_cnv_length);

            forward_scaled_prob[i][DEL] =
                transition_prob(DEL,  DEL, cnv_rate, attenuation)  * prev_del;
            forward_scaled_prob[i][DEL] +=
                transition_prob(NORM, DEL, cnv_rate, attenuation)  * prev_norm;
            forward_scaled_prob[i][DEL] +=
                transition_prob(DUP,  DEL, cnv_rate, attenuation)  * prev_dup;
            forward_scaled_prob[i][DEL] *= E_del;

            forward_scaled_prob[i][NORM] =
                transition_prob(DEL,  NORM, cnv_rate, attenuation) * prev_del;
            forward_scaled_prob[i][NORM] +=
                transition_prob(NORM, NORM, cnv_rate, attenuation) * prev_norm;
            forward_scaled_prob[i][NORM] +=
                transition_prob(DUP,  NORM, cnv_rate, attenuation) * prev_dup;
            forward_scaled_prob[i][NORM] *= E_norm;

            forward_scaled_prob[i][DUP] =
                transition_prob(DEL,  DUP, cnv_rate, attenuation)  * prev_del;
            forward_scaled_prob[i][DUP] +=
                transition_prob(NORM, DUP, cnv_rate, attenuation)  * prev_norm;
            forward_scaled_prob[i][DUP] +=
                transition_prob(DUP,  DUP, cnv_rate, attenuation)  * prev_dup;
            forward_scaled_prob[i][DUP] *= E_dup;

        }

        scale_factor = forward_scaled_prob[i][DEL]  +
                       forward_scaled_prob[i][NORM] +
                       forward_scaled_prob[i][DUP];
        forward_scaled_prob[i][DEL]  /= scale_factor;
        forward_scaled_prob[i][NORM] /= scale_factor;
        forward_scaled_prob[i][DUP]  /= scale_factor;
 
        last_coord = window_start[i];
        prev_del   = forward_scaled_prob[i][DEL];
        prev_norm  = forward_scaled_prob[i][NORM];
        prev_dup   = forward_scaled_prob[i][DUP];
    }

    // compute backward posteriors
    last_chr = 0;
    last_window = n_windows;
    for (i = n_windows-1; i >= 0; i--) {
        if (max_cn[i] < 0) continue;

        if (window_chr[i] != last_chr) {
            last_chr = window_chr[i];
            backward_scaled_prob[i][DEL]  = 1.0;
            backward_scaled_prob[i][NORM] = 1.0;
            backward_scaled_prob[i][DUP]  = 1.0;
        } else {
            E_del  = exp(hmm_state_emission_logp[last_window][DEL]  * model_conf[last_window]);
            E_norm = exp(hmm_state_emission_logp[last_window][NORM] * model_conf[last_window]);
            E_dup  = exp(hmm_state_emission_logp[last_window][DUP]  * model_conf[last_window]);
            attenuation = exp(-((double)(last_coord-window_end[i])) / mean_cnv_length);

            backward_scaled_prob[i][DEL] =
                transition_prob(DEL, DEL,  cnv_rate, attenuation)  * E_del  * prev_del;
            backward_scaled_prob[i][DEL] +=
                transition_prob(DEL, NORM, cnv_rate, attenuation)  * E_norm * prev_norm;
            backward_scaled_prob[i][DEL] +=
                transition_prob(DEL, DUP,  cnv_rate, attenuation)  * E_dup  * prev_dup;

            backward_scaled_prob[i][NORM] =
                transition_prob(NORM, DEL,  cnv_rate, attenuation) * E_del  * prev_del;
            backward_scaled_prob[i][NORM] +=
                transition_prob(NORM, NORM, cnv_rate, attenuation) * E_norm * prev_norm;
            backward_scaled_prob[i][NORM] +=
                transition_prob(NORM, DUP,  cnv_rate, attenuation) * E_dup  * prev_dup;

            backward_scaled_prob[i][DUP] =
                transition_prob(DUP, DEL,  cnv_rate, attenuation)  * E_del  * prev_del;
            backward_scaled_prob[i][DUP] +=
                transition_prob(DUP, NORM, cnv_rate, attenuation)  * E_norm * prev_norm;
            backward_scaled_prob[i][DUP] +=
                transition_prob(DUP, DUP,  cnv_rate, attenuation)  * E_dup  * prev_dup;
        }

        scale_factor = backward_scaled_prob[i][DEL] +
                       backward_scaled_prob[i][NORM] +
                       backward_scaled_prob[i][DUP];
        backward_scaled_prob[i][DEL]  /= scale_factor;
        backward_scaled_prob[i][NORM] /= scale_factor;
        backward_scaled_prob[i][DUP]  /= scale_factor;

        last_coord  = window_end[i];
        last_window = i;
        prev_del    = backward_scaled_prob[i][DEL];
        prev_norm   = backward_scaled_prob[i][NORM];
        prev_dup    = backward_scaled_prob[i][DUP];
    }
}
