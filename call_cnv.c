#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "hmm.h"
#include "utils.h"

// I don't want to have to implement resizeable arrays in C,
// so I just make a big array to store CNVs
// and assume no sample will legitimately have more than this amount.
// (CLAMMS shouldn't be used with cancer data, it's not designed for that).
#define CNV_BUF_SIZE 4096

typedef struct {
    double cnv_rate;
    double mean_cnv_length;
    double min_gc;
    double max_gc;
    char sex;
} Options;

typedef struct {
    unsigned char chr;
    unsigned char type;
    unsigned char ml_copy_number;
    unsigned char max_considered_cn;
    int q_any;
    double model_fit;
    int n_windows;
    int start_window;
    int end_window;
    int start_coord;
    int end_coord;
    unsigned char can_extend_left;
    unsigned char can_extend_right;
    unsigned char can_contract_left;
    unsigned char can_contract_right;
    unsigned char q_extend_left;
    unsigned char q_extend_right;
    unsigned char q_contract_left;
    unsigned char q_contract_right;
    int extend_left_delta;
    int extend_right_delta;
    int contract_left_delta;
    int contract_right_delta;
} CNVCall;

Options parse_args(int argc, char *argv[], int arg_start) {
    Options options;
    options.cnv_rate = 3.0e-8;
    options.mean_cnv_length = 3.5e+4;
    options.min_gc = 0.3;
    options.max_gc = 0.7;
    options.sex = '\0';

    int i;
    for (i = arg_start; i < argc; i += 2) {
        if (strcmp(argv[i], "--cnv_rate") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.cnv_rate = strtod(argv[i+1], NULL);
            if (options.cnv_rate <= 0.0 || options.cnv_rate >= 0.1)
                invalid_value_error(argv[i]);
        } else if (strcmp(argv[i], "--mean_cnv_length") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.mean_cnv_length = strtod(argv[i+1], NULL);
            if(options.mean_cnv_length <= 0.0)
                invalid_value_error(argv[i]);
        } else if (strcmp(argv[i], "--min_gc") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.min_gc = strtod(argv[i+1], NULL);
            if (options.min_gc < 0.0 || options.min_gc > 1.0)
                invalid_value_error(argv[i]);
        } else if (strcmp(argv[i], "--max_gc") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.max_gc = strtod(argv[i+1], NULL);
            if (options.max_gc < 0.0 || options.max_gc > 1.0)
                invalid_value_error(argv[i]);
        } else if (strcmp(argv[i], "--sex") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            sscanf(argv[i+1], "%c", &options.sex);
            if (!(options.sex == 'M' || options.sex == 'F'))
                invalid_value_error(argv[i]);
        } else {
            fprintf(stderr, "Unrecognized argument: %s\n", argv[i]);
            fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
            exit(1);
        }
    }

    return options;
}

void note_cnv(unsigned char chr,
              unsigned char type,
              int n_windows,
              int start_window,
              int end_window,
              int start_coord,
              int end_coord,
              int *n_cnv,
              CNVCall *cnv) {
    if (*n_cnv >= CNV_BUF_SIZE) {
        fprintf(stderr, "ERROR: more than %d CNVs detected!\n", CNV_BUF_SIZE);
        fprintf(stderr, "The reference panel must not be a good fit for this sample.\n");
        fprintf(stderr, "Note: CLAMMS is not meant to be used with cancer samples.\n");
        exit(1);
    }

    cnv->chr = chr;
    cnv->type = type;
    cnv->n_windows = n_windows;
    cnv->start_window = start_window;
    cnv->end_window = end_window;
    cnv->start_coord = start_coord;
    cnv->end_coord = end_coord;
    (*n_cnv)++;
}

int call_cnv(int n_windows,
             unsigned char *window_chr,
             int *window_start,
             int *window_end,
             char *max_cn,
             unsigned char *ml_state_seq,
             CNVCall *cnv) {
    int i;
    int n_cnv = 0;
    unsigned char last_chr = 0;
    unsigned char last_state;
    int last_window;
    int last_end_coord;
    int state_start_coord;
    int state_start_window;
    int state_n_windows;
    for (i = 0; i < n_windows; i++) {
        if (max_cn[i] < 0) continue;
        if (window_chr[i] != last_chr) {
            if (last_chr != 0 && last_state != NORM) {
                note_cnv(last_chr, last_state,
                         state_n_windows, state_start_window, last_window,
                         state_start_coord, last_end_coord,
                         &n_cnv, cnv+n_cnv); }
            last_chr = window_chr[i];
            last_state = ml_state_seq[i];
            state_start_coord = window_start[i];
            state_start_window = i;
            state_n_windows = 1;
        } else if (ml_state_seq[i] != last_state) {
            if (last_state != NORM) {
                note_cnv(last_chr, last_state,
                         state_n_windows, state_start_window, last_window,
                         state_start_coord, last_end_coord,
                         &n_cnv, cnv+n_cnv); }
            last_state = ml_state_seq[i];
            state_start_coord = window_start[i];
            state_start_window = i;
            state_n_windows = 1;
        } else {
            state_n_windows++;
        }
        last_window = i;
        last_end_coord = window_end[i];
    }

    if (last_state != NORM) {
        note_cnv(last_chr, last_state,
                 state_n_windows, state_start_window, last_window,
                 state_start_coord, last_end_coord,
                 &n_cnv, cnv+n_cnv);
    }

    return n_cnv;
}

void estimate_het_or_hom(int n_cnv,
                         CNVCall *cnv,
                         char sex,
                         unsigned char *window_chr,
                         char *max_cn,
                         double **cn_emission_logp) {
    int i, j, k;
    for (i = 0; i < n_cnv; i++) {
        int norm_cn = expected_copy_number(sex, cnv[i].chr);
        if (norm_cn == HAPLOID) {
            cnv[i].max_considered_cn = 2;
            if (cnv[i].type == DEL)
                cnv[i].ml_copy_number = 0;
            else
                cnv[i].ml_copy_number = 2;
        } else {
            cnv[i].max_considered_cn = 3;
            for (j = cnv[i].start_window; j <= cnv[i].end_window; j++) {
                if (max_cn[j] > 3) { cnv[i].max_considered_cn = MAX_CN; break; }
            }
            if (cnv[i].type == DEL) {
                double logp_0 = 0.0;
                double logp_1 = 0.0;
                for (j = cnv[i].start_window; j <= cnv[i].end_window; j++) {
                    if (max_cn[j] < 0) continue;
                    logp_0 += cn_emission_logp[j][0];
                    logp_1 += cn_emission_logp[j][1];
                }
                if (logp_0 > logp_1)
                    cnv[i].ml_copy_number = 0;
                else
                    cnv[i].ml_copy_number = 1;
            } else {
                if (cnv[i].max_considered_cn == 3) {
                    cnv[i].ml_copy_number = 3;
                    continue;
                }

                double cn_logp[MAX_CN+1];
                for (k = 3; k <= MAX_CN; k++) cn_logp[k] = 0.0;
                for (j = cnv[i].start_window; j <= cnv[i].end_window; j++) {
                    if (max_cn[j] < 0) continue;
                    for (k = 3; k <= MAX_CN; k++)
                        cn_logp[k] += cn_emission_logp[j][k];
                }

                cnv[i].ml_copy_number = 3;
                double ml_cn_logp = cn_logp[3];
                for (k = 4; k <= MAX_CN; k++) {
                    if (cn_logp[k] > ml_cn_logp) {
                        cnv[i].ml_copy_number = k;
                        ml_cn_logp = cn_logp[k];
                    }
                }
            }
        }
    }
}

void calc_quality_metrics(int n_cnv,
                          CNVCall *cnv,
                          int n_windows,
                          unsigned char *window_chr,
                          int *window_start,
                          int *window_end,
                          char *max_cn,
                          double *cov,
                          unsigned char *hom_del_flag,
                          double *lambda,
                          double *sigma_dip,
                          double **hmm_state_emission_logp,
                          double cnv_rate,
                          double mean_cnv_length,
                          double **forward_scaled_prob,
                          double **backward_scaled_prob) {
    int i, j, k;
    double log_T_norm_norm = log(1.0 - 2.0*cnv_rate);
    double log_10 = log(10.0);
    double log_rad_4pi = log(sqrt(4.0 * M_PI));
    double log_rad_1_half   = log(sqrt(0.5));
    double log_rad_3_halves = log(sqrt(1.5));
    double log_rad_4_halves = log(sqrt(2.0));
    double log_rad_5_halves = log(sqrt(2.5));
    double log_rad_6_halves = log(sqrt(3.0));

    for (i = 0; i < n_cnv; i++) {
        int start_window = cnv[i].start_window;
        int end_window = cnv[i].end_window;
        unsigned char type = cnv[i].type;
        unsigned char mlcn = cnv[i].ml_copy_number;

        double log_p_norm = log(forward_scaled_prob[start_window][NORM]);
        for (j = start_window+1; j <= end_window; j++) {
            if (max_cn[j] < 0) continue;
            log_p_norm += hmm_state_emission_logp[j][NORM];
            log_p_norm += log_T_norm_norm;
        }

        // problem: suppose backwards scaled probs are 1/3 for DEL, DIP, DUP.
        // this means that future sequence gives no information about present state.
        // if the distribution were not uniform, it would act as a prior distribution.
        // to integrate this prior with the probability of the CNV region being all-diploid
        // would require also computing the marginal probabilities of all other possible
        // state sequences in the CNV region, which are exponentially many.
        // I can't figure out a way to do that, so I use this crappy heuristic instead.
        if (backward_scaled_prob[end_window][NORM] < 1.0 / N_STATES) {
            log_p_norm += log(backward_scaled_prob[end_window][NORM] * N_STATES);
        }

        cnv[i].q_any = (int) fmin(999.0, ((-10.0/log_10) * log_p_norm));

        // compute model goodness-of-fit metric
        // values < 1 fit model less well than would be expected
        // if you took random samples from the model.
        // values > 1 fit better.

        double log_fit_metric = 0.0;

        for (j = start_window; j <= end_window; j++) {
            if (max_cn[j] < 0) continue;
            if (mlcn == 0) {
                log_fit_metric += homozygous_del_log_likelihood(cov[j], hom_del_flag[j], lambda[j]);
                if (hom_del_flag[j])
                    log_fit_metric -= 3.912023; // log(50)
                else
                    log_fit_metric -= log(0.5 * lambda[j]);
            } else {
                log_fit_metric += gaussian_log_likelihood(cov[j], mlcn, sigma_dip[j]);
                log_fit_metric += log(sigma_dip[j]);
                log_fit_metric += log_rad_4pi;
            }
        }

        cnv[i].model_fit = log_fit_metric / cnv[i].n_windows;
             if (mlcn == 1) log_fit_metric += log_rad_1_half;
        else if (mlcn == 3) log_fit_metric += log_rad_3_halves;
        else if (mlcn == 4) log_fit_metric += log_rad_4_halves;
        else if (mlcn == 5) log_fit_metric += log_rad_5_halves;
        else if (mlcn == 6) log_fit_metric += log_rad_6_halves;

        // compute conditional call extension metrics
        // comparing L(call extended by 1 window) / L(call with specified breakpoints)

        int prev_window = get_prev_window(start_window, n_windows, window_chr, max_cn);
        if (prev_window == -1) {
            cnv[i].can_extend_left = 0;
        } else {
            cnv[i].can_extend_left = 1;
            double logp_ratio = 0.0;
            logp_ratio += log(forward_scaled_prob[prev_window][type]);
            logp_ratio -= log(forward_scaled_prob[prev_window][NORM]);
            double attenuation = exp(-((double)(window_start[start_window]-window_start[prev_window])) / mean_cnv_length);
            logp_ratio += log(transition_prob(type, type, cnv_rate, attenuation));
            logp_ratio -= log(transition_prob(NORM, type, cnv_rate, attenuation));
            cnv[i].q_extend_left = (unsigned char) fmin(99.0, fmax(0.0, (-10.0/log_10) * logp_ratio));
            cnv[i].extend_left_delta = window_start[prev_window] - window_start[start_window];
        }

        int next_window = get_next_window(end_window, n_windows, window_chr, max_cn);
        if (next_window == -1) {
            cnv[i].can_extend_right = 0;
        } else {
            cnv[i].can_extend_right = 1;
            double logp_ratio = 0.0;
            logp_ratio += log(backward_scaled_prob[next_window][type]);
            logp_ratio -= log(backward_scaled_prob[next_window][NORM]);
            logp_ratio += hmm_state_emission_logp[next_window][type];
            logp_ratio -= hmm_state_emission_logp[next_window][NORM];
            double attenuation = exp(-((double)(window_start[next_window]-window_start[end_window])) / mean_cnv_length);
            logp_ratio += log(transition_prob(type, type, cnv_rate, attenuation));
            logp_ratio -= log(transition_prob(type, NORM, cnv_rate, attenuation));
            cnv[i].q_extend_right = (unsigned char) fmin(99.0, fmax(0.0, (-10.0/log_10) * logp_ratio));
            cnv[i].extend_right_delta = window_end[next_window] - window_end[end_window];
        }

        // compute conditional call contraction metrics
        // comparing L(call contracted by 1 window) / L(call with specified breakpoints)
        cnv[i].can_contract_left = 0;
        cnv[i].can_contract_right = 0;
        if (cnv[i].n_windows > 1) {
            next_window = get_next_window(start_window, n_windows, window_chr, max_cn);
            if (next_window != -1) {
                cnv[i].can_contract_left = 1;
                double logp_ratio = 0.0;
                logp_ratio += log(forward_scaled_prob[start_window][NORM]);
                logp_ratio -= log(forward_scaled_prob[start_window][type]);
                double attenuation = exp(-((double)(window_start[next_window]-window_start[start_window])) / mean_cnv_length);
                logp_ratio += log(transition_prob(NORM, type, cnv_rate, attenuation));
                logp_ratio -= log(transition_prob(type, type, cnv_rate, attenuation));
                cnv[i].q_contract_left = (unsigned char) fmin(99.0, fmax(0.0, (-10.0/log_10) * logp_ratio));
                cnv[i].contract_left_delta = window_start[next_window] - window_start[start_window];
            }

            prev_window = get_prev_window(end_window, n_windows, window_chr, max_cn);
            if (prev_window != -1) {
                cnv[i].can_contract_right = 1;
                double logp_ratio = 0.0;
                logp_ratio += log(backward_scaled_prob[end_window][NORM]);
                logp_ratio -= log(backward_scaled_prob[end_window][type]);
                logp_ratio += hmm_state_emission_logp[end_window][NORM];
                logp_ratio -= hmm_state_emission_logp[end_window][type];
                double attenuation = exp(-((double)(window_end[end_window]-window_end[prev_window])) / mean_cnv_length);
                logp_ratio += log(transition_prob(type, NORM, cnv_rate, attenuation));
                logp_ratio -= log(transition_prob(type, type, cnv_rate, attenuation));
                cnv[i].q_contract_right = (unsigned char) fmin(99.0, fmax(0.0, (-10.0/log_10) * logp_ratio));
                cnv[i].contract_right_delta = window_end[prev_window] - window_end[end_window];
            }
        }
    }
}

void write_cnv(CNVCall *cnv,
               int n_cnv,
               char *sample_name) {
    int i;
    for (i = 0; i < n_cnv; i++) {
             if (cnv[i].chr == CHR_X) printf("X");
        else if (cnv[i].chr == CHR_Y) printf("Y");
        else if (cnv[i].chr == CHR_M) printf("MT");
        else                          printf("%hhu", cnv[i].chr);

        printf("\t%d\t%d\t",
                cnv[i].start_coord, cnv[i].end_coord);

             if (cnv[i].chr == CHR_X) printf("X");
        else if (cnv[i].chr == CHR_Y) printf("Y");
        else if (cnv[i].chr == CHR_M) printf("MT");
        else                          printf("%hhu", cnv[i].chr);

        printf(":%d-%d\t%s\t%s\t%hhu\t%d\t%d\t%.3g",
                cnv[i].start_coord, cnv[i].end_coord,
                sample_name,
                (cnv[i].type == DEL ? "DEL" : "DUP"),
                cnv[i].ml_copy_number,
                cnv[i].n_windows,
                cnv[i].q_any,
                cnv[i].model_fit);

        if (cnv[i].can_extend_left)
            printf("\t%hhu\t%d", cnv[i].q_extend_left, cnv[i].extend_left_delta);
        else
            printf("\tNA\tNA");

        if (cnv[i].can_extend_right)
            printf("\t%hhu\t%d", cnv[i].q_extend_right, cnv[i].extend_right_delta);
        else
            printf("\tNA\tNA");

        if (cnv[i].can_contract_left)
            printf("\t%hhu\t%d", cnv[i].q_contract_left, cnv[i].contract_left_delta);
        else
            printf("\tNA\tNA");

        if (cnv[i].can_contract_right)
            printf("\t%hhu\t%d", cnv[i].q_contract_right, cnv[i].contract_right_delta);
        else
            printf("\tNA\tNA");

        printf("\n");
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s sample.norm.cov.bed models.out [OPTIONS] >sample.cnvs.bed\n\n", argv[0]);
        fputs("Calls CNVs for a sample.\n", stderr);
        fputs("sample.norm.cov.bed should have been generated by the CLAMMS 'normalize_coverage' program.\n", stderr);
        fputs("models.out should have been generated by the CLAMMS 'fit_models' program.\n\n", stderr);
        fputs("If you want to make calls on sex chromosomes, you must specify\n", stderr);
        fputs("a sex for this sample with the --sex option.\n", stderr);
        fputs("Additionally, when 'fit_models' was run to generate models.out,\n", stderr);
        fputs("you must have specified sexes for each sample.\n\n", stderr);
        fputs("  --cnv_rate           P(DIP->DIP) == P(DIP->DUP)\n", stderr);
        fputs("                       Default = 3.0e-8.\n", stderr);
        fputs("  --mean_cnv_length    Mean of prior distribution for CNV lengths (in b.p.)\n", stderr);
        fputs("                       Default = 3.5e+4.\n", stderr);
        fputs("  --min_gc             If you used a non-default --min_gc for 'normalize_coverage' and 'fit_models'\n", stderr);
        fputs("                       then you must use it again here.\n", stderr);
        fputs("  --max_gc             If you used a non-default --max_gc for 'normalize_coverage' and 'fit_models'\n", stderr);
        fputs("                       then you must use it again here.\n", stderr);
        fputs("  --sex                'M' or 'F'\n", stderr);
        fputs("                       Default = unspecified (no calls on sex chr).\n\n", stderr);
        return 1;
    }

    FILE *coverage = open_file(argv[1]);
    FILE *models   = open_file(argv[2]);

    char sample_name[1024];
    read_sample_name(sample_name, argv[1]);

    Options options = parse_args(argc, argv, 3);

    int i;
    int n_windows = count_lines_in_file(coverage);

    // read all the coverage values and mixture model parameters into memory

    unsigned char *window_chr = (unsigned char *) malloc(n_windows * sizeof(unsigned char));
    int *window_start = (int *) malloc(n_windows * sizeof(int));
    int *window_end = (int *) malloc(n_windows * sizeof(int));
    char *max_cn = (char *) malloc(n_windows * sizeof(char));
    unsigned char *hom_del_flag  = (unsigned char *) malloc(n_windows * sizeof(unsigned char));

    double *window_gc   = (double *) malloc(n_windows * sizeof(double));
    double *cov         = (double *) malloc(n_windows * sizeof(double));
    double *lambda      = (double *) malloc(n_windows * sizeof(double));
    double *mu_dip      = (double *) malloc(n_windows * sizeof(double));
    double *sigma_dip   = (double *) malloc(n_windows * sizeof(double));
    double *model_conf  = (double *) malloc(n_windows * sizeof(double));

    read_model_data(models, n_windows,
                    window_chr, window_start, window_end,
                    max_cn, hom_del_flag, window_gc,
                    lambda, mu_dip, sigma_dip, model_conf);
    read_coverage_data(coverage, n_windows,
                       window_chr, window_start, window_end,
                       cov, mu_dip);
    calc_base_model_conf(n_windows, options.min_gc, options.max_gc,
                         window_chr, window_start, window_end,
                         max_cn, window_gc, model_conf);
    calc_sample_specific_model_conf(n_windows, options.sex, window_chr,
                                    max_cn, cov, hom_del_flag,
                                    lambda, sigma_dip, model_conf);

    free(window_gc);
    free(mu_dip);

    // if the sex of the sample isn't specified, don't make calls for chrX/Y
    if (!options.sex) {
        for (i = 0; i < n_windows; i++) {
            if (window_chr[i] == CHR_X || window_chr[i] == CHR_Y)
                max_cn[i] = -1;
        }
    }

    // calculate the emission log-probabilities for each copy number state
    // and for each HMM state (DEL, DIP, DUP) using Bayes theorem with a uniform prior. 
    //
    //                                     likelihood           uniform prior
    //                                        /                       \
    //                                 P(coverage | CN=2)   *          1
    // example: P(DIP | coverage) = -------------------------------------------------
    //                              sum {k in 0..MAX_CN} P(coverage | CN=k) * 1
    //                                              \
    //                                        normalizing factor ("evidence")
    //

    double **cn_emission_logp        = (double **) malloc(n_windows * sizeof(double *));
    double **hmm_state_emission_logp = (double **) malloc(n_windows * sizeof(double *));
    for (i = 0; i < n_windows; i++) {
        cn_emission_logp[i]        = (double *) malloc((MAX_CN+1) * sizeof(double));
        hmm_state_emission_logp[i] = (double *) malloc( N_STATES  * sizeof(double));
    }

    calc_cn_emission_logp(n_windows, options.sex,
                          window_chr, max_cn, cov, 
                          hom_del_flag, lambda, sigma_dip,
                          cn_emission_logp);
    calc_hmm_state_emission_logp(n_windows, options.sex,
                                 window_chr, max_cn,
                                 cn_emission_logp, hmm_state_emission_logp);

    // We run the Viterbi algorithm in both directions
    // and only call variants when both Viterbi runs predict non-diploid state.
    // This avoids the directionality bias that comes from the transition model
    // where it costs a lot to open an CNV but not much to extend it.
    // It is a bit conservative however, so we increased the default value
    // for p (the transition probability of NORM to DEL or DUP) to compensate.
    unsigned char *ml_seq = viterbi(
        n_windows, FORWARD,
        window_chr, window_start, window_end,
        max_cn, model_conf,
        hmm_state_emission_logp,
        options.cnv_rate, options.mean_cnv_length
    );
    unsigned char *backward_seq = viterbi(
        n_windows, BACKWARD,
        window_chr, window_start, window_end,
        max_cn, model_conf,
        hmm_state_emission_logp,
        options.cnv_rate, options.mean_cnv_length
    );
    mask_sequence(n_windows, max_cn, ml_seq, backward_seq);

    // identify CNVs in the consensus state sequence
    CNVCall *cnv = (CNVCall *) malloc(CNV_BUF_SIZE * sizeof(CNVCall));
    int n_cnv = call_cnv(n_windows,
                         window_chr, window_start, window_end,
                         max_cn, ml_seq, cnv);
    free(ml_seq);
    free(backward_seq);

    // Run the forward-backward algorithm
    // We will use the posterior probabilities it outputs to compute quality metrics
    // such as P(any CNV in called regions) for putative CNV calls
    double **forward_scaled_prob  = (double **) malloc(n_windows * sizeof(double *));
    double **backward_scaled_prob = (double **) malloc(n_windows * sizeof(double *));
    for (i = 0; i < n_windows; i++) {
        forward_scaled_prob[i]  = (double *) malloc(N_STATES * sizeof(double));
        backward_scaled_prob[i] = (double *) malloc(N_STATES * sizeof(double));
    }

    forward_backward(
        n_windows,
        window_chr, window_start, window_end,
        max_cn, model_conf,
        hmm_state_emission_logp,
        options.cnv_rate, options.mean_cnv_length,
        forward_scaled_prob, backward_scaled_prob);
    estimate_het_or_hom(
        n_cnv, cnv, options.sex,
        window_chr, max_cn, cn_emission_logp);
    calc_quality_metrics(
        n_cnv, cnv, n_windows,
        window_chr, window_start, window_end,
        max_cn, cov,
        hom_del_flag, lambda, sigma_dip,
        hmm_state_emission_logp,
        options.cnv_rate, options.mean_cnv_length,
        forward_scaled_prob, backward_scaled_prob);
    
    write_cnv(cnv, n_cnv, sample_name);

    free(window_chr);
    free(window_start);
    free(window_end);
    free(model_conf);
    free(cov);
    free(max_cn);
    free(hom_del_flag);
    free(lambda);
    free(sigma_dip);

    for (i = 0; i < n_windows; i++) {
        free(cn_emission_logp[i]);
        free(hmm_state_emission_logp[i]);
        free(forward_scaled_prob[i]);
        free(backward_scaled_prob[i]);
    }

    free(cn_emission_logp);
    free(hmm_state_emission_logp);
    free(forward_scaled_prob);
    free(backward_scaled_prob);
    free(cnv);

    fclose(coverage);
    fclose(models);

    return 0;
}

