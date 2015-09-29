#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "utils.h"

#define MIN_LAMBDA                              16.0 // lambda is a model parameter; see the MixtureModel struct
#define MAX_EM_ITS                                30 // maximum # of EM iterations to run per window
#define MIN_REL_PARAM_SHIFT                   0.0025 // threshold for early convergence of the EM algorithm
                                                     // if no model parameter shifts by more than this*100%
                                                     // in an iteration, stop iterating
#define GAUSSIAN_OUTLIER_THRESHOLD               2.5 // standard deviations
#define EXPONENTIAL_OUTLIER_THRESHOLD       4.388501 // P(std. exponential <= this) == P(abs(std. normal) <= 2.5-sigma)
#define MAD_ASYMPTOTIC_NORMALITY_SCALE_FACTOR 1.4826 // c.f. the R implementation of 'mad'

typedef struct {
    double min_dip_cov;
    double max_sigma;
    double min_gc;
    double max_gc;
    double min_mappability;
} Options;

typedef struct {
    int max_cn;
    int homozyg_del_is_point;  // if 1, homozygous del is a Dirac distribution (single point) @ -1
                               // if 0, homozygous del is an exponential distribution
                               //       shifted to start at -1
    double lambda;             // parameter for homozygous del exponential distribution
    double mu_dip;             // mean of the component corresponding to CN=2 (diploid)
                               // other gaussian means are constrainted to be a half-integer multiple of this
    double sigma_dip;          // standard deviation of the component corresponding to CN=2 (diploid)
                               // other gaussian standard deviations are fixed ratios of this,
                               // based on empirical estimates of the over/under-dispersion
                               // of the coverage distributions for haploid / duplicated genomic regions
    double *component_weights;
    double *male_comp_weights; // for sex chr, seperate component weights are used for males vs. females
} MixtureModel;

Options parse_args(int argc, char *argv[], int arg_start) {
    Options options;
    options.min_dip_cov     = 0.1;
    options.max_sigma       = 0.25;
    options.min_gc          = 0.3;
    options.max_gc          = 0.7;
    options.min_mappability = 0.75;

    int i;
    for (i = arg_start; i < argc; i += 2) {
        if (strcmp(argv[i], "--min_dip_cov") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.min_dip_cov = strtod(argv[i+1], NULL);
            if (options.min_dip_cov < 0.0)
                invalid_value_error(argv[i]);
        } else if (strcmp(argv[i], "--max_sigma") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.max_sigma = strtod(argv[i+1], NULL);
            if (options.max_sigma <= 0.0)
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
        } else if (strcmp(argv[i], "--min_mappability") == 0) {
            if (i+1 >= argc) missing_value_error(argv[i]);
            options.min_mappability = strtod(argv[i+1], NULL);
            if (options.min_mappability < 0.0 || options.min_mappability > 1.0)
                invalid_value_error(argv[i]);
        } else {
            fprintf(stderr, "Unrecognized argument: %s\n", argv[i]);
            fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
            exit(1);
        }
    }

    return options;
}

int em_iteration(MixtureModel *model,
                 int n_samples,
                 int is_sex_chr,
                 char *sample_sexes,
                 double male_proportion,
                 int *expected_cn,
                 double *raw_cov,  // original values
                 double *cov,      // centered and scaled around mu_dip
                 double **likelihoods,
                 double **soft_assignments,
                 int *is_outlier) {

    int i; // sample index
    int k; // component index
    double *comp_likelihoods; // pointer to current component's likelihood array
    double *comp_assignments; // pointer to current component's soft assignment array
    double min_comp_weight = 0.1 / n_samples; // if a component has weight less than this, skip it

    // this variable will store the largest relative change (abs(new-old)/old)
    // made this iteration for mu_dip, sigma_dip, and the component weights.
    // if the maximum relative change is small, we'll return 1
    // to signal that the EM iterations should be terminated early.
    double max_rel_param_shift = 0.0;

    // compute likelihoods of CN = 0
    comp_likelihoods = likelihoods[0];
    if (model->component_weights[0] < min_comp_weight
        && (!is_sex_chr || model->male_comp_weights[0] < min_comp_weight)) {
        for (i = 0; i < n_samples; i++)
            comp_likelihoods[i] = 0.0;
    } else if (model->homozyg_del_is_point) {
        for (i = 0; i < n_samples; i++) {
            if (is_outlier[i]) continue;
            if (cov[i] < HOM_DEL_THRESHOLD)
                comp_likelihoods[i] = 50.0;
            else
                comp_likelihoods[i] = 0.0;
        }
    } else {
        double log_lambda = log(model->lambda);
        double log_likelihood;
        for (i = 0; i < n_samples; i++) {
            if (is_outlier[i]) continue;
            log_likelihood = log_lambda - model->lambda * (cov[i] + 1.0);
            if (log_likelihood < -100.0)
                comp_likelihoods[i] = 0.0;
            else
                comp_likelihoods[i] = exp(log_likelihood);
        }
    }

    // compute likelihoods of other CN values
    for (k = 1; k <= model->max_cn; k++) {
        comp_likelihoods = likelihoods[k];

        if (model->component_weights[k] < min_comp_weight
            && (!is_sex_chr || model->male_comp_weights[k] < min_comp_weight)) {
            for (i = 0; i < n_samples; i++)
                comp_likelihoods[i] = 0.0;
            continue;
        }

        double mu = -1.0 + 0.5 * k;
        double sigma = model->sigma_dip;
             if (k == 1)                      sigma *= SIGMA_RATIO_CN1;
        else if (k == 3 && model->max_cn > 3) sigma *= SIGMA_RATIO_CN3;
        else if (k == 4 && model->max_cn > 3) sigma *= SIGMA_RATIO_CN4;
        double sigma_sq = sigma * sigma;
        double const_term = -0.5 * log(M_2_PI * sigma_sq);
        double scale_term = 1.0 / (2.0 * sigma_sq);
        double dev, log_likelihood;
        for (i = 0; i < n_samples; i++) {
            if (is_outlier[i]) continue;
            if (k > 2*expected_cn[i]) {
                comp_likelihoods[i] = 0.0;
                continue;
            }
            
            dev = cov[i] - mu;
            log_likelihood = const_term - scale_term * dev * dev;
            if (log_likelihood < -100.0)
                comp_likelihoods[i] = 0.0;
            else
                comp_likelihoods[i] = exp(log_likelihood);
        }
    }

    // compute soft assignments and update cluster weights

    double total_weighted_likelihood[n_samples];
    for (i = 0; i < n_samples; i++)
        total_weighted_likelihood[i] = 0.0;
    for (k = 0; k <= model->max_cn; k++) {
        comp_likelihoods = likelihoods[k];
        if (is_sex_chr) {
            for (i = 0; i < n_samples; i++) {
                if (sample_sexes[i] == 'F') {
                    total_weighted_likelihood[i] += model->component_weights[k] * comp_likelihoods[i];
                } else {
                    total_weighted_likelihood[i] += model->male_comp_weights[k] * comp_likelihoods[i];
                }
            }
        } else {
            for (i = 0; i < n_samples; i++) {
                total_weighted_likelihood[i] += model->component_weights[k] * comp_likelihoods[i];
            }
        }
    }

    double total_component_assignment[model->max_cn+1];
    double total_assignment = 0.0;
    double old_max_comp_weight = 0.0;
    double new_max_comp_weight = 0.0;
    if (is_sex_chr) {
        double tot_f_comp_assignment[model->max_cn+1];
        double tot_m_comp_assignment[model->max_cn+1];
        double tot_f_assignment = 0.0;
        double tot_m_assignment = 0.0;
        for (k = 0; k <= model->max_cn; k++) {
            tot_f_comp_assignment[k] = 0.0;
            tot_m_comp_assignment[k] = 0.0;
            comp_likelihoods = likelihoods[k];
            comp_assignments = soft_assignments[k];

            double score;
            for (i = 0; i < n_samples; i++) {
                if (total_weighted_likelihood[i] == 0.0) {
                    score = 0.0;
                } else if (sample_sexes[i] == 'F') {
                    score = model->component_weights[k] * comp_likelihoods[i] / total_weighted_likelihood[i];
                    tot_f_comp_assignment[k] += score;
                } else {
                    score = model->male_comp_weights[k] * comp_likelihoods[i] / total_weighted_likelihood[i];
                    tot_m_comp_assignment[k] += score;
                }
                comp_assignments[i] = score;
            }

            total_component_assignment[k] = tot_f_comp_assignment[k] + tot_m_comp_assignment[k];
            tot_f_assignment += tot_f_comp_assignment[k];
            tot_m_assignment += tot_m_comp_assignment[k];
        }

        total_assignment = tot_f_assignment + tot_m_assignment;

        for (k = 0; k <= model->max_cn; k++) {
            old_max_comp_weight = fmax(model->component_weights[k], old_max_comp_weight);
            old_max_comp_weight = fmax(model->male_comp_weights[k], old_max_comp_weight);

            model->component_weights[k] = tot_f_comp_assignment[k] / tot_f_assignment;
            model->male_comp_weights[k] = tot_m_comp_assignment[k] / tot_m_assignment;

            new_max_comp_weight = fmax(model->component_weights[k], new_max_comp_weight);
            new_max_comp_weight = fmax(model->male_comp_weights[k], new_max_comp_weight);
        }
    } else {
        for (k = 0; k <= model->max_cn; k++) {
            total_component_assignment[k] = 0.0;
            comp_likelihoods = likelihoods[k];
            comp_assignments = soft_assignments[k];
            double score;
            for (i = 0; i < n_samples; i++) {
                if (total_weighted_likelihood[i] == 0.0)
                    score = 0.0;
                else
                    score = model->component_weights[k] * comp_likelihoods[i] / total_weighted_likelihood[i];
                comp_assignments[i] = score;
                total_component_assignment[k] += score;
            }

            total_assignment += total_component_assignment[k];
        }

        double old_max_comp_weight = 0.0;
        double new_max_comp_weight = 0.0;
        for (k = 0; k <= model->max_cn; k++) {
            old_max_comp_weight = fmax(model->component_weights[k], old_max_comp_weight);
            model->component_weights[k] = total_component_assignment[k] / total_assignment;
            new_max_comp_weight = fmax(model->component_weights[k], new_max_comp_weight);
        }
    }

    max_rel_param_shift = fmax(
        fabs(new_max_comp_weight - old_max_comp_weight) / old_max_comp_weight,
        max_rel_param_shift
    );

    // update homozygous deletion component

    if (model->component_weights[0] < min_comp_weight
        && (!is_sex_chr || model->male_comp_weights[0] < min_comp_weight)) {
        model->lambda = 0.0;
        model->homozyg_del_is_point = 1;
    } else if (!model->homozyg_del_is_point) {
        double mu_hat = 0.0;
        comp_assignments = soft_assignments[0];
        for (i = 0; i < n_samples; i++)
            mu_hat += comp_assignments[i] * cov[i];
        mu_hat /= total_component_assignment[0];

        if (mu_hat < -0.995) {
            model->lambda = 0.0;
            model->homozyg_del_is_point = 1;
        } else if (mu_hat > -1.0 + 1.0/MIN_LAMBDA) {
            model->lambda = MIN_LAMBDA;
        } else {
            model->lambda = 1.0 / (mu_hat + 1.0);
        }
    }

    // update mu_dip
    double old_mu_dip = model->mu_dip;
    double new_mu_dip = 0.0;
    for (k = 1; k <= model->max_cn; k++) {
        if (model->component_weights[k] < min_comp_weight
            && (!is_sex_chr || model->male_comp_weights[k] < min_comp_weight)) {
            continue;
        }
        comp_assignments = soft_assignments[k];
        double mu_hat = 0.0; // estimate of this component's mean
                             // in the centered+scaled coordinates around the old mu_dip
        for (i = 0; i < n_samples; i++)
            mu_hat += comp_assignments[i] * cov[i];
        mu_hat /= total_component_assignment[k];
        if (is_sex_chr) {
            new_mu_dip += ((1.0 - male_proportion) * model->component_weights[k]
                                + male_proportion  * model->male_comp_weights[k])
                          * (2.0/k) * (1.0+mu_hat) * model->mu_dip;
        } else {
            new_mu_dip += model->component_weights[k] * (2.0/k) * (1.0+mu_hat) * model->mu_dip;
        }
    }
    if (is_sex_chr) {
        new_mu_dip /= (1.0 - ((1.0 - male_proportion) * model->component_weights[0]
                                   + male_proportion  * model->male_comp_weights[0]));
    } else {
        new_mu_dip /= (1.0 - model->component_weights[0]);
    }
    max_rel_param_shift = fmax(
        fabs(new_mu_dip - model->mu_dip) / model->mu_dip,
        max_rel_param_shift);
    model->mu_dip = new_mu_dip;

    // re-center/scale coverage values
    for (i = 0; i < n_samples; i++)
        cov[i] = (raw_cov[i] - model->mu_dip) / model->mu_dip;

    // update sigma_dip
    double new_var = 0.0;
    for (k = 1; k <= model->max_cn; k++) {
        if (model->component_weights[k] < min_comp_weight
            && (!is_sex_chr || model->male_comp_weights[k] < min_comp_weight)) continue;
        comp_assignments = soft_assignments[k];
        double mu = -1.0 + 0.5 * k;
        double var_hat = 0.0;
        for (i = 0; i < n_samples; i++) {
            double dev = cov[i] - mu;
            var_hat += comp_assignments[i] * dev * dev;
        }
        var_hat /= total_component_assignment[k];

        if (k == 1)
            var_hat /= (SIGMA_RATIO_CN1 * SIGMA_RATIO_CN1);
        else if (k == 3 && model->max_cn > 3)
            var_hat /= (SIGMA_RATIO_CN3 * SIGMA_RATIO_CN3);
        else if (k == 4 && model->max_cn > 3)
            var_hat /= (SIGMA_RATIO_CN4 * SIGMA_RATIO_CN4);

        if (is_sex_chr) {
            new_var += ((1.0 - male_proportion) * model->component_weights[k]
                             + male_proportion  * model->male_comp_weights[k])
                       * var_hat;
        } else {
            new_var += model->component_weights[k] * var_hat;
        }
    }
    if (is_sex_chr) {
        new_var /= (1.0 - ((1.0 - male_proportion) * model->component_weights[0]
                                + male_proportion  * model->male_comp_weights[0]));
    } else {
        new_var /= (1.0 - model->component_weights[0]);
    }
    double old_var = model->sigma_dip * model->sigma_dip;
    max_rel_param_shift = fmax(fabs(new_var - old_var) / old_var, max_rel_param_shift);
    model->sigma_dip = sqrt(new_var);

    if (max_rel_param_shift < MIN_REL_PARAM_SHIFT)
        return 1;
    else
        return 0;
}

void filter_window(char *window_coords,
                   double window_gc,
                   double window_mappability,
                   double mu_dip,
                   double sigma_dip) {
    int i;
    printf("%s\t-1\t%.3f\t%.3f\t1\t%.4g\t%.4g\t%.4g",
            window_coords, window_gc, window_mappability,
            0.0, mu_dip, sigma_dip);
    for (i = 0; i <= MAX_CN; i++)
        printf("\t%.1f", 0.0);
    printf("\n");
}

// contract: if it is a sex chromosome,
//           sample sexes will be provided
void process_window(char *window_coords,
                    int max_cn,
                    double window_gc,
                    double window_mappability,
                    int n_samples,
                    int n_male,
                    int n_female,
                    char *sample_sexes,
                    int *expected_cn,
                    Options *options,
                    double *cov,
                    double *cov_copy,
                    int *is_outlier,
                    double *component_weights,
                    double *male_comp_weights,
                    double **likelihoods,
                    double **soft_assignments,
                    int *total_n_filtered_windows,
                    int *total_n_outlier_points) {

    int i, j;
    int is_sex_chr = (*window_coords == 'X' || *window_coords == 'Y');

    MixtureModel model;
    model.max_cn = max_cn;
    model.homozyg_del_is_point = 0;
    model.lambda = MIN_LAMBDA;

    // initialize expected_cn for each sample
    if (!sample_sexes) {
        for (i = 0; i < n_samples; i++)
            expected_cn[i] = 2;
    } else {
        for (i = 0; i < n_samples; i++) {
            if (sample_sexes[i] == 'F' && *window_coords == 'Y')
                expected_cn[i] = 0;
            else if ((sample_sexes[i] == 'M' && is_sex_chr) || *window_coords == 'M') {
                expected_cn[i] = 1;
            } else {
                expected_cn[i] = 2;
            }
        }
    }

    double male_proportion;
    if (sample_sexes) male_proportion = (double) n_male / (double) n_samples;

    // initialize model parameter mu_dip (mean coverage of diploid sampless)
    double raw_cov_median;
    memcpy(cov_copy, cov, n_samples * sizeof(double));
    for (i = 0; i < n_samples; i++) {
        if (expected_cn[i] == 0)
            cov_copy[i] = -1.0; // flag
        else if (expected_cn[i] == 1)
            cov_copy[i] *= 2.0;
    }
    qsort(cov_copy, n_samples, sizeof(double), double_comp);
    if (*window_coords == 'Y')
        raw_cov_median = median(cov_copy+n_female, n_male);
    else
        raw_cov_median = median(cov_copy, n_samples);
    model.mu_dip = raw_cov_median;

    // initialize model parameter sigma_dip (standard deviation of gaussian components)
    // note that the centered+scaled coverage values that we compute
    // are always scaled relative to mu_dip, even if the sample is expected to be haploid
    double raw_cov_mad;
    for (i = 0; i < n_samples; i++) {
        if (expected_cn[i] == 0)
            cov_copy[i] = -1.0; // flag
        else if (expected_cn[i] == 1)
            cov_copy[i] = fabs(cov[i] - (model.mu_dip/2.0)) / model.mu_dip;
        else
            cov_copy[i] = fabs(cov[i] - model.mu_dip) / model.mu_dip;
    }
    qsort(cov_copy, n_samples, sizeof(double), double_comp);
    if (*window_coords == 'Y')
        raw_cov_mad = median(cov_copy+n_female, n_male) * MAD_ASYMPTOTIC_NORMALITY_SCALE_FACTOR;
    else
        raw_cov_mad = median(cov_copy, n_samples) * MAD_ASYMPTOTIC_NORMALITY_SCALE_FACTOR;
    model.sigma_dip = raw_cov_mad;

    // if median coverage is too low, don't try to fit a model here
    if (model.mu_dip < options->min_dip_cov) {
        (*total_n_filtered_windows)++;
        if (model.sigma_dip != model.sigma_dip) model.sigma_dip = 0.0;
        filter_window(window_coords, window_gc, window_mappability, model.mu_dip, model.sigma_dip);
        return;
    }

    // center and scale coverage values around mu_dip
    for (i = 0; i < n_samples; i++)
        cov_copy[i] = (cov[i] - model.mu_dip) / model.mu_dip;

    // initialize component weights
    // including a separate set of weights for males if we are in a sex chr
    model.component_weights = component_weights;
    for (i = 0; i <= model.max_cn; i++) {
        model.component_weights[i] = 1.0 / ((double) (model.max_cn + 1));
        for (j = 0; j < n_samples; j++) {
            likelihoods[i][j] = 0.0;
            soft_assignments[i][j] = 0.0;
        }
    }

    if (is_sex_chr) {
        model.male_comp_weights = male_comp_weights;
        for (i = 0; i <= model.max_cn; i++) {
            model.male_comp_weights[i] = 1.0 / ((double) (model.max_cn + 1));
        }
    }

    int n_outliers = 0;
    for (i = 0; i < n_samples; i++)
        is_outlier[i] = 0;

    // fit models parameters using EM algorithm
    int stop = 0;
    for (i = 0; i < MAX_EM_ITS; i++) {
        // cov is raw values, cov_copy is centered/scaled values
        // em_iteration will mutate cov_copy to update it when mu_dip is updated
        stop = em_iteration(&model, n_samples, is_sex_chr,
                            sample_sexes, male_proportion, expected_cn,
                            cov, cov_copy, likelihoods, soft_assignments,
                            is_outlier);
        if (stop) break;
    }

    // flag outlier coverage values
    for (i = 0; i < n_samples; i++) {
        if (expected_cn[i] == 0) continue;

        int flag = 1;
        if (cov_copy[i] <= -0.9999) {
            flag = 0;
        } else if (!model.homozyg_del_is_point
                   && (cov_copy[i] + 1.0) * model.lambda <= EXPONENTIAL_OUTLIER_THRESHOLD) {
            flag = 0;
        }

        for (j = 1; j <= model.max_cn; j++) {
            if (j > 2*expected_cn[i]) break;
            double mu = -1.0 + 0.5 * j;
            double sigma = model.sigma_dip;
                 if (j == 1)                     sigma *= SIGMA_RATIO_CN1;
            else if (j == 3 && model.max_cn > 3) sigma *= SIGMA_RATIO_CN3;
            else if (j == 4 && model.max_cn > 3) sigma *= SIGMA_RATIO_CN4;
            if (fabs((cov_copy[i] - mu) / sigma) <= GAUSSIAN_OUTLIER_THRESHOLD) {
                flag = 0;
                break;
            }
        }

        if (flag) {
            is_outlier[i] = 1;
            n_outliers++;
        }
    }

    if (n_outliers > 0) {
        if (is_sex_chr) {
            int n_non_outlier_samples = n_samples - n_outliers;
            int n_non_outlier_males = 0;
            for (i = 0; i < n_samples; i++) {
                if (!is_outlier[i] && sample_sexes[i] == 'M')
                    n_non_outlier_males++;
            }
            male_proportion = (double) n_non_outlier_males / (double) n_non_outlier_samples;
        }

        model.homozyg_del_is_point = 0;
        model.lambda = MIN_LAMBDA;
        model.mu_dip = raw_cov_median;
        model.sigma_dip = raw_cov_mad;

        for (i = 0; i < n_samples; i++)
            cov_copy[i] = (cov[i] - model.mu_dip) / model.mu_dip;

        for (i = 0; i <= model.max_cn; i++) {
            model.component_weights[i] = 1.0 / ((double) (model.max_cn + 1));
            for (j = 0; j < n_samples; j++) {
                likelihoods[i][j] = 0.0;
                soft_assignments[i][j] = 0.0;
            }
        }

        stop = 0;
        for (i = 0; i < MAX_EM_ITS; i++) {
            stop = em_iteration(&model, n_samples,
                                is_sex_chr, sample_sexes, male_proportion, expected_cn,
                                cov, cov_copy, likelihoods, soft_assignments,
                                is_outlier);
            if (stop) break;
        }
    }

    // filter windows in problematic regions
    if (model.mu_dip < options->min_dip_cov || model.sigma_dip > options->max_sigma) {
        model.max_cn = -1;
        (*total_n_filtered_windows)++;
    } else if (model.component_weights[1] >= 0.5 &&
               model.component_weights[2] <= 0.2 &&
               *window_coords != 'M') {
        model.max_cn = -1;
        (*total_n_filtered_windows)++;
    } else if (n_outliers > 0) {
        (*total_n_outlier_points) += n_outliers;
    }

    double eff_sample_size[model.max_cn+1];
    for (i = 0; i <= model.max_cn; i++) {
        if (is_sex_chr) {
            eff_sample_size[i] = (n_samples - n_outliers) *
                                        ((1.0 - male_proportion) * model.component_weights[i]
                                              + male_proportion  * model.male_comp_weights[i]);
        } else {
            eff_sample_size[i] = (n_samples - n_outliers) * model.component_weights[i];
        }
    }

/*
    double sigma_m = 0.0;
    double sigma_f = 0.0;
    int nm = 0; int nf = 0;
    if (*window_coords == 'X') {
        for (j = 0; j < n_samples; j++) {
            if (is_outlier[j]) continue;
            if (sample_sexes[j] == 'M') {
                sigma_m += (cov_copy[j] + 0.5) * (cov_copy[j] + 0.5);
                nm++;
            } else {
                sigma_f += cov_copy[j] * cov_copy[j];
                nf++;
            }
        }
        sigma_m = sqrt(sigma_m / nm); sigma_f = sqrt(sigma_f / nf);
    }
*/
    printf("%s\t%d\t%.3f\t%.3f\t%d\t%.4g\t%.4g\t%.4g",
            window_coords, model.max_cn,
            window_gc, window_mappability,
            model.homozyg_del_is_point, model.lambda,
            model.mu_dip, model.sigma_dip);
    for (i = 0; i <= model.max_cn; i++)
        printf("\t%.1f", eff_sample_size[i]);
    for (i = model.max_cn+1; i <= MAX_CN; i++)
        printf("\t%.1f", 0.0);
    //if (*window_coords == 'X')
    //    printf("\t%.4g\t%.4g", sigma_f, sigma_m);
    printf("\n");
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s coverage.filelist.txt windows.bed [OPTIONS] >models.out\n\n", argv[0]);
        fputs("Given normalized coverage tracks (BED files) for a reference panel of samples,\n", stderr);
        fputs("fits a mixture model to the coverage distribution at each genomic window\n", stderr);
        fputs("in the tracks. CLAMMS can then use these models to call CNVs\n", stderr);
        fputs("for individual samples (even ones not in the reference panel).\n\n", stderr);
        fputs("We recommend a reference panel of at least 50 samples.\n", stderr);
        fputs("The maximum # of samples this program can handle is equal to\n", stderr);
        fputs("[maximum # of open files on the system] - 2.\n\n", stderr);
        fputs("Relative paths to each input coverage file should be listed\n", stderr);
        fputs("in coverage.filelist.txt on separate lines. An optional second column\n", stderr);
        fputs("on each line specifies the sex of the sample with a character 'M' or 'F'.\n", stderr);
        fputs("You must specify a sex for every sample or omit the second column entirely.\n", stderr);
        fputs("If you omit sample sexes, windows on sex chromosomes will be filtered out.\n\n", stderr);
        fputs("  --min_dip_cov      Windows where the median normalized coverage\n", stderr);
        fputs("                     for diploid samples is less than this fraction of the\n", stderr);
        fputs("                     median normalized coverage for windows with similar GC content\n", stderr);
        fputs("                     will be flagged and ignored during CNV calling.\n", stderr);
        fputs("                     Default = 0.1\n", stderr);
        fputs("  --max_sigma        Windows where the gaussian components in the mixture model\n", stderr);
        fputs("                     have a standard deviation parameter greater than this\n", stderr);
        fputs("                     will be flagged and ignored during CNV calling.\n", stderr);
        fputs("                     Default = 0.25\n", stderr);
        fputs("  --min_gc           If you used a non-default --min_gc for 'normalize_coverage'\n", stderr);
        fputs("                     then you must use it again here.\n", stderr);
        fputs("  --max_gc           If you used a non-default --max_gc for 'normalize_coverage'\n", stderr);
        fputs("                     then you must use it again here.\n", stderr);
        fputs("  --min_mappability  If you used a non-default --min_mappability for 'normalize_coverage'\n", stderr);
        fputs("                     then you must use it again here.\n\n", stderr);
        return 1;
    }

    FILE *filelist = open_file(argv[1]);
    FILE *windows  = open_file(argv[2]);

    Options options = parse_args(argc, argv, 3);

    int i;
    char *line = NULL;
    char *pos;
    size_t line_len;
    ssize_t bytes_read;
    int n_samples = count_lines_in_file(filelist);

    // open all input files
    i = 0;
    FILE *coverage_inputs[n_samples];
    char *sample_sexes = (char *) malloc(n_samples * sizeof(char));
    while ((bytes_read = getline(&line, &line_len, filelist)) != -1) {
        pos = strchr(line, '\t');
        if (pos == NULL) {
            sample_sexes[i] = 0;
            line[strlen(line)-1] = '\0';
        } else {
            pos++;
            sample_sexes[i] = *pos;
            line[pos-line-1] = '\0';
        }

        coverage_inputs[i] = open_file(line);
        i++;
    }

    int sexes_available = 0;
    int n_male = 0;
    int n_female = 0;
    for (i = 0; i < n_samples; i++) {
             if (sample_sexes[i] == 'M') { sexes_available++; n_male++; }
        else if (sample_sexes[i] == 'F') { sexes_available++; n_female++; }
    }
    if (sexes_available > 0 && sexes_available != n_samples) {
        fputs(          "ERROR: sexes must either be specified for all samples\n", stderr);
        fputs(          "       or be left unspecified for all samples.\n", stderr);
        fprintf(stderr, "       input has sex annotations for %d out of %d samples.",
                        sexes_available, n_samples);
        return 1;
    } else if (!sexes_available) {
        free(sample_sexes);
        sample_sexes = NULL;
    }

    int n_windows_total     = 0;
    int n_autosomal_windows = 0;
    int n_windows_processed = 0;
    int n_filtered_windows  = 0;
    int n_outlier_points    = 0;

    char window_coords[256];
    for (i = 0; i < 256; i++) // valgrind complains if we don't initialize this
        window_coords[i] = '\0';
    int window_coord_strlen;

       int  *expected_cn       = (int     *) malloc(n_samples * sizeof(int));
    double  *cov               = (double  *) malloc(n_samples * sizeof(double));
    double  *cov_copy          = (double  *) malloc(n_samples * sizeof(double));
       int  *is_outlier        = (int     *) malloc(n_samples * sizeof(int));
    double  *component_weights = (double  *) malloc(MAX_CN * sizeof(double));
    double  *male_comp_weights = (double  *) malloc(MAX_CN * sizeof(double));
    double **likelihoods       = (double **) malloc(MAX_CN * sizeof(double *));
    double **soft_assignments  = (double **) malloc(MAX_CN * sizeof(double *));
    for (i = 0; i <= MAX_CN; i++) {
        likelihoods[i]      = (double *) malloc(n_samples * sizeof(double));
        soft_assignments[i] = (double *) malloc(n_samples * sizeof(double));
    }

    while ((bytes_read = getline(&line, &line_len, coverage_inputs[0])) != -1) {
        // read first sample's coverage
        pos = strchr(line, '\t');
        pos = strchr(pos+1, '\t');
        pos = strchr(pos+1, '\t');
        window_coord_strlen = pos - line;
        strncpy(window_coords, line, window_coord_strlen);
        window_coords[window_coord_strlen] = '\0';
        cov[0] = strtof(pos+1, &pos);

        // read subsequent samples' coverage
        for (i = 1; i < n_samples; i++) {
            bytes_read = getline(&line, &line_len, coverage_inputs[i]);
            if (bytes_read == -1) {
                fprintf(stderr, "ERROR: all input coverage files must have values for the exact same windows\n");
                return 1;
            }
            pos = strchr(line, '\t');
            pos = strchr(pos+1, '\t');
            pos = strchr(pos+1, '\t');
            cov[i] = strtof(pos+1, NULL);
        }

        // read window gc content and mappability
        double window_gc;
        double window_mappability;
           int max_mm_cn;
        bytes_read = getline(&line, &line_len, windows);
        if (strncmp(window_coords, line, window_coord_strlen) != 0) {
            fprintf(stderr, "ERROR: input file %s does not appear to be the same file ", argv[2]);
            fprintf(stderr, "that the 'normalize_coverage' program was run with.");
            return 1;
        } else {
            pos = strchr(line, '\t');
            pos = strchr(pos+1, '\t');
            pos = strchr(pos+1, '\t');
            pos = strchr(pos+1, '\t');
            pos = strchr(pos+1, '\t');
            window_gc = strtod(pos+1, &pos);
            window_mappability = strtod(pos+1, &pos);
            sscanf(++pos, "%d", &max_mm_cn);
        }
/*
        if (*window_coords != 'X') {
            filter_window(window_coords, window_gc, window_mappability, 0.0, 0.0);
            continue;
        }
*/
        if ( window_gc < options.min_gc ||
             window_gc > options.max_gc ||
             window_mappability < options.min_mappability ||
             max_mm_cn < 0 ||
            (!sample_sexes && (window_coords[0] == 'X' || window_coords[0] == 'Y'))) {
            filter_window(window_coords, window_gc, window_mappability, 0.0, 0.0);
        } else {
            process_window(window_coords, max_mm_cn,
                           window_gc, window_mappability,
                           n_samples, n_male, n_female, 
                           sample_sexes, expected_cn, &options,
                           cov, cov_copy, is_outlier,
                           component_weights, male_comp_weights,
                           likelihoods, soft_assignments,
                           &n_filtered_windows, &n_outlier_points);
            n_windows_processed++;
        }

        n_windows_total++;
        if (!(window_coords[0] == 'X' || window_coords[0] == 'Y'))
            n_autosomal_windows++;
    }

    if (sample_sexes) {
        fprintf(stderr, "%d windows in total\n", n_windows_total);
        fprintf(stderr, "%d windows pass gc/mappability filter (%.1f%%)\n",
            n_windows_processed,
            100.0 * ((double) n_windows_processed) / n_windows_total);
        fprintf(stderr, "%d windows pass model-fitting filters (%.1f%%)\n",
                n_windows_processed - n_filtered_windows,
                100.0 * ((double) (n_windows_processed - n_filtered_windows)) / n_windows_total);
        fprintf(stderr, "%.2e outlier data points (%.3g per sample per window)\n",
                (double) n_outlier_points,
                ((double) n_outlier_points) / (n_samples * n_windows_processed));
    } else {
        fprintf(stderr, "%d windows in total (excluding sex chromosomes)\n",
                n_autosomal_windows);
        fprintf(stderr, "%d windows pass gc/mappability filter (%.1f%%)\n",
            n_windows_processed,
            100.0 * ((double) n_windows_processed) / n_autosomal_windows);
        fprintf(stderr, "%d windows pass model-fitting filters (%.1f%%)\n",
                n_windows_processed - n_filtered_windows,
                100.0 * ((double) (n_windows_processed - n_filtered_windows)) / n_autosomal_windows);
        fprintf(stderr, "%.2e outlier data points (%.3g per sample per window)\n",
                (double) n_outlier_points,
                ((double) n_outlier_points) / (n_samples * n_windows_processed));
    }

    // cleanup

    if (line) free(line);
    if (sample_sexes) free(sample_sexes);
    free(expected_cn);
    free(cov);
    free(cov_copy);
    free(is_outlier);
    free(component_weights);
    free(male_comp_weights);

    for (i = 0; i <= MAX_CN; i++) {
        free(likelihoods[i]);
        free(soft_assignments[i]);
    }

    free(likelihoods);
    free(soft_assignments);

    for (i = 0; i < n_samples; i++)
        fclose(coverage_inputs[i]);
    fclose(windows);
    fclose(filelist);

    return 0;
}

