#ifndef PREDICTIVE_ERRORS_HPP_
#define PREDICTIVE_ERRORS_HPP_

#include <vector>      // For std::vector usage throughout

struct bb_cri_t {
    std::vector<double> bb_Ey; // central location of the Bayesian bootstrap estimates
    std::vector<double> cri_L; // credible intervals lower limit
    std::vector<double> cri_U; // credible intervals upper limit

};

bb_cri_t bb_cri(const std::vector<std::vector<double>>& bb_Ey,
                bool use_median,
                double p);

std::vector<double> compute_bbwasserstein_errors(const std::vector<std::vector<double>>& bb_Ey,
                                                 const std::vector<double>& y);

std::vector<double> compute_bbwasserstein_errors(const std::vector<std::vector<double>>& bb_Ey,
                                                 const std::vector<double>& y,
                                                 const std::vector<std::vector<double>>& bb_y);

double compute_bbwasserstein_error(const std::vector<std::vector<double>>& bb_Ey,
                                   const std::vector<double>& y);

std::pair<double, double> compute_bbcov_error(const std::vector<std::vector<double>>& bb_Ey,
                                              const std::vector<double>& y,
                                              bool use_median);

#if 0
// in-class member initializers (C++11 and later)
struct predictive_error_t {
    double bb_cov_errorA = -1;
    double bb_cov_errorB = -1;
    double bb_wasserstein_error = -1;
    double cv_error = -1;
    double true_error = -1;
    std::vector<double> Ey;
    std::pair<std::vector<double>, std::vector<double>> Ey_cri;
};
#endif

enum class error_methods_t {
    ALL,
    BB_COV_ERRORS,
    BB_WASSERSTEIN_ERROR,
    CV_ERROR
};

#endif // PREDICTIVE_ERRORS_HPP_
