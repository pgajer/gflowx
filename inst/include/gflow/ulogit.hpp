#ifndef ULOGIT_HPP
#define ULOGIT_HPP

#include <R.h>
#include <Rinternals.h>

#include <Eigen/Dense>
#include <vector>
#include <string>

/**
* @struct eigen_ulogit_t
* @brief Structure containing results from Eigen-based logistic regression fit
*
* @details This structure holds the fitted model parameters and diagnostic information
* from logistic regression using the Eigen linear algebra library. It includes both
* model fit parameters and various leave-one-out cross-validation (LOOCV) error
* measures for model validation, as well as numerical stability warnings.
*
* The structure contains three types of LOOCV prediction errors:
* 1. Deviance errors: Based on the binomial deviance
*    $$-[y_i\log(\hat{p}_{(-i)}) + (1-y_i)\log(1-\hat{p}_{(-i)})]$$
*
* 2. Brier score errors: Squared difference between prediction and actual value
*    $$(\hat{p}_{(-i)} - y_i)^2$$
*
* 3. Absolute errors: Absolute difference between prediction and actual value
*    $$|y_i - \hat{p}_{(-i)}|$$
*
* where $$\hat{p}_{(-i)}$$ represents the predicted probability for observation i
* using a model fitted without that observation.
*
* The fitting process includes checks for numerical stability and perfect separation,
* with relevant warnings stored in the warnings vector.
*/
struct eigen_ulogit_t {
   /** @brief Vector of fitted probabilities (μᵢ) for each observation
    *  @details Values are bounded between tolerance and 1-tolerance */
    std::vector<double> predictions;

   /** @brief Vector of LOOCV deviance-based prediction errors
    *  @details Calculated using binomial deviance formula for each leave-one-out prediction */
    std::vector<double> loocv_deviance_errors;

   /** @brief Vector of LOOCV Brier score prediction errors
    *  @details Calculated as squared differences between true values and leave-one-out predictions */
    std::vector<double> loocv_brier_errors;

   /** @brief Vector of LOOCV absolute prediction errors
    *  @details Calculated as absolute differences between true values and leave-one-out predictions */
    //std::vector<double> loocv_abs_errors;

   /** @brief Model coefficients vector
    *  @details Contains:
    *   - beta[0]: Intercept term
    *   - beta[1]: Linear term coefficient
    *   - beta[2]: Quadratic term coefficient (if fit_quadratic=true)
    */
    Eigen::VectorXd beta;

   /** @brief Convergence status of the fitting algorithm
    *  @details True if algorithm converged within max_iterations and tolerance */
    bool converged;

   /** @brief Number of iterations used in model fitting
    *  @details Will be ≤ max_iterations if converged=true,
    *           equal to max_iterations if converged=false */
    int iterations;

   /** @brief Vector of Rf_warning messages from the fitting process
    *  @details Contains warnings about numerical issues such as:
    *   - Perfect or quasi-complete separation
    *   - Numerical stability issues in IWLS steps
    *   - Extreme fitted probabilities (0 or 1)
    *   - Additional ridge penalties applied for stability
    */
    std::vector<std::string> warnings;

    bool fit_quadratic;

public:
    std::vector<double> predict(const std::vector<double>& x) const {
        // Check for empty input
        if (x.empty()) {
            return std::vector<double>();
        }

        // Verify beta has correct dimensions based on model type
        int expected_size = fit_quadratic ? 3 : 2;
        if (beta.size() != expected_size) {
            Rf_error("Model coefficients vector has incorrect size");
        }

        int n = x.size();
        std::vector<double> predictions(n);

        for (int i = 0; i < n; ++i) {
            double eta = beta(0) + beta(1) * x[i];
            if (fit_quadratic) {
                eta += beta(2) * x[i] * x[i];
            }
            predictions[i] = 1.0 / (1.0 + std::exp(-eta));
        }

        return predictions;
    }
};

eigen_ulogit_t eigen_ulogit_fit(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    bool fit_quadratic,  //  = false
    int max_iterations,  //  = 100
    double ridge_lambda, // = 0.002
    double tolerance,    // = 1e-8,
    bool with_errors);   //  = false

struct ulogit_t {
    std::vector<double> w;        // weight values for the x_values
    int x_min_index = -1;             // index of smallest x value in dataset
    int x_max_index = -1;             // index of largest x value in dataset

    // Model evaluation
    std::vector<double> predictions;   // predicted probabilities
    std::vector<double> errors;        // Leave-One-Out Cross-Validation log-loss errors

    // debugging info
    int iteration_count = 0;
    bool converged = false;
};

ulogit_t ulogit(const double* x,
                const double* y,
                const std::vector<double>& w,
                int max_iterations,   // = 100
                double ridge_lambda,  // = 0.002
                double max_beta,      // = 100.0
                double tolerance,     // = 1e-8
                bool verbose);        // = false

std::vector<double> ulogit_predict(
    const double* x,
    const double* y,
    const std::vector<double>& w,
    int max_iterations,    //  = 100
    double ridge_lambda,   //  = 0.002
    double max_beta,       //  = 100
    double tolerance);     //  = 1e-8

#endif // ULOGIT_HPP
