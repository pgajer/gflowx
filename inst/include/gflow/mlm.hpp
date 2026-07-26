#ifndef MLM_HPP
#define MLM_HPP

#include "error_utils.h"  // REPORT_WARNING etc.

#include <cstddef>             // For std::size_t
#include <vector>              // For std::vector
#include <algorithm>           // For std::min
#include <map>                 // For std::map
#include <string>              // For std::string
#include <fstream>             // For std::ofstream

#include <Eigen/Dense>         // For Eigen::VectorXd, Eigen::MatrixXd

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>       // For Rprintf

using std::size_t;

/**
 * @brief Enhanced linear model structure for spectral embedding
 *
 * This structure holds the results of fitting a weighted linear model
 * in the spectral embedding space, including:
 * - Model coefficients
 * - Predictions for vertices in the model's domain
 * - Leave-one-out cross-validation errors
 * - Vertex indices and weights
 * - Summary statistics for model evaluation
 */
struct lm_t {
    // Vertex information
    std::vector<size_t> vertices;         ///< Original vertex indices in the model domain
    std::vector<double> weights;          ///< Kernel weights for each vertex

    // Model results
    std::vector<double> predictions;      ///< Predicted values for each vertex
    std::vector<double> errors;           ///< LOOCV squared errors for each vertex
    double mean_error;                    ///< Mean LOOCV error across all vertices

    // Model parameters
    Eigen::VectorXd coefficients;         ///< Fitted coefficients in the spectral embedding
    double intercept;                     ///< Intercept term

    /**
     * @brief Default constructor
     */
    lm_t() : mean_error(0.0), intercept(0.0) {}

    /**
     * @brief Make predictions for new points in the spectral embedding
     *
     * @param X_new Matrix where each row is a point in the spectral embedding
     * @return Vector of predictions
     */
    // Make predictions on new data
    std::vector<double> predict(const Eigen::MatrixXd& X_new) const {
        std::vector<double> new_predictions(X_new.rows());
        for (int i = 0; i < X_new.rows(); ++i) {
            new_predictions[i] = intercept + X_new.row(i).dot(coefficients);
        }
        return new_predictions;
    }

    // Debug method to print model details
    void print(size_t max_num_coefs = 5) const {
        Rprintf("Linear Model Summary:\n");
        Rprintf("  Intercept: %.6f\n", intercept);
        Rprintf("  Coefficients (%zu):", coefficients.size());
        for (size_t i = 0; i < std::min(max_num_coefs, (size_t)coefficients.size()); i++) {
            Rprintf(" %.6f", coefficients(i));
        }
        if (coefficients.size() > static_cast<Eigen::Index>(max_num_coefs)) {
            Rprintf(" ...");
        }
        Rprintf("\n");

        Rprintf("  Sample sizes: predictions=%zu, errors=%zu, weights=%zu\n",
                predictions.size(), errors.size(), weights.size());

        Rprintf("  Mean error: %.6f\n", mean_error);
    }
    
    // Export model data to CSV for debugging/validation
    void write_model_data_to_csv(
        const Eigen::MatrixXd& X,
        const std::vector<double>& y,
        const std::map<size_t, double>& vertex_map,
        const std::string& filename) const {

        std::ofstream file(filename);
        if (!file.is_open()) {
            REPORT_WARNING("Failed to open file for writing model data: %s\n", filename.c_str());
            return;
        }

        // Write header
        file << "vertex,distance,y,prediction,error,weight";
        for (int j = 1; j < X.cols(); j++) {
            file << ",dim" << j;
        }
        file << "\n";

        // Write data rows
        int row_idx = 0;
        for (const auto& [v_idx, distance] : vertex_map) {
            file << (v_idx + 1) << ","
                 << distance << ","
                 << y[v_idx] << ","
                 << predictions[row_idx] << ","
                 << errors[row_idx] << ","
                 << weights[row_idx];

            // Write embedding coordinates
            for (int j = 1; j < X.cols(); j++) {
                file << "," << X(row_idx, j);
            }
            file << "\n";
            row_idx++;
        }

        file.close();
        Rprintf("Model data written to: %s\n", filename.c_str());
    }
};

lm_t fit_linear_model(
    const Eigen::MatrixXd& embedding,
    const std::vector<double>& y,
    const std::map<size_t, double>& vertex_map,
    double dist_normalization_factor);

lm_t cleveland_fit_linear_model(
    const Eigen::MatrixXd& embedding,
    const std::vector<double>& y,
    const std::map<size_t, double>& vertex_map,
    double dist_normalization_factor,
    size_t n_iterations = 3,
    double robust_scale = 6.0);

#endif // MLM_HPP
