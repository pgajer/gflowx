#ifndef ONED_LINEAR_MODELS_H_
#define ONED_LINEAR_MODELS_H_

#include <algorithm>  // for std::find, std::clamp
#include <vector>
#include <utility>   // for std::pair
#include "error_utils.h"
#include "cpp_utils.hpp"

struct lm_fit_t {
    double predicted_value;
    double r_squared;
    double mae;      // Mean Absolute Error
    double rmse;     // Root Mean Square Error
    double aic;      // Akaike Information Criterion
};

struct lm_mae_t {
    double predicted_value;
    double mae;      // Mean Absolute Error
};

struct lm_dmudy_t {
    double prediction;
    double derivative;
};

struct lm_loocv_t {
    // Model components
    double slope;
    double y_wmean;    // weighted mean of y
    double x_wmean;    // weighted mean of x (after centering around ref_point)
    double x_ref;      // reference point value
    int ref_index;     // index of the reference point in x and y vectors
    std::vector<int> vertex_indices;    // indices of vertices in the training set
    std::vector<double> x_values;       // original x values for these vertices
    std::vector<double> w_values;       // weight values for the x_values

    // Model evaluation
    double predicted_value;     // prediction at reference point
    double loocv;              // Leave-One-Out Cross-Validation MSE
    double loocv_at_ref_vertex; // LOOCV squared error at reference vertex

    // Method to compute prediction at any vertex from the training set
    double predict(int vertex_idx) const {
        #if 0
        Rprintf("predict() called with vertex_idx: %d\n", vertex_idx);
        Rprintf("this->vertex_indices contents: ");
        for(const auto& idx : this->vertex_indices) {
            Rprintf("%d, ", idx);
        }
        Rprintf("\n");
        #endif
        // Find position of vertex_idx in vertex_indices
        auto it = std::find(vertex_indices.begin(), vertex_indices.end(), vertex_idx);
        if (it == vertex_indices.end()) {
            Rprintf("In predict(): vertex_idx: %d\n", vertex_idx);
            print_vect(vertex_indices,"vertex_indices");
            REPORT_ERROR("Vertex index %d not found in training set", vertex_idx);
        }
        int pos = std::distance(vertex_indices.begin(), it);
        return predict_at_x(x_values[pos]);
    }

private:
    // Internal helper to compute prediction at an x value
    double predict_at_x(double x) const {
        double x_adj = x - x_ref;
        return y_wmean + slope * (x_adj - x_wmean);
    }
};

#endif // ONED_LINEAR_MODELS_H_
