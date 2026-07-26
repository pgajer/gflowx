#ifndef GRID_VERTEX_PATH_MODEL_MODEL_PRIORITY_QUEUE_H_
#define GRID_VERTEX_PATH_MODEL_MODEL_PRIORITY_QUEUE_H_

#include <unordered_map>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <vector>
#include <cstddef>

#include <R.h>  // For Rprintf

using std::size_t;

struct grid_vertex_path_model_t {
    size_t grid_vertex;
    std::vector<size_t> vertices;              ///< path throught grid_vertex consisting of the original graph vertices
    std::unordered_set<size_t> core_vertices;  ///< a set of vertices close by vertices; determined using diff_threshold and/or min_path_size
    std::vector<double> predictions;           ///< model predictions
    std::vector<double> errors;                ///< model's prediction error at each vertex
    double mean_error;                         ///< model's mean prediction error
};

// Custom comparator for min heap based on mean_error
struct compare_models {
    bool operator()(const grid_vertex_path_model_t& a, const grid_vertex_path_model_t& b) {
        return a.mean_error > b.mean_error; // For min heap
    }
};

class grid_vertex_path_model_priority_queue {
private:
    std::priority_queue<
        grid_vertex_path_model_t,
        std::vector<grid_vertex_path_model_t>,
        compare_models
    > pq;

    // Helper function to check if two models have intersecting core vertices
    bool has_intersection(const std::unordered_set<size_t>& set1,
                         const std::unordered_set<size_t>& set2) {
        for (const auto& vertex : set1) {
            if (set2.find(vertex) != set2.end()) {
                return true;
            }
        }
        return false;
    }

    // Helper function to find all models that intersect with given model
    std::vector<grid_vertex_path_model_t> find_intersecting_models(
        const grid_vertex_path_model_t& model) {
        std::vector<grid_vertex_path_model_t> result;
        std::vector<grid_vertex_path_model_t> temp_models;

        // Move all elements to temporary vector
        while (!pq.empty()) {
            auto current_model = pq.top();
            pq.pop();

            if (has_intersection(current_model.core_vertices, model.core_vertices)) {
                result.push_back(current_model);
            }
            temp_models.push_back(current_model);
        }

        // Restore queue
        for (const auto& temp_model : temp_models) {
            pq.push(temp_model);
        }

        return result;
    }

public:
    void custom_insert(const grid_vertex_path_model_t& model) {
        // Find all models that intersect with the new model
        auto intersecting_models = find_intersecting_models(model);

        if (intersecting_models.empty()) {
            // Case 1: No intersecting models, simply insert
            pq.push(model);
            return;
        }

        // Find model with minimum mean_error among intersecting models
        auto min_error_model = *std::min_element(
            intersecting_models.begin(),
            intersecting_models.end(),
            [](const auto& a, const auto& b) {
                return a.mean_error < b.mean_error;
            }
        );

        if (model.mean_error < min_error_model.mean_error) {
            // Case 2: New model has smaller error than all intersecting models
            // Remove all intersecting models and insert new model
            std::vector<grid_vertex_path_model_t> temp_models;
            while (!pq.empty()) {
                auto current_model = pq.top();
                pq.pop();

                bool should_keep = true;
                for (const auto& intersecting_model : intersecting_models) {
                    if (current_model.grid_vertex == intersecting_model.grid_vertex) {
                        should_keep = false;
                        break;
                    }
                }

                if (should_keep) {
                    temp_models.push_back(current_model);
                }
            }

            // Restore non-intersecting models and add new model
            for (const auto& temp_model : temp_models) {
                pq.push(temp_model);
            }
            pq.push(model);
        }
        // Case 3: There exists an intersecting model with smaller error
        // Do nothing (model is not inserted)
    }

    bool empty() const { return pq.empty(); }

    grid_vertex_path_model_t top() const { return pq.top(); }

    void pop() { pq.pop(); }

    size_t size() const { return pq.size(); }

    void print() const {
        std::vector<grid_vertex_path_model_t> temp;
        auto pq_copy = pq;  // Make a copy so we don't modify the original

        Rprintf("Queue contents (size=%zu):\n", pq.size());
        while (!pq_copy.empty()) {
            auto model = pq_copy.top();
            pq_copy.pop();

            Rprintf("Model %d: mean_error=%f, core_vertices={", 
                    model.grid_vertex, model.mean_error);

            // Print core vertices
            bool first = true;
            for (const auto& v : model.core_vertices) {
                if (!first) Rprintf(",");
                Rprintf("%d", v);
                first = false;
            }
            Rprintf("}\n");
        }
        Rprintf("\n");
    }
};

#endif // GRID_VERTEX_PATH_MODEL_MODEL_PRIORITY_QUEUE_H_
