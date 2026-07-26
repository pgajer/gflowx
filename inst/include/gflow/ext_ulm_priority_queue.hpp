#ifndef EXT_ULM_PRIORITY_QUEUE_H_
#define EXT_ULM_PRIORITY_QUEUE_H_

#include "ulm.hpp"

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>
#include <algorithm>
#include <cstddef>

#include <R.h>  // For Rprintf

using std::size_t;

// Custom comparator for min heap based on mean_error
struct compare_ext_ulm {
    bool operator()(const ext_ulm_t& a, const ext_ulm_t& b) {
        return a.mean_error > b.mean_error; // For min heap
    }
};

class ext_ulm_priority_queue {
private:
    std::priority_queue<
        ext_ulm_t,
        std::vector<ext_ulm_t>,
        compare_ext_ulm
    > pq;

    // Helper function to check if two models have intersecting vertices
    bool has_intersection(const std::vector<size_t>& vertices1,
                         const std::vector<size_t>& vertices2) {
        for (const auto& vertex : vertices1) {
            if (std::find(vertices2.begin(), vertices2.end(), vertex) != vertices2.end()) {
                return true;
            }
        }
        return false;
    }

    #if 0
    // optimized version
    bool has_intersection(const std::vector<size_t>& vertices1,
                          const std::vector<size_t>& vertices2) {
        std::unordered_set<size_t> set2(vertices2.begin(), vertices2.end());
        for (const auto& vertex : vertices1) {
            if (set2.find(vertex) != set2.end()) {
                return true;
            }
        }
        return false;
    }
    #endif

    // Helper function to find all models that intersect with given model
    std::vector<ext_ulm_t> find_intersecting_models(const ext_ulm_t& model) {
        std::vector<ext_ulm_t> result;
        std::vector<ext_ulm_t> temp_models;

        // Move all elements to temporary vector
        while (!pq.empty()) {
            auto current_model = pq.top();
            pq.pop();

            if (has_intersection(current_model.vertices, model.vertices)) {
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
    void custom_insert(const ext_ulm_t& model) {
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
            std::vector<ext_ulm_t> temp_models;
            while (!pq.empty()) {
                auto current_model = pq.top();
                pq.pop();

                bool should_keep = true;
                for (const auto& intersecting_model : intersecting_models) {
                    // Compare using vertices instead of grid_vertex
                    if (has_intersection(current_model.vertices, intersecting_model.vertices)) {
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

    ext_ulm_t top() const { return pq.top(); }

    void pop() { pq.pop(); }

    size_t size() const { return pq.size(); }

    #if 0
    void print() const {
        auto pq_copy = pq;  // Make a copy so we don't modify the original

        Rprintf("Queue contents (size=%zu):\n", pq.size());
        while (!pq_copy.empty()) {
            auto model = pq_copy.top();
            pq_copy.pop();

            Rprintf("Model mean_error=%f, vertices={", model.mean_error);

            // Print vertices
            bool first = true;
            for (const auto& v : model.vertices) {
                if (!first) Rprintf(",");
                Rprintf("%d", v);
                first = false;
            }
            Rprintf("}\n");
        }
        Rprintf("\n");
    }
    #endif

     void print() const {
        auto pq_copy = pq;  // Make a copy so we don't modify the original

        Rprintf("Queue contents (size=%zu):\n", pq.size());
        while (!pq_copy.empty()) {
            auto model = pq_copy.top();
            pq_copy.pop();

            Rprintf("Model mean_error=%.6f, vertices={", model.mean_error);

            // Print vertices
            for (size_t i = 0; i < model.vertices.size(); ++i) {
                if (i > 0) Rprintf(",");
                Rprintf("%zu", model.vertices[i]);
            }
            Rprintf("}\n");
        }
        Rprintf("\n");
    }
};

#endif // EXT_ULM_PRIORITY_QUEUE_H_
