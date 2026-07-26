#ifndef PGLM_H_
#define PGLM_H_

#include "1D_linear_models.h"

#include <cstddef>
#include <vector>
#include <map>

using std::size_t;

/**
 * @brief Triplet of integers identifying a unique path model
 *
 * This struct represents a path by its start and end vertices, along with
 * a reference vertex that must lie on the path.
 */
struct itriplet_t {
    int start;     ///< Starting vertex of the path
    int end;       ///< Ending vertex of the path
    int ref_index; ///< Reference vertex on the path

    /**
     * @brief Comparison operator for ordering in std::map
     *
     * Implements lexicographical comparison of triplets.
     * Required for use as key in std::map.
     */
    bool operator<(const itriplet_t& other) const {
        if (start != other.start) return start < other.start;
        if (end != other.end) return end < other.end;
        return ref_index < other.ref_index;
    }
};

/**
 * @brief Cache structure for storing and retrieving path linear models
 *
 * This class maintains a cache of linear models computed for paths in a graph.
 * Each model is uniquely identified by a triplet of (start vertex, end vertex, reference vertex)
 * where the reference vertex must lie on the path between start and end vertices.
 *
 * The cache is designed for efficient lookup of precomputed models and includes
 * error checking to ensure all requested models exist in the cache.
 */
struct path_lm_cache_t {
private:
    /** Map storing computed linear models indexed by vertex triplets */
    std::map<itriplet_t, lm_loocv_t> model_cache;

public:
    /**
     * @brief Retrieves a linear model from the cache
     *
     * @param start Starting vertex of the path
     * @param end Ending vertex of the path
     * @param ref_index Reference vertex index on the path
     * @return lm_loocv_t The cached linear model
     * @throws Runtime error if the requested model is not found in cache
     *
     * @note Time complexity: O(log N) where N is the number of cached models
     */
    lm_loocv_t get_model(int start, int end, int ref_index) const {
        const itriplet_t key_triplet{start, end, ref_index}; // Use aggregate initialization
        auto it = model_cache.find(key_triplet);
        if (it != model_cache.end()) {
            return it->second;
        }
        REPORT_ERROR("Model not found in cache for path (%d->%d) with reference vertex %d",
                     start, end, ref_index);
        // Adding an explicit return after error to satisfy compiler
        return lm_loocv_t{}; // This line will never be reached if REPORT_ERROR terminates
    }

    /**
     * @brief Checks if a model exists in the cache
     *
     * @param start Starting vertex of the path
     * @param end Ending vertex of the path
     * @param ref_index Reference vertex index on the path
     * @return bool True if model exists in cache, false otherwise
     */
    bool has_model(int start, int end, int ref_index) const {
        return model_cache.count(itriplet_t{start, end, ref_index}) > 0;
    }

    /**
     * @brief Clears all cached models
     */
    void clear() {
        model_cache.clear();
    }

    /**
     * @brief Returns the number of models in the cache
     *
     * @return size_t Number of cached models
     */
    size_t size() const {
        return model_cache.size();
    }
};



#endif // PGLM_H_
