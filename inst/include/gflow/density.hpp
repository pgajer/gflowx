#ifndef DENSITY_H
#define DENSITY_H

#include <Rinternals.h> // Needed for SEXP type definition

double kernel_specific_bandwidth(const std::vector<double>& x, int kernel_type);

/**
* @brief Structure to hold density estimation results
*/
struct density_t {
   double density;        // Estimated density value
   double bandwidth;      // Bandwidth used (either provided or automatically selected)
   bool auto_selected;    // Flag indicating if bandwidth was automatically selected
};

density_t estimate_local_density(
   const std::vector<double>& x,
   int center_idx,
   double pilot_bandwidth,
   int kernel_type,
   bool verbose = false);

struct gdensity_t {
    int grid_size;        // 1d grid size
    double offset;        // offset of the start and end from min(x) and max(x)
    double start;         // start point of the grid
    double end;           // end point of the grid
    std::vector<double> density;  // data density estimate over a uniform grid
    double bandwidth;     // Bandwidth used (either provided or automatically selected)
    bool auto_selected;   // Flag indicating if bandwidth was automatically selected

    // Constructor can ensure the vector has the right size
    gdensity_t(int size) : grid_size(size), density(size) {}
};

gdensity_t estimate_local_density_over_grid(
    const std::vector<double>& x,
    int grid_size,
    double poffset,
    double pilot_bandwidth,
    int kernel_type,
    bool verbose);

#endif // DENSITY_H
