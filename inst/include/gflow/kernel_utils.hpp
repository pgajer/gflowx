#ifndef KERNEL_UTILS_H_
#define KERNEL_UTILS_H_

#include <vector>

/**
 * @brief Compute normalized kernel weights from a distance vector
 */
std::vector<double> get_weights(
	std::vector<double>& dists,
	double dist_normalization_factor
	);

#endif // KERNEL_UTILS_H_
