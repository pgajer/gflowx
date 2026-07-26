#ifndef SORT_UTILS_H_
#define SORT_UTILS_H_

#include <vector>
#include <cstddef>
using std::size_t;

/// Sort x (and y in parallel), and produce an “inverse” index so that
///   x_sorted[order[i]] == x[i]
void sort_by_x_keep_y_and_order(
	const std::vector<double>& x,
	const std::vector<double>& y,
	std::vector<double>& x_sorted,
	std::vector<double>& y_sorted,
	std::vector<std::size_t>& order
);

#endif // SORT_UTILS_H_
