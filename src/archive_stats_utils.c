#include <R.h>
#include <stdlib.h>

static int compare_double(const void *lhs, const void *rhs)
{
    const double a = *(const double *) lhs;
    const double b = *(const double *) rhs;
    return (a > b) - (a < b);
}

double median(double *values, int n)
{
    if (n <= 0) {
        Rf_error("median cannot be computed for an empty array");
    }
    qsort(values, (size_t) n, sizeof(double), compare_double);
    if (n % 2 == 0) {
        return (values[n / 2 - 1] + values[n / 2]) / 2.0;
    }
    return values[n / 2];
}
