#ifndef ARRAY_UTILS_H
#define ARRAY_UTILS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl_sort.h>

#include "math_utils.h"

int array_scalar_multiply(double *data, size_t n, double value);

int array_scalar_sum(double *data, size_t n, double value);


int array_sum(double *data1, const double *data2, size_t n);

int array_substract(double *data1, const double *data2, size_t n);

int array_dotproduct(double *data1, const double *data2, size_t n, double *res);


int array_log(double *values, size_t n);

int array_log10(double *values, size_t n);

int array_log_base(double *values, size_t n, double base);


int array_accum(const double *values, size_t n, double *res);

int array_accum_range(const double *values, size_t begin, size_t end, double *res);


int order(const double *values, size_t n, int desc, size_t *indices);

int ordered(const double *values, size_t n, const int *indices, double *ordered);


#endif
