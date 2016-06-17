#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <cstdlib>


typedef unsigned char u8;
typedef unsigned int  u32;

typedef double f64;

typedef const unsigned int  cu32;
typedef const std::size_t csize_t;
typedef const double cf64;


f64  calculateCenteredStandardDeviation(cf64 *weights, csize_t &numEntries);

f64* centerMean(cf64 *array, csize_t &size);

f64  covariance(cf64 *left,  cf64 leftSum, cf64 *right, cf64 rightSum, 
                                                  csize_t arrayLength);

f64  getMean(cf64 *array, csize_t &size);

f64  getStandardDeviation(cf64 *array, csize_t &size);

f64  getSum(cf64 *array, csize_t &size);

f64  getSumOfSquares(cf64 *array, csize_t &size);

f64  getSumOfMultipliedArrays(cf64 *left, cf64 *right, csize_t &size);

f64  getCenteredCorrelation(cf64 *left, cf64 *right, csize_t &size);

void inplaceCenterMean(f64 *array, csize_t &size);

#endif