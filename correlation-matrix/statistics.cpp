
#include <cerrno>
#include <cstdlib>
#include <math.h>
#include <pthread.h>


#include "statistics.hpp"


//WARNING: This function assumes a mean of 0
f64 calculateCenteredStandardDeviation(cf64 *weights, csize_t &numEntries){
  return sqrt(getSumOfSquares(weights, numEntries) / numEntries);
}


f64 getStandardDeviation(cf64 *array, csize_t &size){
  f64 tr = 0;
  cf64 mean = getMean(array, size);
  
  for(size_t i = 0; i < size; i++){
    cf64 tmp = array[i] - mean;
    tr += (tmp * tmp);
  }
  tr = sqrt(tr / size);
  
  return tr;
}


f64 calculateSigma(cf64 *array, csize_t size){
  f64 meanOfSquareDifferences;

  cf64 mean = getSum(array, size) / size;

  meanOfSquareDifferences = 0;
  for(size_t i = 0; i < size; i++){
    f64 difference = array[i] - mean;
    meanOfSquareDifferences += (difference * difference);
  }
  meanOfSquareDifferences /= size;

  return sqrt(meanOfSquareDifferences);
}


f64 getMean(cf64 *array, csize_t &size){
  return getSum(array, size) / size;
}


f64 getSum(cf64 *array, csize_t &size){
  f64 tr = 0;
  for(size_t i = 0; i < size; i++) tr += array[i];
  return tr;
}


f64 getSumOfSquares(cf64 *array, csize_t &size){
  return getSumOfMultipliedArrays(array, array, size);
}


f64 getSumOfMultipliedArrays(cf64 *left, cf64 *right, csize_t &size){
  f64 tr = 0;
  
  for(size_t i = 0; i < size; i++)
    tr += (left[i] * right[i]);
  
  return tr;
}


//WARNING: This functions assumes left and right have a mean of 0
f64 getCenteredCorrelation(cf64 *left, cf64 *right, csize_t &size){
  
  cf64 aSumOfSquares = getSumOfSquares(left, size);
  cf64 bSumOfSquares = getSumOfSquares(right, size);
  cf64 abCrossSum = getSumOfMultipliedArrays(left, right, size);
  
  return abCrossSum / sqrt(aSumOfSquares * bSumOfSquares);
}


f64* centerMean(cf64 *array, csize_t &size){
  f64 *tr;
  
  cf64 mean = getMean(array, size);
  tr = (f64*) malloc(sizeof(*array) * size);
  
  for(size_t i = 0; i < size; i++)
    tr[i] = array[i] - mean;
  
  return tr;
}


void inplaceCenterMean(f64 *array, csize_t &size){
  cf64 mean = getMean(array, size);
  
  for(size_t i = 0; i < size; i++)
    array[i] -= mean;
}


f64 covariance(cf64 *left, cf64 leftSum, cf64 *right, cf64 rightSum, 
                                                  csize_t arrayLength){
  return (getSumOfMultipliedArrays(left, right, arrayLength) / arrayLength) / ((leftSum / arrayLength)*(rightSum / arrayLength));
}