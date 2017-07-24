/*Copyright 2016-2017 Josh Marshall************************************/

/***********************************************************************
    This file is part of TF-Cluster.

    TF-Cluster is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TF-Cluster.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

#include "statistics.h"


#ifdef __cplusplus
extern "C" {
#endif

////////////////////////////////////////////////////////////////////////
//PUBLIC FUNCTION DEFINITIONS///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//WARNING: This function assumes a mean of 0
f64 calculateCenteredStandardDeviation(cf64 *array, csize_t size){
  f64 tr = 0;

  for(size_t i = 0; i < size; i++){
    tr += (array[i] * array[i]);
  }
  tr = sqrt(tr / (size-1));

  return tr;
}


f64 getStandardDeviation(cf64 *array, csize_t size){
  f64 tr = 0;
  cf64 mean = getMean(array, size);

  for(size_t i = 0; i < size; i++){
    cf64 tmp = array[i] - mean;
    tr += (tmp * tmp);
  }
  tr = sqrt(tr / (size-1));

  return tr;
}


f64 getMean(cf64 *array, csize_t size){
  return getSum(array, size) / size;
}


f64 getSum(cf64 *array, csize_t size){
  f64 tr = 0;
  for(size_t i = 0; i < size; i++) tr += array[i];
  return tr;
}


size_t getIntegerSum(csize_t *array, csize_t size){
  size_t tr = 0;
  for(size_t i = 0; i < size; i++) tr += array[i];
  return tr;
}


f64 getSumOfSquares(cf64 *array, csize_t size){
  return getSumOfMultipliedArrays(array, array, size);
}


f64 getSumOfMultipliedArrays(cf64 *left, cf64 *right, csize_t size){
  f64 tr = 0;

  for(size_t i = 0; i < size; i++)
    tr += (left[i] * right[i]);

  return tr;
}


//WARNING: This functions assumes left and right have a mean of 0
f64 getCenteredCorrelation(cf64 *left, cf64 *right, csize_t size){

  cf64 aSumOfSquares = getSumOfSquares(left, size);
  cf64 bSumOfSquares = getSumOfSquares(right, size);
  cf64 abCrossSum = getSumOfMultipliedArrays(left, right, size);

  return getCenteredCorrelationBasic(aSumOfSquares, bSumOfSquares,
                                                            abCrossSum);
}


//WARNING: This functions assumes left and right have a mean of 0
f64 getCenteredCorrelationBasic(cf64 aSumOfSquares, cf64 bSumOfSquares,
                                                      cf64 abCrossSum){
  return abCrossSum / sqrt(aSumOfSquares * bSumOfSquares);
}


f64* centerMean(cf64 *array, const size_t size){
  f64 *tr;

  cf64 mean = getMean(array, size);
  tr = (f64*) malloc(sizeof(*array) * size);

  for(size_t i = 0; i < size; i++)
    tr[i] = array[i] - mean;

  return tr;
}


void inplaceCenterMean(f64 *array, const size_t size){
  cf64 mean = getMean(array, size);

  for(size_t i = 0; i < size; i++)
    array[i] -= mean;
}


f64 covarianceOne(cf64 *left, cf64 leftSum, cf64 *right, cf64 rightSum,
                                                  csize_t arrayLength){
  cf64 *l = left, m = leftSum, *r = right, t = rightSum;
  csize_t s = arrayLength;
  return (getSumOfMultipliedArrays(l, r, s) / s) / ((m / s) * (t / s));
}


f64 covarianceTwo(cf64 *left, cf64 *right, csize_t arrayLength){
  cf64 *l = left, *r = right;
  csize_t s = arrayLength;
  return getSumOfMultipliedArrays(l, r, s) / (2.0 * s * s);
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif
