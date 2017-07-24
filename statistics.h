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

#ifndef STATISTICS_HPP
#define STATISTICS_HPP


////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <stdlib.h>


#ifdef __cplusplus
extern "C" {
#endif

////////////////////////////////////////////////////////////////////////
//TYPE DEFINITIONS//////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


typedef unsigned char u8;
typedef unsigned int  u32;

typedef double f64;

typedef const unsigned int  cu32;
typedef const size_t csize_t;
typedef const double cf64;

////////////////////////////////////////////////////////////////////////
//PUBLIC FUNCTION DECLARATIONS//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  Calculate the centered standard deviation of an array
 *
 * @param[in] array Array with a mean of 0.  Must have more than 1
                    element.
 * @param[in] size Number of elements in array.
 **********************************************************************/
f64  calculateCenteredStandardDeviation(cf64 *array, csize_t size);


/*******************************************************************//**
 *  Return a new array from the variable array with each element from
 * the inputted array having the mean of the passed in array from it, so
 * the returned array has a mean of 0.
 *
 * @param[in] array Array of values.
 * @param[in] size Number of elements in array.
 **********************************************************************/
f64* centerMean(cf64 *array, csize_t size);


/*******************************************************************//**
 *  Calculate the covariance between 2 equal length arrays.
 *
 * @param[in] left An array of values.
 * @param[in] leftSum Sum of values in left.
 * @param[in] right An array of values.
 * @param[in] rightSum Sum of values in right.
 * @param[in] number of elements in left and right.
 **********************************************************************/
f64  covarianceOne(cf64 *left,  cf64 leftSum, cf64 *right, cf64 rightSum,
                                                  csize_t arrayLength);



f64  covarianceTwo(cf64 *left, cf64 *right, csize_t arrayLength);


/*******************************************************************//**
 *  Get the mean of an array.
 *
 * @param[in] array Array of values.
 * @param[in] size number of elements in array.
 **********************************************************************/
f64  getMean(cf64 *array, csize_t size);


/*******************************************************************//**
 *  Calculate the standard deviation of array.
 *
 * @param[in] array Array of values.
 * @param[in] size Number of elements in array.
 **********************************************************************/
f64  getStandardDeviation(cf64 *array, csize_t size);


/*******************************************************************//**
 *  Calculate sum of all values in array.
 *
 * @param[in] array Array of values.
 * @param[in] size Number of elements in array.
 **********************************************************************/
f64  getSum(cf64 *array, csize_t size);


/*******************************************************************//**
 *  Calculate the sum of each value in array squared.
 *
 * @param[in] array Array of values.
 * @param[in] size Number of elements in array.
 **********************************************************************/
f64  getSumOfSquares(cf64 *array, csize_t size);


/*******************************************************************//**
 *  Calculate sum of multiplied same index values in 2 arrays.
 *
 * @param[in] left Array of values.
 * @param[in] right Array of values.
 * @param[in] size Number of elements in each array.
 **********************************************************************/
f64  getSumOfMultipliedArrays(cf64 *left, cf64 *right, csize_t size);


/*******************************************************************//**
 *  Calculate correlation of two arrays, each of which having a mean of
 * 0.  This is using Spearman's Correlation Coefficient formula.
 *
 * @param[in] left Array of values.
 * @param[in] right Array of values.
 * @param[in] size number of elements in each array.
 **********************************************************************/
f64  getCenteredCorrelation(cf64 *left, cf64 *right, csize_t size);


/*******************************************************************//**
 *  Calculate correlation of two arrays from the sum of squares of each
 * array and the cross sum over those two arrays.  This is using
 * Spearman's Correlation Coefficient formula.
 *
 * @param[in] aSumOfSquares Value of sum of squares of an array
 * @param[in] bSumOfSquares Value of sum of squares of an array
 * @param[in] abCrossSum Value of the cross sum of two arrays
 **********************************************************************/
f64 getCenteredCorrelationBasic(cf64 aSumOfSquares, cf64 bSumOfSquares,
                                                      cf64 abCrossSum);


/*******************************************************************//**
 *  Change an array so that it has a mean of 0.
 *
 * @param[in,out] array Array to center about 0.
 * @param[in] size Number of elements in array.
 **********************************************************************/
void inplaceCenterMean(f64 *array, csize_t size);


//TODO: add doc

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

#endif
