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

#ifndef AUXILLARY_UTILITIES_HPP
#define AUXILLARY_UTILITIES_HPP

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <queue>
#include <vector>
#include <string>
#include <utility>

#include "correlation-matrix.hpp"
#include "geneData.hpp"
#include "graph.hpp"
#include "upper-diagonal-square-matrix.t.hpp"

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::pair;
using std::size_t;
using std::string;
using std::queue;
using std::vector;


typedef unsigned char u8;
typedef unsigned int  u32;

typedef double f64;

typedef const unsigned char cu8;
typedef const unsigned int  cu32;
typedef const std::size_t csize_t;
typedef const double cf64;

////////////////////////////////////////////////////////////////////////
//STRUCTS///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  Describe relevant configuration information for a run of TF-cluster
 **********************************************************************/
struct config{  
  s8 *tflist;
  s8 *exprData;
  s8 *corrMethod;
  
  double threeSigma, twoSigma, oneSigma;
  u8 threeSigmaAdj, twoSigmaAdj, oneSigmaAdj;
  
  u8 keepTopN;
};


////////////////////////////////////////////////////////////////////////
//PUBLIC/ FUNCTION DECLARATIONS/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



/*******************************************************************//**
 *  Make a mostly usable graph for triple-link clustering given
 * triple-link's related constraints.
 *
 * @param[in,out] protoGraph Holds gene names and an Upper Diagonal
                             matrix.  During this function, UDMatrix is
                             free'd.
 * @param[in] threeSigma Value which each row much and an edge greater
                         or equal to.
 * @param[in] oneSigma Value which any edge that is less than is not
                       used for triple-link.
 * @param[in] manNumEdges Upper bound on the number of edges a vertex
                          can have.
 **********************************************************************/
//TODO: update doc
UpperDiagonalSquareMatrix<u8>* constructCoincidenceMatrix(const CMF &protoGraph, 
                                              struct config &settings);
                                          

//TODO: add doc
graph<geneData, u8>* constructGraph(UpperDiagonalSquareMatrix<u8> *SCCM,
                        const CMF &protoGraph, struct config &settings);


/*******************************************************************//**
 *  Fail-proof (though slow) way of limiting the number of edges for
 * each vertex does not exceed the maximum value specified in the
 * configuration file.
 *
 * @param[in,out] corrData Graph which will have it's excess edges
                           removed.
 * @param[in] keepTopN Maximum number of edges a vertex may have, as
                       specified in the settings file.
 **********************************************************************/
//void pruneGraph(graph<geneData, f64> *corrData, cu8 keepTopN);


/*******************************************************************//**
 *  Using the quick-merge algorithm, sort an array of pairs by it's
 * first value.  Sorts high to low.
 *
 * @param[in,out] toSort Array of pairs to sort, and contains the sorted
                         result.
 * @param[in] size Number of pairs in toSort.
 **********************************************************************/
void sortDoubleSizeTPairHighToLow(pair<f64, size_t> *toSort,
                                                          csize_t size);


/*******************************************************************//**
 *  Using the quick-merge algorithm, sort an array of pairs by it's
 * first value.  Sorts low to high.
 *
 * @param[in,out] toSort Array of pairs to sort, and contains the sorted
                         result.
 * @param[in] size Number of pairs in toSort.
 **********************************************************************/
void sortDoubleSizeTPairLowToHigh(pair<f64, size_t> *toSort,
                                                          csize_t size);


/*******************************************************************//**
 *  Set all elements in an array to their respective absolute values.
 *
 * @param[in,out] array Array of values to be changed to their absolute
 *                      values.
 * @param[in] size Number of elements in array.
 **********************************************************************/
void inPlaceAbsoluteValue(f64 *array, csize_t size);


//TODO add doc
pair<u8, size_t>* countingSortHighToLow(pair<u8, size_t> *toSort, 
                                                            csize_t n);


////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif
