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

#include <graph.tpp>
#include <kendall-correlation-matrix.hpp>
#include <pearson-correlation-matrix.hpp>
#include <queue>
#include <spearman-correlation-matrix.hpp>
#include <string>
#include <upper-diagonal-square-matrix.tpp>
#include <utility>
#include <tf-cluster-files.hpp>
#include <vector>

#include "geneData.hpp"

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
  s8 *tflistFile;//TODO: cascade name change
  s8 *tflistFileFormat;//TODO: add runtime option in argp
  s8 *exprDataFile;//TODO: cascade name change
  s8 *exprDataFileFormat;//TODO: add runtime option in argp

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
UpperDiagonalSquareMatrix<u8>* constructCoincidenceMatrix(const double **fullMatrix,
           const ED fileData, struct config settings);


//TODO: add doc
graph<geneData, u8>* constructGraph(UpperDiagonalSquareMatrix<u8> *SCCM, const ED fileData,
                                                      struct config &settings);

//TODO: add doc
double** generateMatrixFromFile(ED fileData, unordered_map<string, char> TFCheck,
                                                  cs8 *correlationType);



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
