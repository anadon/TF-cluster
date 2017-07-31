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

#ifndef CORRELATION_MATRIX_HPP
#define CORRELATION_MATRIX_HPP

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <utility>
#include <unordered_map>

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::pair;
using std::size_t;
using std::string;
using std::unordered_map;
using std::vector;

////////////////////////////////////////////////////////////////////////
//PULIC STRUCTS DECLARATIONS////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

struct geneExprRecord; //gene expression record

struct geneCorrRecord;

struct expressionData;

struct UDCorrelationMatrix;

////////////////////////////////////////////////////////////////////////
//TYPE DEFINITIONS//////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

typedef char s8;

typedef double f64;

typedef const s8 cs8;

typedef const double cf64;

////////////////////////////////////////////////////////////////////////
//PULIC STRUCTS DEFINITIONS/////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

struct geneExprRecord{ //gene expression record
  size_t nameIndex;
  f64 *exprData;
};

typedef struct geneExprRecord GER;


struct geneCorrRecord{
  size_t nameIndex;
  f64 corrCoeff;
};

typedef struct geneCorrRecord GCR;


/*******************************************************************//**
 *  Holds a table representing values to be used later to construct a
 * correlation matrix.
 **********************************************************************/
struct expressionData{
  vector< size_t > TFCorrData;//indexes into geneCorrData
  vector< GER > geneCorrData;
  size_t corrVecLeng;
  vector<string> genes, TFs;
};

typedef struct expressionData ED;


/*******************************************************************//**
 *  TODO -- update this documentation
 **********************************************************************/
struct correlationMatrixFull{
  double **fullMatrix;
  vector<string> TFLabels;
  vector<string> GeneLabels;
  inline size_t numRows() const {return TFLabels.size();}
  inline size_t numCols() const {return GeneLabels.size();}
  //inline const size_t numRows(){return TFLabels.size();}
  //inline const size_t numCols(){return GeneLabels.size();}
};

typedef struct correlationMatrixFull CMF;


/*******************************************************************//**
 *  TODO -- update this documentation
 **********************************************************************/
struct correlationMatrixDistilled{
  GCR **sparseMatrix;
  size_t numRows;
  size_t numCols;
  vector<string> TFLabels;
  vector<string> GeneLabels;
};

typedef struct correlationMatrixDistilled CMD;

////////////////////////////////////////////////////////////////////////
//PUBLIC FUNCTION DECLARATIONS//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  Read in and load rows of expression data.  Each row must be of the
 * following form:
 *    <NAME> <FLOAT 1> ... <FLOAT N>
 *
 * @param[in] expressionFile File path to correlation table
 **********************************************************************/
ED loadExpressionData(const string &expressionFile);


/*******************************************************************//**
 *  Correlation data in UDMatrix is represented in a value in [-1, 1],
 * but often it is useful to use standard deviations.  This function
 * converts all values in the UDMatrix array from a coefficient to a
 * corosponding sigma value.
 *
 * @param[in,out] input Assumed UDMatrix of correlation coefficients and
                        convert to sigma values.
 * \pre This function has not been called on input.
 **********************************************************************/
//void convertCoeffToSigmaValue(CM &input);


/*******************************************************************//**
 *  Heavily optimized method to calculate a correlation matrix.
 * Correlation is calculated using Pearson's correlation coefficient.
 *
 * @param[in] input Read in correlation data which can be clculated on.
 * @param[in] sumsOfSquares Some pre-calculated sum of squares for each
 *                          input row.
 **********************************************************************/
CMF calculatePearsonCorrelationMatrix(const ED &input,
                                                  cf64 *sumsOfSquares);


/*******************************************************************//**
 *  Heavily optimized method to calculate a correlation matrix.
 * Correlation is calculated using Spearman's correlation coefficient.
 *
 * @param[in] input Read in correlation data which can be clculated on.
 * @param[in] sumsOfSquares Some pre-calculated sum of squares for each
 *                          input row.
 **********************************************************************/
CMF calculateRankCorrelationMatrix(const ED &input);


/*******************************************************************//**
 *  Heavily optimized method to calculate a correlation matrix.
 * Correlation is calculated using Kendall's tau correlation
 * coefficient.
 *
 * @param[in] input Read in correlation data which can be clculated on.
 * @param[in] sumsOfSquares Some pre-calculated sum of squares for each
 *                          input row.
 **********************************************************************/
CMF calculateKendallsTauCorrelationCorrelationMatrix(const ED &input);


/*******************************************************************//**
 *  A wrapper function to handle loading of correlation data from a
 * file, handle sinple bad data, and in an efficient way generate a
 * correlation matrix with coefficient values (not sigma).
 *
 * @param[in] expressionFile File path to file with correlation data.
 **********************************************************************/
CMF generateMatrixFromFile(cs8 *expressionFile, cs8 *geneList,
                                                  cs8 *correlationType);


/*******************************************************************//**
 *  TODO: doc
 **********************************************************************/
bool isCorrelationAvailable(const char *corrType);

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif
