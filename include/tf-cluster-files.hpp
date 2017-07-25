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
#include <unistd.h>
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
  union{
    size_t corrVecLeng;//exprData[][corrVecLeng]
    size_t numCols;
  };
  union{
    size_t numGenes;   //exprData[numGenes][]
    size_t numRows;
  };
  double **exprData;    //[numGenes][corrVecLeng]

  size_t TFCorrDataLength;
  size_t *TFCorrData;//indexes into geneCorrData, [TFCorrDataLength]

  vector<string> genes;
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
 *
 * TODO: update documentation
 **********************************************************************/
ED loadExpressionData(const char *expressionFile,
                                           unordered_map<string, char> TFCheck);


/*******************************************************************//**
 *  Read in a whitespace seperated list of gene names to be cross references
 * against gene names in the expression data file.
 *
 * @param[in] expressionFile File path to correlation table
 *
 * TODO: update documentation
 **********************************************************************/
unordered_map<string, char> readInTFList(const char *TFListFile);

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif
