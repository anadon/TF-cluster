/*******************************************************************//**
         FILE:  auxillaryUtilities.hpp

  DESCRIPTION:  Miscelaneous functions used in TF-cluser

         BUGS:  ---
        NOTES:  ---
       AUTHOR:  Josh Marshall <jrmarsha@mtu.edu>
      COMPANY:  Michigan technological University
      VERSION:  See git log
      CREATED:  See git log
     REVISION:  See git log
     LISCENSE:  GPLv3
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
  string expressionFile;
  string geneList;
  u8 keepTopN;
  size_t kickSize;
  f64 threeSigma;
  f64 twoSigma;
  f64 oneSigma;
};


/*******************************************************************//**
 *  Struct to ease parameter passing from addTopEdges() to
 * addTopEdgesHelper()
 **********************************************************************/
struct addTopEdgesHelperStruct{
  struct correlationMatrix *protoGraph;
  u8 keepTopN;
  size_t startIndex;
  size_t endIndex;
};

////////////////////////////////////////////////////////////////////////
//PUBLIC/ FUNCTION DECLARATIONS/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  Verify passed arguments are valid
 **********************************************************************/
int verifyInput(int argc, char **argv, struct config settings);


/*******************************************************************//**
 *  Load configuration file data into a struct config .
 **********************************************************************/
struct config loadConfig(const char *configFilePath);

//struct loadFromFileReturn loadFromFile(const string geneListFile);


/*******************************************************************//**
 *  Sort edges high to low in a a number of vertexes.
 *
 * @param[in,out] toPrune array of vertexes to sort their edges.
 * @param[in] size Number of vertexes in toPrune.
 **********************************************************************/
//void quickMergeEdges(vertex<geneData, f64> *toPrune, csize_t size);


/*******************************************************************//**
 *  Print clusters of (presumably) genes.
 *
 * @param[in] clusters Hold name indexes.
 * @param[in] names Hold the names to be printed.
 **********************************************************************/
void printClusters(queue< queue<size_t> > clusters,
                                          const vector<string> &names);


/*******************************************************************//**
 * \deprecated {No longer meaningful for most uses.}
 * Print the edges and the vertexes each is connected to in a graph.
 *
 * @param[in] corrData Graph containing edges, and vertexes to be
 * printed.
 **********************************************************************/
void printEdges(graph<geneData, f64> *corrData);


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
graph<geneData, f64>* constructGraph(const CMF &protoGraph, 
                                        cf64 oneSigma, cu8 maxNumEdges);


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


////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif