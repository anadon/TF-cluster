#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP


////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <queue>
#include <string>
#include <vector>

#include "auxillaryUtilities.hpp"
#include "upper-diagonal-square-matrix.t.hpp"

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::queue;
using std::string;
using std::vector;

////////////////////////////////////////////////////////////////////////
//GLOBAL FUNCTION DEFINITIONS///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  Print clusters of (presumably) genes.
 *
 * @param[in] clusters Hold name indexes.
 * @param[in] names Hold the names to be printed.
 **********************************************************************/
void printClusters(queue< queue<size_t> > clusters,
                                            const vector<string> &TFs);


void printCoincidenceMatrix(const UpperDiagonalSquareMatrix<u8> matrix, 
                                cu8 maxMatch, const vector<string> TFs);



/*******************************************************************//**
 * \deprecated {No longer meaningful for most uses.}
 * Print the edges and the vertexes each is connected to in a graph.
 *
 * @param[in] corrData Graph containing edges, and vertexes to be
 * printed.
 **********************************************************************/
//void printEdges(graph<geneData, f64> *corrData);


//TODO: add doc
void printCorrelationMatrix(const CMF &protoGraph);


//TODO: add doc
bool isProtoGraphValid(const struct correlationMatrix &protoGraph);


//TODO: add doc
void graphvizRepresentation(graph<geneData, f64> *corrData, vector<string> labels);


//TODO: add doc
void printEdgeWeights(graph<geneData, u8> *corrData);


//TODO: add doc
void printProtoGraph(const CMF &toPrint);


/*******************************************************************//**
 * print graph to make sense of it's contents to stderr.  The graph is
 * printed with the vertex name on a line, followed by a pair of square
 * brackets ('[ ' ' ]')  which contain curly bracket pairs which hold
 * the name of a connected vertex followed by the edge weight.  Each of
 * these entries is followed by a blank line.
 *
 * @param[in] toPrint The graph structure to print
 * @param[in] labels
 **********************************************************************/
void printGraph(graph<geneData, f64> *toPrint,
                                          const vector<string> &labels);


////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


#endif
