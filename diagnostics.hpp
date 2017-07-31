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

#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP


////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <queue>
#include <string>
#include <vector>

#include "auxiliaryUtilities.hpp"
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


void printCoincidenceMatrix(UpperDiagonalSquareMatrix<u8> *matrix,
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
void printGraph(graph<geneData, f64> *toPrint, const vector<string> &labels);



void printGraphTopOneHundred(graph<geneData, f64> *toPrint, const vector<string> &labels);


////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


#endif
