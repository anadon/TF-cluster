#ifndef TRIPLELINK_HPP
#define TRIPLELINK_HPP

#include <vector>

#include "geneData.hpp"
#include "graph.hpp"

using std::vector;

graph<struct geneData, double>* tripleLinkIteration(graph<struct geneData, double> *geneNetwork,
            double tripleLink1, double tripleLink2, double tripleLink3, 
                                                size_t targetEdgeIndex);

vector< graph<struct geneData, double>* > tripleLink(graph<struct geneData, double> *geneNetwork,
            double tripleLink1, double tripleLink2, double tripleLink3);

#endif