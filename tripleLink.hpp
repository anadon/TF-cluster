#ifndef TRIPLELINK_HPP
#define TRIPLELINK_HPP

#include <vector>

#include "geneData.hpp"
#include "graph.hpp"

using std::vector;

graph<geneData, double>* tripleLinkIteration(graph<geneData, double> *geneNetwork,
                    const double tripleLink1, const double tripleLink2, 
                const double tripleLink3, const size_t targetEdgeIndex);

vector< graph<geneData, double>* > tripleLink(graph<geneData, double> *geneNetwork,
                    const double tripleLink1, const double tripleLink2, 
                                              const double tripleLink3);

#endif