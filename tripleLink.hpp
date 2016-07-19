#ifndef TRIPLELINK_HPP
#define TRIPLELINK_HPP

#include <vector>

#include "geneData.hpp"
#include "graph.hpp"

using std::vector;

queue<size_t> tripleLinkIteration(graph<geneData, double> *geneNetwork,
                    cf64 tripleLink1, cf64 tripleLink2, 
                                          const size_t targetEdgeIndex);

queue< queue<size_t> > tripleLink(graph<geneData, double> *geneNetwork,
                    cf64 tripleLink1, cf64 tripleLink2);

#endif