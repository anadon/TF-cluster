/***********************************************************************

***********************************************************************/
#ifndef GNODE_HPP
#define GNODE_HPP

#include <unistd.h>
#include <string>

#include "auxillaryUtilities.hpp"


class nodeNetwork;

class gEdge;


class gNode{
  public:
  size_t nodeID;
  size_t numEdges;
  size_t clusterGroup;
  gEdge **edges;
  string geneID;
  nodeNetwork *net;

  gNode(nodeNetwork *belongTo, size_t ID);

  ~gNode();

  size_t addEdge(gEdge *toRegister);

  void removeEdge(gEdge *toRemove);

  void prune(u8 keepTopN);
  
  private:
  void quickMergeEdges(gEdge **toSort, const size_t size);
  
  void mergeHelper(gEdge **toSort, const size_t leftIndex, 
                        const size_t rightIndex, const size_t endIndex);
  
};

#endif