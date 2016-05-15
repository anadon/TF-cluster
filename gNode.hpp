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

  void removeEdge(size_t edgeIndex);

  void prune(u8 keepTopN);
  
};

#endif