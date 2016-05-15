#ifndef NODENETWORK_HPP
#define NODENETWORK_HPP

#include <stdlib.h>
#include <unistd.h>
#include <unordered_map>

#include "gEdge.hpp"
#include "gNode.hpp"

using std::unordered_map;

class nodeNetwork{
  private:
  gNode **nodeArray;
  gEdge **edgeArray;
  size_t numNodes, numEdges;
  unordered_map<string, size_t> geneNameToNodeID;
  
  public:
  
  void addNode(string geneID);
  
  void addEdge(gNode *left, gNode *right, double coeff);
  
  void removeEdge(size_t edgeIndex);
  
  void removeNode(size_t nodeIndex);
  
  
  nodeNetwork();
  ~nodeNetwork();
  
  int loadFromFile(FILE *input);
  
  void prune(u8 keepTopN = 100);
};

#endif