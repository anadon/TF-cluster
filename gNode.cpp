
#include "gNode.hpp"
#include "gEdge.hpp"
#include "nodeNetwork.hpp"


gNode::gNode(nodeNetwork *belongTo, size_t ID){
  nodeID = ID;
  net = belongTo;
  numEdges = 0;
  clusterGroup = 0;
  edges = (gEdge**) NULL;
}


gNode::~gNode(){
  for(size_t i = 0; i < numEdges; i++)
    net->removeEdge(edges[i]->edgeID);
}


size_t gNode::addEdge(gEdge *toRegister){
  size_t toReturn, newCount;

  toReturn = numEdges;
  newCount = numEdges + 1;

  edges = (gEdge**) realloc(edges, newCount * sizeof(*edges));
  edges[numEdges] = toRegister;
  numEdges = newCount;

  return toReturn;
}


void gNode::removeEdge(size_t edgeIndex){
  numEdges--;
  gEdge *tmp;

  tmp = edges[numEdges];
  edges[numEdges] = edges[edgeIndex];
  edges[edgeIndex] = tmp;
  if(this == edges[edgeIndex]->left)
    edges[edgeIndex]->leftEdgeIndex = edgeIndex;
  else
    edges[edgeIndex]->rightEdgeIndex = edgeIndex;

  edges = (gEdge**) realloc(edges, numEdges * sizeof(*edges));
}


void gNode::prune(u8 keepTopN){
  for(size_t i = keepTopN; i < numEdges; i++){
    net->removeEdge(edges[i]->edgeID);
    //edges[i] = (gNode*) NULL;
  }
}