#include "gEdge.hpp"
#include "gNode.hpp"


gEdge::gEdge(gNode *newLeft, gNode *newRight, size_t myID){
  left = newLeft;
  right = newRight;

  leftEdgeIndex = left->addEdge(this);
  rightEdgeIndex = right->addEdge(this);

  edgeID = myID;
}


gEdge::~gEdge(){
  left->removeEdge(this);
  right->removeEdge(this);
}