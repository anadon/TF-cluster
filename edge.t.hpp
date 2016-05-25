#include "edge.hpp"
#include "vertex.hpp"


template <typename T, typename U> edge<T, U>::edge(vertex<T, U> *newLeft, vertex<T, U> *newRight, U newWeight){
  left = newLeft;
  right = newRight;

  leftEdgeIndex = left->addEdge(this);
  rightEdgeIndex = right->addEdge(this);

  weight = newWeight;
}


template <typename T, typename U> edge<T, U>::~edge(){
  left->removeEdge(this);
  right->removeEdge(this);
}