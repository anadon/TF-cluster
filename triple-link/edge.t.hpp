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


template <typename T, typename U> vertex<T, U>* edge<T, U>::other(const vertex<T, U> *side){
  if(*side == *left){
    return right;
  }else if(*side == *right){
    return left;
  }
  
  return ( vertex<T, U>* ) NULL;
}