/***********************************************************************

***********************************************************************/
#ifndef GEDGE_HPP
#define GEDGE_HPP

#include <unistd.h>

#include "auxillaryUtilities.hpp"


template <typename T, typename U> class vertex;

template <typename T, typename U> class edge{
  public:
  vertex<T, U> *left, *right;
  U weight;
  size_t edgeID, leftEdgeIndex, rightEdgeIndex;

  edge(vertex<T, U> *newLeft, vertex<T, U> *newRight, U newWeight, 
                                                      size_t edgeIndex);

  ~edge();
  
  vertex<T, U>* other(const vertex<T, U> *side);
};

#endif