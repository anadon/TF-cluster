/***********************************************************************

***********************************************************************/
#ifndef GEDGE_HPP
#define GEDGE_HPP

#include <unistd.h>

#include "auxillaryUtilities.hpp"


template <typename T, typename U> class vertex;

template <typename T, typename U> class edge{
  public:
  size_t leftEdgeIndex, rightEdgeIndex, edgeID;
  vertex<T, U> *left, *right;
  U weight;

  edge(vertex<T, U> *newLeft, vertex<T, U> *newRight, U newWeight);

  ~edge();
  
  vertex<T, U>* other(const vertex<T, U> *side);
};

#endif