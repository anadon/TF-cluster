/***********************************************************************

***********************************************************************/
#ifndef GNODE_HPP
#define GNODE_HPP

#include <unistd.h>
#include <string>


template <typename T, typename U> class graph;

template <typename T, typename U> class edge;


template <typename T, typename U> class vertex{
  public:
  size_t nodeID;
  size_t numEdges;
  edge<T, U> **edges;
  
  T value;

  vertex(size_t ID);

  vertex(size_t ID, T value);

  ~vertex();

  size_t addEdge(edge<T, U> *toRegister);

  void removeEdge(edge<T, U> *toRemove);
  
  bool operator==(const vertex<T, U> &other) const;
  
};

#endif