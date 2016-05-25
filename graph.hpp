/***********************************************************************

***********************************************************************/
#ifndef NODENETWORK_HPP
#define NODENETWORK_HPP

#include <stdlib.h>
#include <unistd.h>
#include <unordered_map>

#include "edge.hpp"
#include "vertex.hpp"

using std::unordered_map;


template <typename T, typename U> class graph{
  public:
  vertex<T, U> **vertexArray;
  edge<T, U> **edgeArray;
  size_t numNodes, numEdges;
  unordered_map<T, size_t> geneNameToNodeID;


  graph();
  ~graph();

  void addVertex(T newValue);

  void addEdge(vertex<T, U> *left, vertex<T, U> *right, U newWeight);

  void removeEdge(edge<T, U> *toRemove);

  void removeVertex(vertex<T, U> *toRemove);

  const edge<T, U>** getEdges();

  const size_t getNumEdges();
};

#endif