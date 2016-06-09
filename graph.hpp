/***********************************************************************

***********************************************************************/
#ifndef NODENETWORK_HPP
#define NODENETWORK_HPP

#include <mutex>
#include <stdlib.h>
#include <unistd.h>
#include <unordered_map>

#include "edge.hpp"
#include "vertex.hpp"

using std::unordered_map;
using std::mutex;


template <typename T, typename U> class graph{
  public:
  vertex<T, U> *vertexArray;
  edge<T, U> *edgeArray;
  size_t numVertexes, numEdges;
  unordered_map<T, size_t> geneNameToNodeID;
  mutex edgeGuard, vertexGuard;


  graph();
  ~graph();

  void addVertex(T newValue);

  void addVertex(const vertex<T, U> &newVertex);

  void addEdge(vertex<T, U> *left, vertex<T, U> *right, U newWeight);

  void addEdge(const edge<T, U> &toAdd);

  void removeEdge(edge<T, U> *toRemove);

  void removeVertex(vertex<T, U> *toRemove);
  
  vertex<T, U>* getVertexForValue(const T testValue);

  const edge<T, U>* getEdges();
  
  const vertex<T, U>* getVertexes();

  size_t getNumEdges() const;

  size_t getNumVertexes() const;
  
  graph<T, U>& operator=(const graph<T, U> &other);

  vertex<T, U>* getVertexForValue(const T testValue) const;

};

#endif