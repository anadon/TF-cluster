/***********************************************************************

***********************************************************************/
#ifndef NODENETWORK_HPP
#define NODENETWORK_HPP

#include <mutex>
#include <stdlib.h>
#include <unistd.h>
#include <unordered_map>

#include "edge.hpp"
#include "geneData.hpp"
#include "vertex.hpp"

using std::unordered_map;
using std::mutex;


template <typename T, typename U> class graph{
  private:
  vertex<T, U> **vertexArray;
  edge<T, U> **edgeArray;
  size_t numVertexes, numEdges;
  size_t vertexArraySize, edgeArraySize;
  unordered_map<T, size_t, std::hash<T>> geneNameToNodeID;
  //std::mutex edgeGuard, vertexGuard;

/*GRAPH OPERATIONS*****************************************************/
  public:
  
  graph();
  ~graph();
  
  graph<T, U>& operator=(const graph<T, U> &other);
  
/*EDGE OPERATIONS******************************************************/
  public:
  edge<T, U>* addEdge(vertex<T, U> *left, vertex<T, U> *right, U newWeight);

  edge<T, U>* addEdge(edge<T, U> &toAdd);

  edge<T, U>** getEdges();

  const edge<T, U>** getEdges() const;

  size_t getNumEdges() const;
  
  void hintNumEdges(const size_t suggestSize);

  U removeEdge(edge<T, U> *toRemove);
  
  void shrinkEdgeCapacityToFit();
  
  private:
  void ensureEdgeCapacity(const size_t size);
  
/*VERTEX OPERATIONS****************************************************/
  public:
  vertex<T, U>* addVertex(T newValue);

  vertex<T, U>* addVertex(vertex<T, U> *newVertex);

  vertex<T, U>* getVertexForValue(const T &testValue);
  
  vertex<T, U>** getVertexes();

  size_t getNumVertexes() const;
  
  void hintNumVertexes(const size_t suggestSize);

  T removeVertex(const vertex<T, U> *toRemove);
  
  //T removeVertex(T value);
  
  void shrinkVertexCapacityToFit();
  
  private:
  void ensureVertexCapacity(const size_t size);
  
};

#endif