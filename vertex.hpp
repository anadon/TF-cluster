/***********************************************************************

***********************************************************************/
#ifndef GNODE_HPP
#define GNODE_HPP

#include <unistd.h>
#include <string>


template <typename T, typename U> class graph;

template <typename T, typename U> class edge;


template <typename T, typename U> class vertex{
  private:
  edge<T, U> **edges;
  size_t numEdges, edgesSize;
  
  public:
  size_t vertexIndex;
  T value;
  
  vertex(size_t index, T data);

  ~vertex();

  size_t addEdge(edge<T, U> *toRegister);
  
  edge<T, U>** getEdges();
  
  const edge<T, U>** getEdges() const;
  
  size_t getNumEdges() const;

  void removeEdge(edge<T, U> *toRemove);
  
  //void updateIndex(size_t newIndex);
  
  void clear();
  
  void hintNumEdges(const size_t suggestSize);
  
  void shrinkEdgeCapacityToFit();
  
  bool operator==(const vertex<T, U> &other) const;
  
  private:
  
  void ensureEdgeCapacity(const size_t size);
  
};

#endif