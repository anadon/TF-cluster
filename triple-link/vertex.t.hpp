#include <cstring>

#include "vertex.hpp"
#include "edge.hpp"


template <typename T, typename U> vertex<T, U>::vertex(size_t ID){
  nodeID = ID;
  numEdges = 0;
  edges = (edge<T, U>**) NULL;
}


template <typename T, typename U> vertex<T, U>::vertex(size_t ID, T newValue){
  nodeID = ID;
  numEdges = 0;
  edges = (edge<T, U>**) NULL;
  value = newValue;
}


template <typename T, typename U> vertex<T, U>::~vertex(){
  for(size_t i = numEdges-1; i < numEdges; i--)
    removeEdge(edges[i]);
}


template <typename T, typename U> size_t vertex<T, U>::addEdge(edge<T, U> *toRegister){
  size_t toReturn, newCount;
  edge<T, U> **memCheck;

  toReturn = numEdges;
  newCount = numEdges + 1;

  memCheck = (edge<T, U>**) realloc(edges, newCount * sizeof(*edges));
  
  if(NULL == memCheck){
    exit(1);
  }
  edges = memCheck;
  
  edges[numEdges] = toRegister;
  numEdges = newCount;

  return toReturn;
}


template <typename T, typename U> void vertex<T, U>::removeEdge(edge<T, U> *toRemove){
  edge<T, U> *tmp, **memCheck;
  size_t targetEdgeIndex;
  
  
  numEdges--;

  //check if this is the right or left of the edge, error if edge
  //doesn't connect this node/vertex
  if(toRemove == edges[toRemove->leftEdgeIndex]){
    targetEdgeIndex = toRemove->leftEdgeIndex;
  }else if(toRemove == edges[toRemove->rightEdgeIndex]){
    targetEdgeIndex = toRemove->rightEdgeIndex;
  }else{
    fprintf(stderr, "vertex::removeEdge: passed edges to remove is not"
        " connected to this vertex\n");
    fflush(stderr);
    exit(1);
  }

  //swap edge this is removing to the end
  tmp = edges[numEdges];
  edges[numEdges] = edges[targetEdgeIndex];
  edges[targetEdgeIndex] = tmp;

  //update the swapped edge's location in this structure so it still
  //knows where it is in this vertex/node
  if(this == edges[targetEdgeIndex]->left)
    edges[targetEdgeIndex]->leftEdgeIndex = targetEdgeIndex;
  else
    edges[targetEdgeIndex]->rightEdgeIndex = targetEdgeIndex;

  //deallocate the last edge pointer (but don't actually delete --
  //that's the network's job).
  memCheck = (edge<T, U>**) realloc(edges, numEdges * sizeof(*edges));
  
  if(NULL == memCheck){
    exit(1);
  }
  edges = memCheck;
}


template <typename T, typename U> inline bool vertex<T, U>::operator==(const vertex<T, U> &other) const{
  return (value == other.value);
}
