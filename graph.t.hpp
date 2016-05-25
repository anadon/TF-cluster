#include <unistd.h>
#include <stdlib.h>

#include "graph.hpp"

using std::unordered_map;


template <typename T, typename U> graph<T, U>::graph(){
  numNodes = numEdges = 0;
  vertexArray = (vertex<T, U>**) NULL;
  edgeArray = (edge<T, U>**) NULL;
}


template <typename T, typename U> graph<T, U>::~graph(){
  for(size_t i = 0; i < numEdges; i++)
    removeEdge(edgeArray[i]);
  for(size_t i = 0; i < numNodes; i++)
    removeVertex(vertexArray[i]);
}


template <typename T, typename U> void graph<T, U>::removeEdge(edge<T, U> *toRemove){
  const size_t edgeIndex = toRemove->edgeID;

  if(toRemove != edgeArray[edgeIndex]){
    fprintf(stderr, "Error: trying to remove edge that doesn't belong "
        "to this network!\n");
    exit(1);
  }

  delete toRemove;
  edgeArray[edgeIndex] = edgeArray[--numEdges];
  edgeArray[edgeIndex]->edgeID = edgeIndex;

  edgeArray = (edge<T, U>**) realloc(edgeArray, sizeof(*edgeArray) * numEdges);
}


template <typename T, typename U> void graph<T, U>::removeVertex(vertex<T, U> *toRemove){
  const size_t nodeIndex = toRemove->nodeID;

  if(toRemove != vertexArray[nodeIndex]){
    fprintf(stderr, "Error: trying to remove vertex that doesn't belong"
        "to this network!\n");
    exit(1);
  }

  for(size_t i = vertexArray[nodeIndex]->numEdges - 1; i < vertexArray[nodeIndex]->numEdges; i--){
    removeEdge(vertexArray[nodeIndex]->edges[i]);
  }


  delete toRemove;
  vertexArray[nodeIndex] = vertexArray[--numNodes];
  vertexArray[nodeIndex]->nodeID = nodeIndex;

  vertexArray = (vertex<T, U>**) realloc(vertexArray, sizeof(*vertexArray) * numNodes);
}


template <typename T, typename U> void graph<T, U>::addVertex(T data){
  size_t newSize = numNodes + 1;
  vertexArray = (vertex<T, U>**) realloc(vertexArray, sizeof(*vertexArray) * newSize);
  geneNameToNodeID.emplace(data, numNodes);
  vertexArray[numNodes] = new vertex<T, U>(numNodes);

  numNodes = newSize;
}


template <typename T, typename U> void graph<T, U>::addEdge(vertex<T, U> *left, vertex<T, U> *right, U newWeight){

  typename unordered_map<T, size_t>::const_iterator potentialFind;
  potentialFind = geneNameToNodeID.find(right->value);
  if(potentialFind == geneNameToNodeID.end()) return;
  //TODO if found, exit early

  size_t newSize = numEdges + 1;
  edgeArray = (edge<T, U>**) realloc(vertexArray, sizeof(*vertexArray) * newSize);
  geneNameToNodeID.emplace(left->value, numEdges);
  edgeArray[numNodes] = new edge<T, U>(left, right, newWeight);
}


template <typename T, typename U> const edge<T, U>** graph<T, U>::getEdges(){
  return (const edge<T, U>**) edgeArray;
}

template <typename T, typename U> const size_t graph<T, U>::getNumEdges(){
  return numEdges;
}