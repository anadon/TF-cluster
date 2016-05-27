#include <unistd.h>
#include <stdlib.h>

#include "graph.hpp"

using std::unordered_map;


template <typename T, typename U> graph<T, U>::graph(){
  numVertexes = numEdges = 0;
  vertexArray = (vertex<T, U>**) NULL;
  edgeArray = (edge<T, U>**) NULL;
}


template <typename T, typename U> graph<T, U>::~graph(){
  for(size_t i = 0; i < numEdges; i++)
    removeEdge(edgeArray[i]);
  for(size_t i = 0; i < numVertexes; i++)
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
  vertexArray[nodeIndex] = vertexArray[--numVertexes];
  vertexArray[nodeIndex]->nodeID = nodeIndex;

  vertexArray = (vertex<T, U>**) realloc(vertexArray, sizeof(*vertexArray) * numVertexes);
}


template <typename T, typename U> void graph<T, U>::addVertex(T data){
  size_t newSize = numVertexes + 1;
  vertexArray = (vertex<T, U>**) realloc(vertexArray, sizeof(*vertexArray) * newSize);
  geneNameToNodeID.emplace(data, numVertexes);
  vertexArray[numVertexes] = new vertex<T, U>(numVertexes);

  numVertexes = newSize;
}


template <typename T, typename U> void graph<T, U>::addEdge(vertex<T, U> *left, vertex<T, U> *right, U newWeight){
  typename unordered_map<T, size_t>::const_iterator potentialFind;
  potentialFind = geneNameToNodeID.find(right->value);
  if(potentialFind == geneNameToNodeID.end()) return;
  
  potentialFind = geneNameToNodeID.find(left->value);
  if(potentialFind == geneNameToNodeID.end()) return;

  size_t newSize = numEdges + 1;
  edgeArray = (edge<T, U>**) realloc(vertexArray, sizeof(*vertexArray) * newSize);
  geneNameToNodeID.emplace(left->value, numEdges);
  geneNameToNodeID.emplace(right->value, numEdges);
  edgeArray[numVertexes] = new edge<T, U>(left, right, newWeight);
}


template <typename T, typename U> void graph<T, U>::addEdge(const edge<T, U> &toAdd){
  typename unordered_map<T, size_t>::const_iterator findLeft, findRight, endItr;
  findLeft = geneNameToNodeID.find(toAdd.left->value);
  findRight = geneNameToNodeID.find(toAdd.right->value);
  endItr = geneNameToNodeID.end();
  
  if(findLeft == endItr || findRight == endItr) exit(1);

  size_t newSize = numEdges + 1;
  edgeArray = (edge<T, U>**) realloc(vertexArray, sizeof(*vertexArray) * newSize);
  edgeArray[numVertexes] = new edge<T, U>(vertexArray[findLeft.second], vertexArray[findRight.second], toAdd.weight);
}


template <typename T, typename U> const edge<T, U>** graph<T, U>::getEdges(){
  return (const edge<T, U>**) edgeArray;
}


template <typename T, typename U> const size_t graph<T, U>::getNumEdges() const{
  return numEdges;
}


template <typename T, typename U> const vertex<T, U>** graph<T, U>::getVertexes(){
  return (const vertex<T, U>**) vertexArray;
}


template <typename T, typename U> const size_t graph<T, U>::getNumVertexes() const{
  return numVertexes;
}


template <typename T, typename U> void graph<T, U>::addVertex(const vertex<T, U> &newVertex){
  const size_t newSize = numVertexes + 1;
  vertexArray = (vertex<T, U>**) realloc(vertexArray, sizeof(*vertexArray) * newSize);
  geneNameToNodeID.emplace(newVertex.value, numVertexes);
  vertexArray[numVertexes] = new vertex<T, U>(numVertexes, newVertex.value);
  
  numVertexes = newSize;
}


template <typename T, typename U> graph<T, U>& graph<T, U>::operator=(const graph<T, U> &other) const{
  graph<T, U> toReturn;
  for(size_t i = 0; i < this->numVertexes; i++)
    toReturn.addVertex(*(this->vertexArray[i]));
  
  for(size_t i = 0; i < this->numEdges; i++){
    toReturn.addEdge(*(this->edgeArray[i]));
  }
  
  return toReturn;
}


template <typename T, typename U> vertex<T, U>* graph<T, U>::getVertexForValue(const T testValue){
  return vertexArray[geneNameToNodeID.at(testValue)];
}
