
#include <iostream>
#include <stdlib.h>
#include <unistd.h>

#include "geneData.hpp"
#include "graph.hpp"

using std::unordered_map;
using std::cout;
using std::endl;


template <typename T, typename U> graph<T, U>::graph(){
  numVertexes = numEdges = 0;
  vertexArray = (vertex<T, U>*) NULL;
  edgeArray = (edge<T, U>*) NULL;
}


template <typename T, typename U> graph<T, U>::~graph(){
  for(size_t i = 0; i < numEdges; i++)
    removeEdge(&edgeArray[i]);
  for(size_t i = 0; i < numVertexes; i++)
    removeVertex(&vertexArray[i]);
}


template <typename T, typename U> void graph<T, U>::removeEdge(edge<T, U> *toRemove){
  const size_t edgeIndex = toRemove->edgeID;
  edge<T, U> *memCheck;

  edgeGuard.lock();
  if(toRemove != &edgeArray[edgeIndex]){
    fprintf(stderr, "Error: trying to remove edge that doesn't belong "
        "to this network!\n");
    exit(1);
  }

  delete toRemove;
  edgeArray[edgeIndex] = edgeArray[--numEdges];
  edgeArray[edgeIndex].edgeID = edgeIndex;

  memCheck = (edge<T, U>*) realloc(edgeArray, sizeof(*edgeArray) * numEdges);
  if(NULL == memCheck){
    exit(1);
  }
  edgeArray = memCheck;
  edgeGuard.unlock();
}


template <typename T, typename U> void graph<T, U>::removeVertex(vertex<T, U> *toRemove){
  const size_t nodeIndex = toRemove->nodeID;
  vertex<T, U> *memCheck;

  vertexGuard.lock();
  if(toRemove != &vertexArray[nodeIndex]){
    fprintf(stderr, "Error: trying to remove vertex that doesn't belong"
        "to this network!\n");
    exit(1);
  }

  for(size_t i = vertexArray[nodeIndex].numEdges - 1; i < vertexArray[nodeIndex].numEdges; i--){
    removeEdge(vertexArray[nodeIndex].edges[i]);
  }


  vertexArray[nodeIndex] = vertexArray[--numVertexes];
  vertexArray[nodeIndex].nodeID = nodeIndex;
  delete toRemove;

  memCheck = (vertex<T, U>*) realloc(vertexArray, sizeof(*vertexArray) * numVertexes);
  
  if(NULL == memCheck){
    exit(1);
  }
  
  vertexArray = memCheck;
  vertexGuard.unlock();
}


template <typename T, typename U> void graph<T, U>::addVertex(T data){
  vertex<T, U> *memCheck;
  
  vertexGuard.lock();
  size_t newSize = numVertexes + 1;
  
  memCheck = (vertex<T, U>*) realloc(vertexArray, sizeof(*vertexArray) * newSize);
  
  if(NULL == memCheck){
    exit(1);
  }
  vertexArray = memCheck;
  
  geneNameToNodeID.emplace(data, numVertexes);
  vertexArray[numVertexes] = vertex<T, U>(numVertexes, data);

  numVertexes = newSize;
  vertexGuard.unlock();
}


template <typename T, typename U> void graph<T, U>::addEdge(vertex<T, U> *left, vertex<T, U> *right, U newWeight){
  typename unordered_map<T, size_t>::const_iterator potentialFind;
  edge<T, U> *memCheck;
  
  potentialFind = geneNameToNodeID.find(left->value);
  if(potentialFind == geneNameToNodeID.end()){
    cout << "Could not locate left vertex in graph" << endl;
    return;
  }
  
  potentialFind = geneNameToNodeID.find(right->value);
  if(potentialFind == geneNameToNodeID.end()){
    cout << "Could not locate right vertex in graph" << endl;
    return;
  }
  
  edgeGuard.lock();
  size_t newSize = numEdges + 1;
  memCheck = (edge<T, U>*) realloc(edgeArray, sizeof(*edgeArray) * newSize);
  
  if(NULL == memCheck){
    exit(1);
  }
  
  edgeArray = memCheck;
  
  edgeArray[numEdges] = edge<T, U>(left, right, newWeight);
  
  numEdges = newSize;
  edgeGuard.unlock();
}


template <typename T, typename U> void graph<T, U>::addEdge(const edge<T, U> &toAdd){
  typename unordered_map<T, size_t>::const_iterator findLeft, findRight, endItr;
  edge<T, U> memCheck;
  
  findLeft = geneNameToNodeID.find(toAdd.left->value);
  findRight = geneNameToNodeID.find(toAdd.right->value);
  endItr = geneNameToNodeID.end();
  
  if(findLeft == endItr || findRight == endItr) exit(1);
  
  edgeGuard.lock();
  size_t newSize = numEdges + 1;
  memCheck = (edge<T, U>*) realloc(edgeArray, sizeof(*edgeArray) * newSize);
  
  if(NULL == memCheck){
    exit(1);
  }
  
  edgeArray = memCheck;
  
  edgeArray[numEdges] = edge<T, U>(&vertexArray[findLeft.second], &vertexArray[findRight.second], toAdd.weight);
  
  numEdges = newSize;
  edgeGuard.unlock();
}


template <typename T, typename U> const edge<T, U>* graph<T, U>::getEdges(){
  return (const edge<T, U>*) edgeArray;
}


template <typename T, typename U> size_t graph<T, U>::getNumEdges() const{
  return numEdges;
}


template <typename T, typename U> const vertex<T, U>* graph<T, U>::getVertexes(){
  return (const vertex<T, U>*) vertexArray;
}


template <typename T, typename U> size_t graph<T, U>::getNumVertexes() const{
  return numVertexes;
}


template <typename T, typename U> void graph<T, U>::addVertex(const vertex<T, U> &newVertex){
  vertex<T, U> *memCheck;
  
  vertexGuard.lock();
  const size_t newSize = numVertexes + 1;
  
  memCheck = (vertex<T, U>*) realloc(vertexArray, sizeof(*vertexArray) * newSize);
  
  if(NULL == memCheck){
    exit(1);
  }
  
  vertexArray = memCheck;
  
  geneNameToNodeID.emplace(newVertex.value, numVertexes);
  vertexArray[numVertexes] = vertex<T, U>(numVertexes, newVertex.value);
  
  numVertexes = newSize;
  vertexGuard.unlock();
}


template <typename T, typename U> graph<T, U>& graph<T, U>::operator=(const graph<T, U> &other){
  graph<T, U> toReturn;
  
  edgeGuard.lock();
  vertexGuard.lock();
  
  for(size_t i = 0; i < this->numVertexes; i++)
    toReturn.addVertex(this->vertexArray[i]);
  
  for(size_t i = 0; i < this->numEdges; i++){
    toReturn.addEdge(this->edgeArray[i]);
  }
  
  edgeGuard.unlock();
  vertexGuard.unlock();
  
  return toReturn;
}


template <typename T, typename U> vertex<T, U>* graph<T, U>::getVertexForValue(const T testValue){
  typename unordered_map<T, size_t>::const_iterator found;
  vertex<T, U>* tr;
  
  tr = (vertex<T, U>*) NULL;
  
  vertexGuard.lock();
  found = geneNameToNodeID.find(testValue);
  
  if(geneNameToNodeID.end() != found) tr = &vertexArray[found->second];
  vertexGuard.unlock();
  
  return tr;
}
