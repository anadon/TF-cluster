#include <csignal>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>

#include "geneData.hpp"
#include "graph.hpp"

using std::unordered_map;
using std::cout;
using std::endl;


template <typename T, typename U> graph<T, U>::graph(){
  numVertexes = numEdges = vertexArraySize = edgeArraySize = 0;
  vertexArray = (vertex<T, U>**) NULL;
  edgeArray = (edge<T, U>**) NULL;
}


template <typename T, typename U> graph<T, U>::~graph(){
  for(size_t i = 0; i < numEdges; i++)
     delete edgeArray[i];
  free(edgeArray);
  
  for(size_t i = 0; i < numVertexes; i++){
    vertexArray[i]->clear();//TODO once debugging is complete, wrap clear into ~?
    delete vertexArray[i];
  }
  free(vertexArray);
  
  geneNameToNodeID.clear();
}


template <typename T, typename U> U graph<T, U>::removeEdge(edge<T, U> *toRemove){
  const size_t edgeIndex = toRemove->edgeID;
  void *memCheck;
  U tr;

  if(edgeIndex >= numEdges) raise(SIGABRT);
  if(toRemove != edgeArray[edgeIndex])  raise(SIGABRT);
  
  tr = edgeArray[edgeIndex]->weight;
  --numEdges;
  delete edgeArray[edgeIndex];
  edgeArray[edgeIndex] = edgeArray[numEdges];
  edgeArray[edgeIndex]->edgeID = edgeIndex;

  memCheck = realloc(edgeArray, sizeof(*edgeArray) * numEdges);
  edgeArray = (edge<T, U>**) memCheck;
  
  return tr;
}


template <typename T, typename U> T graph<T, U>::removeVertex(vertex<T, U> *toRemove){
  const size_t nodeIndex = toRemove->vertexIndex;
  void *memCheck;
  T tr;
  
  if(nodeIndex >= numVertexes)  raise(SIGABRT);
  if(toRemove != vertexArray[nodeIndex])  raise(SIGABRT);

  while(toRemove->getNumEdges())
    removeEdge(toRemove->getEdges()[toRemove->getNumEdges()-1]);

  tr = vertexArray[nodeIndex]->value;
  
  geneNameToNodeID.erase(tr);
  delete vertexArray[nodeIndex];
  vertexArray[nodeIndex] = vertexArray[numVertexes-1];
  vertexArray[nodeIndex]->updateIndex(nodeIndex);

  --numVertexes;
  memCheck = realloc(vertexArray, sizeof(*vertexArray) * numVertexes);
  vertexArray = (vertex<T, U>**) memCheck;
  
  return tr;
}


template <typename T, typename U> T graph<T, U>::removeVertex(T value){
  return removeVertex(getVertexForValue(value));
}


template <typename T, typename U> vertex<T, U>* graph<T, U>::addVertex(T data){
  vertex<T, U> *tr;
  
  tr = getVertexForValue(data);
  if(NULL != tr) return tr;
  
  ensureVertexCapacity(numVertexes + 1);
  
  geneNameToNodeID.emplace(data, numVertexes);
  vertexArray[numVertexes] = new vertex<T, U>(numVertexes, data);
  
  numVertexes++;
  return vertexArray[numVertexes-1];
}


template <typename T, typename U> edge<T, U>* graph<T, U>::addEdge(
                  vertex<T, U> *left, vertex<T, U> *right, U newWeight){
  
  if(0 == geneNameToNodeID.count(left->value))  raise(SIGABRT);  
  if(0 == geneNameToNodeID.count(right->value)) raise(SIGABRT);
  
  ensureEdgeCapacity(numEdges + 1);
  
  edgeArray[numEdges] = new edge<T, U>(left, right, newWeight, numEdges);
  
  numEdges++;
  return edgeArray[numEdges-1];
}


template <typename T, typename U> edge<T, U>* graph<T, U>::addEdge(
                                                    edge<T, U> &toAdd){
  
  return addEdge(vertexArray[geneNameToNodeID[toAdd.left->value]], 
                      vertexArray[geneNameToNodeID[toAdd.right->value]], 
                                                          toAdd.weight);
}


template <typename T, typename U> edge<T, U>** graph<T, U>::getEdges(){
  return (edge<T, U>**) edgeArray;
}


template <typename T, typename U> size_t graph<T, U>::getNumEdges() const{
  return numEdges;
}


template <typename T, typename U> vertex<T, U>** graph<T, U>::getVertexes(){
  return (vertex<T, U>**) vertexArray;
}


template <typename T, typename U> size_t graph<T, U>::getNumVertexes() const{
  return numVertexes;
}


template <typename T, typename U> vertex<T, U>* graph<T, U>::addVertex(vertex<T, U> *newVertex){
  return addVertex(newVertex->value);
}


template <typename T, typename U> graph<T, U>& graph<T, U>::operator=(const graph<T, U> &other){

  raise(SIGABRT);

  graph<T, U> toReturn;
  
  for(size_t i = 0; i < this->numVertexes; i++)
    toReturn.addVertex(*this->vertexArray[i]);
  
  for(size_t i = 0; i < this->numEdges; i++)
    toReturn.addEdge(*this->edgeArray[i]);
  
  return toReturn;
}


template <typename T, typename U> vertex<T, U>* graph<T, U>::getVertexForValue(const T testValue){
  if(geneNameToNodeID.count(testValue)) 
    return vertexArray[geneNameToNodeID[testValue]];
  return NULL;
}


template <typename T, typename U> void graph<T, U>::hintNumEdges(const size_t suggestSize){
  void *memCheck;
  
  if(suggestSize <= numEdges) return;
  
  if(suggestSize == edgeArraySize) return;
    
  memCheck = realloc(edgeArray, suggestSize * sizeof(*edgeArray));
  if(NULL != memCheck){
    edgeArray = (edge<T, U>**) memCheck;
    edgeArraySize = suggestSize;
  }else{
    raise(SIGABRT);
  }
}


template <typename T, typename U> void graph<T, U>::hintNumVertexes(const size_t suggestSize){
  void *memCheck;
  
  if(suggestSize <= numVertexes) return;
  
  if(suggestSize == vertexArraySize) return;
  
  geneNameToNodeID.reserve(suggestSize);
    
  memCheck = realloc(vertexArray, suggestSize * sizeof(*vertexArray));
  if(NULL != memCheck){
    vertexArray = (vertex<T, U>**) memCheck;
    vertexArraySize = suggestSize;
  }else{
    raise(SIGABRT);
  }
}


template <typename T, typename U> void graph<T, U>::shrinkEdgeCapacityToFit(){
  hintNumEdges(numEdges);
}
  
  
template <typename T, typename U> void graph<T, U>::shrinkVertexCapacityToFit(){
  hintNumVertexes(numVertexes);
}


template <typename T, typename U> void graph<T, U>::ensureEdgeCapacity(const size_t size){
  void *memCheck;
  size_t nextSize;
  
  if(size <= edgeArraySize) return;
  
  nextSize = 1 + (edgeArraySize << 1);
  
  memCheck = realloc(edgeArray, nextSize * sizeof(*edgeArray));
  if(NULL == memCheck){
    raise(SIGABRT);
  }
  edgeArray = (edge<T, U>**) memCheck;
  edgeArraySize = nextSize;
  
}
  
  
template <typename T, typename U> void graph<T, U>::ensureVertexCapacity(const size_t size){
  void *memCheck;
  size_t nextSize;
  if(size <= vertexArraySize) return;
  
  nextSize = 1 + (vertexArraySize << 1);
  
  memCheck = realloc(vertexArray, nextSize * sizeof(*vertexArray));
  if(NULL == memCheck){
    raise(SIGABRT);
  }
  vertexArray = (vertex<T, U>**) memCheck;
  vertexArraySize = nextSize;
  
}