/*Copyright 2016-2017 Josh Marshall************************************/

/***********************************************************************
    This file is part of TF-Cluster.

    TF-Cluster is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TF-Cluster.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#ifndef STD_GRAPH_T_HPP
#define STD_GRAPH_T_HPP

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <csignal>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>

#include "geneData.hpp"
#include "graph.hpp"

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::unordered_map;
using std::cout;
using std::endl;

////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

template <typename T, typename U> graph<T, U>::graph(){
  numVertexes = numEdges = vertexArraySize = edgeArraySize = 0;
  vertexArray = (vertex<T, U>**) NULL;
  edgeArray = (edge<T, U>**) NULL;
}


template <typename T, typename U> graph<T, U>::~graph(){

  //while(numEdges)  removeEdge(edgeArray[numEdges-1]);

  while(numVertexes) removeVertex(vertexArray[numVertexes-1]);

  geneNameToNodeID.clear();
}


template <typename T, typename U> U graph<T, U>::removeEdge(
                                                  edge<T, U> *toRemove){
  const size_t edgeIndex = toRemove->edgeID;
  void *memCheck;
  U tr;

  if(edgeIndex >= numEdges) raise(SIGABRT);
  if(toRemove != edgeArray[edgeIndex])  raise(SIGABRT);

  tr = edgeArray[edgeIndex]->weight;
  --numEdges;
  edgeArray[numEdges]->edgeID = edgeIndex;
  delete edgeArray[edgeIndex];
  edgeArray[edgeIndex] = edgeArray[numEdges];
  if(0 < numEdges){
    memCheck = realloc(edgeArray, sizeof(*edgeArray) * numEdges);
    edgeArray = (edge<T, U>**) memCheck;
  }else{
    free(edgeArray);
    edgeArray = NULL;
  }

  return tr;
}


//TODO: update vertex address in hash
template <typename T, typename U> T graph<T, U>::removeVertex(
                                          const vertex<T, U> *toRemove){
  const size_t nodeIndex = toRemove->vertexIndex;
  vertex<T, U> *target = vertexArray[nodeIndex];
  void *memCheck;


  if(NULL == toRemove)  raise(SIGABRT);
  if(nodeIndex >= numVertexes)  raise(SIGABRT);
  if(toRemove != target)  raise(SIGABRT);
  
  while(target->getNumEdges())
    removeEdge(target->getEdges()[target->getNumEdges()-1]);

  T tr = vertexArray[nodeIndex]->value;

  //Always remember, update THEN remove.
  geneNameToNodeID[vertexArray[numVertexes-1]->value] = nodeIndex;
  vertexArray[numVertexes-1]->vertexIndex = nodeIndex;
  vertexArray[nodeIndex] = vertexArray[numVertexes-1];

  geneNameToNodeID.erase(tr);
  delete target;

  --numVertexes;
  memCheck = realloc(vertexArray, sizeof(*vertexArray) * numVertexes);
  vertexArray = (vertex<T, U>**) memCheck;

  return tr;
}


/*template <typename T, typename U> T graph<T, U>::removeVertex(T value){
  return removeVertex(getVertexForValue(value));
}*/


template <typename T, typename U> vertex<T, U>* graph<T, U>::addVertex(
                                                                T data){
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

  //if(0 == geneNameToNodeID.count(left->value))  raise(SIGABRT);
  //if(0 == geneNameToNodeID.count(right->value)) raise(SIGABRT);
  //if(left == right) raise(SIGABRT);

  ensureEdgeCapacity(numEdges + 1);

  edgeArray[numEdges] = new edge<T, U>(left, right, newWeight,
                                                              numEdges);

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


template <typename T, typename U>const edge<T, U>**
                                          graph<T, U>::getEdges() const{
  return (const edge<T, U>**) edgeArray;
}


template <typename T, typename U> size_t graph<T, U>::getNumEdges() const{
  return numEdges;
}


template <typename T, typename U> vertex<T, U>**
                                            graph<T, U>::getVertexes(){
  return (vertex<T, U>**) vertexArray;
}


template <typename T, typename U> size_t graph<T, U>::getNumVertexes()
                                                                  const{
  return numVertexes;
}


template <typename T, typename U> vertex<T, U>*
                        graph<T, U>::addVertex(vertex<T, U> *newVertex){
  return addVertex(newVertex->value);
}


template <typename T, typename U> graph<T, U> graph<T, U>::operator=(
                                              const graph<T, U> &other){

  raise(SIGABRT);

  graph<T, U> toReturn;

  for(size_t i = 0; i < this->numVertexes; i++)
    toReturn.addVertex(*this->vertexArray[i]);

  for(size_t i = 0; i < this->numEdges; i++)
    toReturn.addEdge(*this->edgeArray[i]);

  return toReturn;
}


template <typename T, typename U> vertex<T, U>*
                    graph<T, U>::getVertexForValue(const T &testValue){
  if(geneNameToNodeID.count(testValue))
    return vertexArray[geneNameToNodeID[testValue]];
  return NULL;
}


template <typename T, typename U> void graph<T, U>::hintNumEdges(
                                                  csize_t suggestSize){
  void *memCheck;

  if(suggestSize <= numEdges) return;

  if(suggestSize == edgeArraySize) return;

  memCheck = realloc(edgeArray, suggestSize * sizeof(*edgeArray));
  if(NULL != memCheck){
    edgeArray = (edge<T, U>**) memCheck;
    edgeArraySize = suggestSize;
  }else{
    fprintf(stderr, "ERROR: Could not allocate edges\n"); fflush(stderr);
    raise(SIGABRT);
  }
}


template <typename T, typename U> void graph<T, U>::hintNumVertexes(
                                                  csize_t suggestSize){
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


template <typename T, typename U> void
                                graph<T, U>::shrinkEdgeCapacityToFit(){
  hintNumEdges(numEdges);
}


template <typename T, typename U> void
                              graph<T, U>::shrinkVertexCapacityToFit(){
  hintNumVertexes(numVertexes);
}


template <typename T, typename U> void graph<T, U>::shrinkToFit(){
  shrinkEdgeCapacityToFit();
  shrinkVertexCapacityToFit();
}


template <typename T, typename U> void graph<T, U>::ensureEdgeCapacity(
                                                          csize_t size){
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


template <typename T, typename U> void
                  graph<T, U>::ensureVertexCapacity(const size_t size){
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

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif
