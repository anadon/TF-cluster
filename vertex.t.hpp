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

#ifndef STD_VERTEX_T_HPP
#define STD_VERTEX_T_HPP

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <csignal>
#include <cstring>

#include "vertex.hpp"
#include "edge.hpp"

////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

template <typename T, typename U> vertex<T, U>::vertex(size_t index,
                              T data):  vertexIndex(index), value(data){
  numEdges = edgesSize = 0;
  edges = (edge<T, U>**) NULL;
}


template <typename T, typename U> inline size_t
                                      vertex<T, U>::getNumEdges() const{
  return numEdges;
}


template <typename T, typename U> vertex<T, U>::~vertex(){
  if(0 != numEdges) raise(SIGABRT);
}


template <typename T, typename U> size_t vertex<T, U>::addEdge(
                                                edge<T, U> *toRegister){
  ensureEdgeCapacity(numEdges+1);
  edges[numEdges] = toRegister;
  pair<vertex<T, U>*, edge<T, U>* > toInsert;
  toInsert = pair<vertex<T, U>*, edge<T, U>* >(toRegister->other(this), toRegister);
  connected.insert(toInsert);
  numEdges++;

  return numEdges-1;
}


template <typename T, typename U> void vertex<T, U>::removeEdge(
                                                  edge<T, U> *toRemove){
  void *memCheck;
  edge<T, U> *tmp;
  size_t targetEdgeIndex = 0;


  numEdges--;

  //check if this is the right or left of the edge, error if edge
  //doesn't connect this node/vertex
  if(toRemove->leftEdgeIndex <= numEdges
                        && toRemove == edges[toRemove->leftEdgeIndex]){
    targetEdgeIndex = toRemove->leftEdgeIndex;
  }else if(toRemove->rightEdgeIndex <= numEdges
                        && toRemove == edges[toRemove->rightEdgeIndex]){
    targetEdgeIndex = toRemove->rightEdgeIndex;
  }else{
    raise(SIGABRT);
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
  if(0 < numEdges){
    memCheck = realloc(edges, numEdges * sizeof(*edges));
    edges = (edge<T, U>**) memCheck;
  }else{
    free(edges);
    edges = NULL;
  }
  edgesSize = numEdges;
  
  connected.erase(toRemove->other(this));
}


template <typename T, typename U> inline bool vertex<T, U>::operator==(
                                      const vertex<T, U> &other) const{
  return (value == other.value);
}


template <typename T, typename U> edge<T, U>** vertex<T, U>::getEdges(){
  return edges;
}


template <typename T, typename U> const edge<T, U>**
                                        vertex<T, U>::getEdges() const{
  return (const edge<T, U>**) edges;
}


template <typename T, typename U> void vertex<T, U>::hintNumEdges(
                                              const size_t suggestSize){
  void *tmpPtr;

  if(suggestSize <= numEdges) return;

  tmpPtr = realloc(edges, sizeof(*edges) * suggestSize);
  if(NULL == tmpPtr) raise(SIGABRT);
  edges = (edge<T, U>**) tmpPtr;

  edgesSize = suggestSize;
}


template <typename T, typename U> void vertex<T, U>::shrinkToFit(){
  hintNumEdges(numEdges);
}


template <typename T, typename U> void vertex<T, U>::ensureEdgeCapacity(
                                                    const size_t size){
  while(size > edgesSize)
    hintNumEdges(1 + (edgesSize << 1));
}


template <typename T, typename U> bool vertex<T, U>::areConnected(
                                            vertex<T, U> *other) const{
  return connected.count(other);
}


////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif
