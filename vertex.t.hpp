/*******************************************************************//**
         FILE:  vertex.t.hpp

  DESCRIPTION:  Implementation for a somewhats STL quality vertex

         BUGS:  ---
        NOTES:  ---
       AUTHOR:  Josh Marshall <jrmarsha@mtu.edu>
      COMPANY:  Michigan technological University
      VERSION:  See git log
      CREATED:  See git log
     REVISION:  See git log
     LISCENSE:  GPLv3
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
