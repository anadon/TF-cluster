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
  numLeftEdges = numRightEdges = leftEdgesSize = rightEdgesSize = 0;
  leftEdges = rightEdges = (edge<T, U>**) NULL;
}


template <typename T, typename U> inline size_t
                                      vertex<T, U>::getNumLeftEdges() const{
  return numLeftEdges;
}


template <typename T, typename U> inline size_t
                                      vertex<T, U>::getNumRightEdges() const{
  return numRightEdges;
}


template <typename T, typename U> vertex<T, U>::~vertex(){
  if(0 != numLeftEdges) raise(SIGABRT);
  if(0 != numRightEdges) raise(SIGABRT);
}


template <typename T, typename U> size_t vertex<T, U>::addLeftEdge(
                                                edge<T, U> *toRegister){
  ensureEdgeCapacity(numLeftEdges+1);
  leftEdges[numLeftEdges] = toRegister;
  pair<vertex<T, U>*, edge<T, U>* > toInsert;
  toInsert = pair<vertex<T, U>*, edge<T, U>* >(toRegister->other(this), toRegister);
  connected.insert(toInsert);
  numLeftEdges++;

  return numLeftEdges-1;
}


template <typename T, typename U> size_t vertex<T, U>::addRightEdge(
                                                edge<T, U> *toRegister){
  ensureEdgeCapacity(numRightEdges+1);
  rightEdges[numRightEdges] = toRegister;
  pair<vertex<T, U>*, edge<T, U>* > toInsert;
  toInsert = pair<vertex<T, U>*, edge<T, U>* >(toRegister->other(this), toRegister);
  connected.insert(toInsert);
  numRightEdges++;

  return numRightEdges-1;
}


template <typename T, typename U> void vertex<T, U>::removeEdge(
                                                  edge<T, U> *toRemove){
  void *memCheck;
  edge<T, U> *tmp;
  size_t targetEdgeIndex = 0;
  bool isleft;

  if(this == toRemove->left){
    int found = 0;
    for(size_t i = 0; i < numLeftEdges; i++){
      if(leftEdges[i] == toRemove){
        found++;
      }
    }
    if(found == 0) raise(SIGABRT);
    if(found != 1) raise(SIGABRT);
    isleft = true;
  }else if(this == toRemove->right){
    int found = 0;
    for(size_t i = 0; i < numRightEdges; i++)
      if(rightEdges[i] == toRemove) found++;
    if(0 == found) raise(SIGABRT);
    if(1 != found) raise(SIGABRT);
    isleft = false;
  }else raise(SIGABRT);



  //check if this is the right or left of the edge, error if edge
  //doesn't connect this node/vertex
  if(isleft){
    if(toRemove->leftEdgeIndex >= numLeftEdges)
      raise(SIGABRT);
    if(toRemove != leftEdges[toRemove->leftEdgeIndex]){
      raise(SIGABRT);
    }
    targetEdgeIndex = toRemove->leftEdgeIndex;
    isleft = true;
  }else{
    if(toRemove->rightEdgeIndex >= numRightEdges)
      raise(SIGABRT);
    if(toRemove != rightEdges[toRemove->rightEdgeIndex]){
      raise(SIGABRT);
    }
    targetEdgeIndex = toRemove->rightEdgeIndex;
    isleft = false;
  }

  if(isleft){
    numLeftEdges--;

    //swap edge this is removing to the end
    tmp = leftEdges[numLeftEdges];
    leftEdges[numLeftEdges] = leftEdges[targetEdgeIndex];
    leftEdges[targetEdgeIndex] = tmp;
  }else{
    numRightEdges--;

    //swap edge this is removing to the end
    tmp = rightEdges[numRightEdges];
    rightEdges[numRightEdges] = rightEdges[targetEdgeIndex];
    rightEdges[targetEdgeIndex] = tmp;
  }

  //update the swapped edge's location in this structure so it still
  //knows where it is in this vertex/node
  if(isleft && this != leftEdges[targetEdgeIndex]->left)
    raise(SIGABRT);
  if(!isleft && this != rightEdges[targetEdgeIndex]->right)
    raise(SIGABRT);

  if(isleft)
    leftEdges[targetEdgeIndex]->leftEdgeIndex = targetEdgeIndex;
  else
    rightEdges[targetEdgeIndex]->rightEdgeIndex = targetEdgeIndex;

  //deallocate the last edge pointer (but don't actually delete --
  //that's the network's job).
  if(isleft){
    if(0 < numLeftEdges){
      memCheck = realloc(leftEdges, numLeftEdges * sizeof(*leftEdges));
      leftEdges = (edge<T, U>**) memCheck;
    }else{
      free(leftEdges);
      leftEdges = NULL;
    }
    leftEdgesSize = numLeftEdges;
  }else{
    if(0 < numRightEdges){
      memCheck = realloc(rightEdges, numRightEdges * sizeof(*rightEdges));
      rightEdges = (edge<T, U>**) memCheck;
    }else{
      free(rightEdges);
      rightEdges = NULL;
    }
    rightEdgesSize = numRightEdges;
  }

  connected.erase(toRemove->other(this));
}


template <typename T, typename U> inline bool vertex<T, U>::operator==(
                                      const vertex<T, U> &other) const{
  return (value == other.value);
}


template <typename T, typename U> edge<T, U>** vertex<T, U>::getLeftEdges(){
  return leftEdges;
}


template <typename T, typename U> edge<T, U>** vertex<T, U>::getRightEdges(){
  return rightEdges;
}


template <typename T, typename U> const edge<T, U>**
                                        vertex<T, U>::getLeftEdges() const{
  return (const edge<T, U>**) leftEdges;
}


template <typename T, typename U> const edge<T, U>**
                                        vertex<T, U>::getRightEdges() const{
  return (const edge<T, U>**) rightEdges;
}


template <typename T, typename U> void vertex<T, U>::hintNumEdges(
                                              const size_t suggestSize){
  void *tmpPtr;

  if(suggestSize <= numLeftEdges) return;
  if(suggestSize <= numRightEdges) return;

  tmpPtr = realloc(leftEdges, sizeof(*leftEdges) * suggestSize);
  if(NULL == tmpPtr) raise(SIGABRT);
  leftEdges = (edge<T, U>**) tmpPtr;

  tmpPtr = realloc(rightEdges, sizeof(*rightEdges) * suggestSize);
  if(NULL == tmpPtr) raise(SIGABRT);
  rightEdges = (edge<T, U>**) tmpPtr;

  leftEdgesSize = suggestSize;
  rightEdgesSize = suggestSize;
}


template <typename T, typename U> void vertex<T, U>::shrinkToFit(){
  hintNumEdges(numLeftEdges);
  hintNumEdges(numRightEdges);
}


template <typename T, typename U> void vertex<T, U>::ensureEdgeCapacity(
                                                    const size_t size){
  while(size > leftEdgesSize)
    hintNumEdges(1 + (leftEdgesSize << 1));
  while(size > rightEdgesSize)
    hintNumEdges(1 + (rightEdgesSize << 1));
}


template <typename T, typename U> bool vertex<T, U>::areConnected(
                                            vertex<T, U> *other) const{
  return connected.count(other);
}


////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif
