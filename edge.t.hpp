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

#ifndef EDGE_T_HPP
#define EDGE_T_HPP
#include "edge.hpp"
#include "vertex.hpp"

#include <csignal>


template <typename T, typename U> edge<T, U>::edge(vertex<T, U> *newLeft,
                vertex<T, U> *newRight, U newWeight, size_t edgeIndex):
                left(newLeft), right(newRight), weight(newWeight),
                                                      edgeID(edgeIndex){
  leftEdgeIndex = left->addLeftEdge(this);
  rightEdgeIndex = right->addRightEdge(this);
}


template <typename T, typename U> edge<T, U>::~edge(){
  left->removeEdge(this);
  right->removeEdge(this);
}


template <typename T, typename U> vertex<T, U>* edge<T, U>::other(const vertex<T, U> *side){
  if(side == left)        return right;
  else if(side == right)  return left;
  else                    raise(SIGABRT);
  return ( vertex<T, U>* ) NULL;
}

#endif
