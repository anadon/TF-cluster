/***********************************************************************
         FILE:  edge.t.hpp

  DESCRIPTION:  Early representation of an edge.

         BUGS:  ---
        NOTES:  ---
       AUTHOR:  Josh Marshall <jrmarsha@mtu.edu>
      COMPANY:  Michigan technological University
      VERSION:  See git log
      CREATED:  See git log
     REVISION:  See git log
     LISCENSE:  GPLv3
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
  leftEdgeIndex = left->addEdge(this);
  rightEdgeIndex = right->addEdge(this);
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