/***********************************************************************
         FILE:  edge.hpp

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
#ifndef EDGE_HPP
#define EDGE_HPP

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <unistd.h>

////////////////////////////////////////////////////////////////////////
//FORWARDS CLASS DECLARATION////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

template <typename T, typename U> class vertex;

////////////////////////////////////////////////////////////////////////
//CLASS DEFINITION//////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  Moderately complete edge implementation.
 **********************************************************************/
template <typename T, typename U> class edge{
  public:
  vertex<T, U> *left, *right;
  U weight;
  size_t edgeID, leftEdgeIndex, rightEdgeIndex;


/*******************************************************************//**
 *  Make a new edge between two verticies with a given weight.
 * 
 * @param[in,out] newLeft Left vertex to connect.
 * @param[in,out] newRight Right vertex to connect.
 * @param[in] newWeight Weight value.
 * @param[in] edgeIndex Index of this edge in the containing graph.
 **********************************************************************/
  edge(vertex<T, U> *newLeft, vertex<T, U> *newRight, U newWeight, 
                                                const size_t edgeIndex);


/*******************************************************************//**
 *  Basic deconstructor.
 **********************************************************************/
  ~edge();
  
  
/*******************************************************************//**
 *  Given a connected vertex, tell what the other connected vertex is.
 **********************************************************************/
  vertex<T, U>* other(const vertex<T, U> *side);
};

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif