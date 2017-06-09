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
