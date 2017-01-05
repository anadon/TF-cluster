/*******************************************************************//**
         FILE:  vertex.hpp

  DESCRIPTION:  Public interface for a somewhats STL quality graph

         BUGS:  ---
        NOTES:  ---
       AUTHOR:  Josh Marshall <jrmarsha@mtu.edu>
      COMPANY:  Michigan technological University
      VERSION:  See git log
      CREATED:  See git log
     REVISION:  See git log
     LISCENSE:  GPLv3
***********************************************************************/
#ifndef STD_VERTEX_HPP
#define STD_VERTEX_HPP

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <unistd.h>
#include <string>

////////////////////////////////////////////////////////////////////////
//TYPE FORWARD DECLARATIONS/////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

template <typename T, typename U> class graph;

template <typename T, typename U> class edge;

////////////////////////////////////////////////////////////////////////
//CLASS DECLARATION/////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  A moderately complete and well formed vertex
 **********************************************************************/
template <typename T, typename U> class vertex{
  private:
  edge<T, U> **edges;
  size_t numEdges, edgesSize;
  unordered_map<vertex<T, U>*, edge<T, U>*> connected;

  public:
  size_t vertexIndex;
  T value;

/*******************************************************************//**
 *  Make a vertex with basic data about itself
 *
 * @param[in] index Index the vertex will be located at in graph
 * @param[in] data Data the vertex is meant to represent
 **********************************************************************/
  vertex(size_t index, T data);


/*******************************************************************//**
 *  Basic destructor.  All edges must have been removed by graph prior
 * to deletion.
 *
 * \pre All edges must have been removed already.
 **********************************************************************/
  ~vertex();


/*******************************************************************//**
 *  Register an edge on called vertex.
 *
 * @param[in,out] toRegister Edge to register and register with by
                             returning to it what index it is in the
                             called vertex.
 **********************************************************************/
  size_t addEdge(edge<T, U> *toRegister);


/*******************************************************************//**
 *  Get array of stored edges.
 **********************************************************************/
  edge<T, U>** getEdges();


/*******************************************************************//**
 *  Get array of stored edges.
 **********************************************************************/
  const edge<T, U>** getEdges() const;


/*******************************************************************//**
 *  Get number of stored edges.
 **********************************************************************/
  size_t getNumEdges() const;


/*******************************************************************//**
 *  Given an edge, remove references from it in called vertex
 *
 * @param[in,out] toRemove Edge to de-register in called vertex.
 **********************************************************************/
  void removeEdge(edge<T, U> *toRemove);


/*******************************************************************//**
 *  Suggest number of edges to be able to store.  Use this to optimize
 * memory management.
 *
 * @param[in] suggestSize Suggested number of edges to accomidate.  Can
                          be used to increase or decrease allocated
                          size.
 **********************************************************************/
  void hintNumEdges(const size_t suggestSize);


/*******************************************************************//**
 *  Minimize allocated space to fit current number of edges.
 **********************************************************************/
  void shrinkToFit();


/*******************************************************************//**
 *  Tell if contents of vertexes are the same, but not nessicarily the
 * same vertex from a single graph.
 **********************************************************************/
  bool operator==(const vertex<T, U> &other) const;
  
  
/*******************************************************************//**
 * Tell if there is an edge connecting this vertex to another vertex
 **********************************************************************/
  bool areConnected(vertex<T, U> *other) const;


  private:

/*******************************************************************//**
 *  Make sure than when adding an edge there is enough space.
 *
 * @param[in] size make sure vertex can accomidate at least size edges.
 **********************************************************************/
  void ensureEdgeCapacity(const size_t size);

};



//namespace std{
  //template <> struct hash< *vertex > {

///*******************************************************************//**
 //*  Add a hash specialization to handle geneData in the standard
 //* namespace.  Hash is given as the name index, since it can be assumed
 //* this is unique within a given graph.
 //*
 //* @param[in] toHash geneData that will give it's nameIndex as a hash.
 //**********************************************************************/
    //std::size_t operator()(vertex const& *toHash) const {
      //return (size_t) toHash;
    //}
  //};
//}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif
