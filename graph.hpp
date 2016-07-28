/***********************************************************************
         FILE:  graph.hpp

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
#ifndef STD_GRAPH_HPP
#define STD_GRAPH_HPP

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <unistd.h>
#include <unordered_map>

#include "edge.hpp"
#include "geneData.hpp"
#include "vertex.hpp"

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::unordered_map;

////////////////////////////////////////////////////////////////////////
//CLASS DEFINITION//////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  
 **********************************************************************/
template <typename T, typename U> class graph{
  private:
  vertex<T, U> **vertexArray;
  edge<T, U> **edgeArray;
  size_t numVertexes, numEdges;
  size_t vertexArraySize, edgeArraySize;
  unordered_map<T, size_t, std::hash<T>> geneNameToNodeID;

/*GRAPH OPERATIONS*****************************************************/
  public:
  
/*******************************************************************//**
 *  Basic constructor.
 **********************************************************************/
  graph();
  
  
/*******************************************************************//**
 *  Handles the deallocation of edges and vertexes (but not nesicarily
 * the user defined values they hold) before self deletion.
 **********************************************************************/
  ~graph();
  
  
/*******************************************************************//**
 *  Copy the contents of one graph to another.
 **********************************************************************/
  graph<T, U>& operator=(const graph<T, U> &other);
  
  
/*******************************************************************//**
 *  Minimize the memory used to hold vertexes and edges.
 **********************************************************************/
  void shrinkToFit();
  
/*EDGE OPERATIONS******************************************************/
  public:
  
/*******************************************************************//**
 *  Create and register a new edge.
 * 
 * @param[in,out] left A vertex to connect with an edge.
 * @param[in,out] right A vertex to connect with an edge.
 * @param[in] newWeight Value new edge will hold.
 **********************************************************************/
  edge<T, U>* addEdge(vertex<T, U> *left, vertex<T, U> *right, 
                                                          U newWeight);


/*******************************************************************//**
 *  Copy and register a new edge from another edge (possibly in another
 * graph).
 * 
 * @param[in] toAdd edge to replicate.
 * 
 * TODO -- can toAdd be constant?
 **********************************************************************/
  edge<T, U>* addEdge(edge<T, U> &toAdd);


/*******************************************************************//**
 *  Get array of edges in graph.
 **********************************************************************/
  edge<T, U>** getEdges();


/*******************************************************************//**
 *  Get array of edges in graph.
 **********************************************************************/
  const edge<T, U>** getEdges() const;


/*******************************************************************//**
 *  Get number of edges in graph.
 **********************************************************************/
  size_t getNumEdges() const;


/*******************************************************************//**
 *  Suggest number of edges gaph should be able to hold.  Used to help
 * memory allocation and management.
 * 
 * @param[in,out] suggestSize Suggested number of edges graph should 
 *                            have capacity for.
 **********************************************************************/
  void hintNumEdges(const size_t suggestSize);


/*******************************************************************//**
 *  Safely remove an edge from the graph.  All other methods of edge
 * removal will corrupt graph data.
 * 
 * @param[in,out] toRemove Edge to de-register from it's connected 
 *                         verticies and the graph.
 **********************************************************************/
  U removeEdge(edge<T, U> *toRemove);


/*******************************************************************//**
 *  Minimize the space needed to hold current number of edges.
 **********************************************************************/
  void shrinkEdgeCapacityToFit();
  
  private:
  
/*******************************************************************//**
 *  Ensure that the graph has capacity for size number of edges.
 * 
 * @param[in] size Number of edges the graph should have capacity to 
 *                 hold.
 **********************************************************************/
  void ensureEdgeCapacity(const size_t size);
  
/*VERTEX OPERATIONS****************************************************/
  public:
  
/*******************************************************************//**
 *  Construct and register vertex in graph.
 * 
 * @param[in] newValue Value the new vertex will hold.
 **********************************************************************/
  vertex<T, U>* addVertex(T newValue);


/*******************************************************************//**
 *  Copy and register vertex in graph.
 * 
 * @param[in,out] newVertex Vertex to copy contents of, but not edges.
 **********************************************************************/
  vertex<T, U>* addVertex(vertex<T, U> *newVertex);


/*******************************************************************//**
 *  Get the edges which hold the contents equal to testValue.
 * 
 * @param[in] testValue Data which a vertex in the called graph should
 *                      hold.
 **********************************************************************/
  vertex<T, U>* getVertexForValue(const T &testValue);


/*******************************************************************//**
 *  Get array of vertexes in graph;
 **********************************************************************/
  vertex<T, U>** getVertexes();


/*******************************************************************//**
 *  Get number of vertexes in graph;
 **********************************************************************/
  size_t getNumVertexes() const;


/*******************************************************************//**
 *  Suggest number of vertexes graph should have capacity to hold.
 * 
 * @param[in] suggestSize Number of vertexes graph should have capacity 
 *                        to hold.
 **********************************************************************/
  void hintNumVertexes(const size_t suggestSize);


/*******************************************************************//**
 *  Remove the registered vertex from graph after removing all of its
 * edges, and return the data to held.
 * 
 * @param[in,out] toRemove Vertex in called graph which should be 
 *                         removed.
 **********************************************************************/
  T removeVertex(const vertex<T, U> *toRemove);
  
  
/*******************************************************************//**
 *  Minimize the memory used to hold vertexes in called graph.
 **********************************************************************/
  void shrinkVertexCapacityToFit();
  
  private:

/*******************************************************************//**
 *  Ensure the graph has capacity for size vertexes.
 * 
 * @param[in] size Number of vertexes graph should be able to hold.
 **********************************************************************/
  void ensureVertexCapacity(const size_t size);
  
};

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#endif