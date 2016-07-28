/*******************************************************************//**
         FILE:  tripleLink.cpp

  DESCRIPTION:  Implementation of tripleLink clustering

         BUGS:  ---
        NOTES:  ---
       AUTHOR:  Josh Marshall <jrmarsha@mtu.edu>
      COMPANY:  Michigan technological University
      VERSION:  See git log
      CREATED:  See git log
     REVISION:  See git log
     LISCENSE:  GPLv3
***********************************************************************/

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <queue>
#include <vector>

#include "auxillaryUtilities.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"
#include "graph.t.hpp"

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::queue;
using std::vector;

////////////////////////////////////////////////////////////////////////
//PRIVATE FUNCTION DECLARATIONS/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  Expand a cluster using the triple-link heuristic, and return a
 * queue of name indexes found, while removing those vertexes the 
 * indexes came from from the graph.
 * 
 * @param[in,out] geneNetwork Graph of interconnected genes
 * @param[in] threeSigma High connection value for edges.
 * @param[in] twoSigma Medium connection value for edges.
 **********************************************************************/
queue<size_t> tripleLinkIteration(graph<geneData, double> *geneNetwork,
                                        cf64 threeSigma, cf64 twoSigma);


/*******************************************************************//**
 *  Mark vertex as reached as appropriate and return true if edge can be
 * safely removed.
 * 
 * @param[in] edgeWeight The weight of the edge connecting the outgoing 
 *                       vertex.
 * @param[in,out] toMark The outgoing vertex to mark.
 * @param[in] high High value edge weight cutoff.
 * @param[in] med medium value edge weight cutoff.
 * @param[out] toProcess Processing queue should toMark becomes well
 *                       connected.
 **********************************************************************/
bool markConnectedVertex(cf64 edgeWeight, 
                  vertex<geneData, double> *toMark, cf64 high, cf64 med, 
                                            queue<geneData> &toProcess);


/*******************************************************************//**
 *  Appropriately mark all vertexes connected to markFrom, and add ones
 * who are well connected to the current process queue.  Also remove
 * traversed edges on markFrom, as they are only needed once.
 * 
 * @param[in,out] markFrom Vertex to mark all connecting vertexes from.
 * @param[in] high High value edge weight cutoff.
 * @param[in] med Medium value edge weight cutoff.
 * @param[in,out] geneNetwork Network of connected genes which markFrom
 *                            needs to operate through and be 
 *                            subsequently removed.
 * @param[in,out] toProcess Record of vertexes which need to be 
 *                          processed in this iteration of triple-link.
 **********************************************************************/
void markConnectedVertexes(vertex<geneData, double> *markFrom, cf64 high, 
                        cf64 med, graph<geneData, double> *geneNetwork, 
                                            queue<geneData> &toProcess);


/*******************************************************************//**
 *  Reset connection markers on verticex.
 * 
 * @param[out] toReset Remove connection marks from passed vertex.
 **********************************************************************/
inline void untouchVertex(vertex<geneData, double> *toReset);


/*******************************************************************//**
 *  Remove verticies which are apparent that they can no longer be 
 * included in any future cluster.
 * 
 * @param[in,out] geneNetwork Graph to search through and prune.
 * @param[in] high High value edge weight cutoff.
 * @param[in] med medium value edge weight cutoff.
 **********************************************************************/
void removeWeakVerticies(graph<geneData, f64> *geneNetwork, cf64 high, 
                                                              cf64 med);

////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

bool markConnectedVertex(cf64 edgeWeight, 
                  vertex<geneData, double> *toMark, cf64 high, cf64 med, 
                                            queue<geneData> &toProcess){
  geneData *value = &toMark->value;
  
  if(value->threeSigmaLink & value->twoSigmaLink & value->oneSigmaLink)
    return true;
  
  if(!value->threeSigmaLink && edgeWeight >= high){    
                                        value->threeSigmaLink = true; }
  else if(!value->twoSigmaLink && edgeWeight >= med){  
                                        value->twoSigmaLink   = true; }
  else if(!value->oneSigmaLink){                       
                                        value->oneSigmaLink   = true; }
  
  if(value->threeSigmaLink & value->twoSigmaLink & value->oneSigmaLink){
    toProcess.push(toMark->value.nameIndex);
    return true;
  }else{
    return false;
  }
}


void markConnectedVertexes(vertex<geneData, double> *markFrom, 
              cf64 high, cf64 med, graph<geneData, double> *geneNetwork, 
                                            queue<geneData> &toProcess){
  for(size_t i = 0; i < markFrom->getNumEdges(); i++)
    if(markConnectedVertex(markFrom->getEdges()[i]->weight,
            markFrom->getEdges()[i]->other(markFrom), high, med, 
                                                            toProcess))
      geneNetwork->removeEdge(markFrom->getEdges()[i]);
}


inline void untouchVertex(vertex<geneData, double> *toReset){
  toReset->value.threeSigmaLink = false;
  toReset->value.twoSigmaLink = false;
  toReset->value.oneSigmaLink = false;
}


//NOTE: This can be optimized
void removeWeakVerticies(graph<geneData, f64> *geneNetwork, cf64 high, 
                                                              cf64 med){
  //first, we need to remove all nodes which do not have 3 available
  //links.  We make the assumption that all edges below tripleLink3 are
  //removed.
  bool disconnectedVerticiesFound;
  do{
    disconnectedVerticiesFound = false;
    for(size_t i = 0; i < geneNetwork->getNumVertexes(); i++){
      const vertex<geneData, f64> *target = 
                                          geneNetwork->getVertexes()[i];
      if(2 > target->getNumEdges()){
        geneNetwork->removeVertex(target);
        disconnectedVerticiesFound = true;
        break;
      }
      
      bool highFound, medFound;
      highFound = medFound = false;
      for(size_t j = 0; j < target->getNumEdges() 
                                    && (!highFound || !medFound); j++){
        if(target->getEdges()[j]->weight >= high && !highFound)
          highFound = true;
        else if(target->getEdges()[j]->weight >= med)
          medFound = true;
      }
      if(!highFound || !medFound){
        geneNetwork->removeVertex(target);
        disconnectedVerticiesFound = true;
        break;
      }
    }
  }while(disconnectedVerticiesFound);
}


queue<size_t> tripleLinkIteration(graph<geneData, double> *geneNetwork,
                                        cf64 threeSigma, cf64 twoSigma){
  queue<size_t> toReturn;
  edge<geneData, double> *initialEdge;
  vertex<geneData, double> *firstVertex, *secondVertex;
  vertex<geneData, double> *connectedVertex;
  queue<geneData> toProcess;
  size_t targetEdgeIndex = 0;
  double maxFoundValue = geneNetwork->getEdges()[0]->weight;

  //Find strongest edge, use this to grow the tree
  for(size_t i = 1; i < geneNetwork->getNumEdges(); i++){
    f64 weight = geneNetwork->getEdges()[i]->weight;
    if(weight >= maxFoundValue){
      maxFoundValue = weight;
      targetEdgeIndex = i;
    }
  }
  
  //Reset all verticies to untouched
  for(size_t i = 0; i < geneNetwork->getNumVertexes(); i++)
    untouchVertex(geneNetwork->getVertexes()[i]);

  //Add the highest weighted edge's verticies to processing queue
  initialEdge = geneNetwork->getEdges()[targetEdgeIndex];
  firstVertex = initialEdge->left;
  secondVertex = initialEdge->right;
  geneNetwork->removeEdge(geneNetwork->getEdges()[targetEdgeIndex]);

  toReturn.push(firstVertex->value.nameIndex);
  toReturn.push(secondVertex->value.nameIndex);

  //Primer connections for triple link
  markConnectedVertexes(firstVertex, threeSigma, twoSigma, geneNetwork, 
                                                            toProcess);
  geneNetwork->removeVertex(firstVertex);
  
  markConnectedVertexes(secondVertex, threeSigma, twoSigma, geneNetwork, 
                                                            toProcess);
  geneNetwork->removeVertex(secondVertex);

  //Main tiple link tree expantion loop
  while(!toProcess.empty()){
    connectedVertex = geneNetwork->getVertexForValue(toProcess.front());
    toProcess.pop();
    if(NULL == connectedVertex){
    //  printf("ignoring redundantly included vertex\n");
      continue;
    }
    
    toReturn.push(connectedVertex->value.nameIndex);

    markConnectedVertexes(connectedVertex, threeSigma, twoSigma, 
                                                geneNetwork, toProcess);
    geneNetwork->removeVertex(connectedVertex);
  }

  return toReturn;
}


queue< queue<size_t> > tripleLink(graph<geneData, double> *geneNetwork,
                                        cf64 threeSigma, cf64 twoSigma){
  queue< queue<size_t> > toReturn;

  //remove verticies which lack connection strength for inclusion
  removeWeakVerticies(geneNetwork, threeSigma, twoSigma);

  while(geneNetwork->getNumVertexes() > 0 
                                    && geneNetwork->getNumEdges() > 0){
    toReturn.push(tripleLinkIteration(geneNetwork, threeSigma, 
                                                            twoSigma));
    removeWeakVerticies(geneNetwork, threeSigma, twoSigma);
  }

  return toReturn;
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
