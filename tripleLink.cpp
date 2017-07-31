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

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <queue>
#include <vector>

#include "auxiliaryUtilities.hpp"
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
queue<size_t> tripleLinkIteration(graph<geneData, u8> *geneNetwork,
                                        cu8 threeSigma, cu8 twoSigma);


/*******************************************************************//**
 *  Mark vertex as reached as appropriate and return true if edge can be
 * safely removed.  Requires only a single connection to be considered
 * reached.
 *
 * @param[in] edgeWeight The weight of the edge connecting the outgoing
 *                       vertex.
 * @param[in,out] toMark The outgoing vertex to mark.
 * @param[in] high High value edge weight cutoff.
 * @param[in] med medium value edge weight cutoff.
 * @param[out] toProcess Processing queue should toMark becomes well
 *                       connected.
 **********************************************************************/
bool markConnectedVertexSingle(cu8 edgeWeight,
                  vertex<geneData, u8> *toMark, cu8 high, cu8 med,
                                            queue<geneData> &toProcess);


/*******************************************************************//**
 *  Mark vertex as reached as appropriate and return true if edge can be
 * safely removed.  Requires two stronger connection to be considered
 * reached.
 *
 * @param[in] edgeWeight The weight of the edge connecting the outgoing
 *                       vertex.
 * @param[in,out] toMark The outgoing vertex to mark.
 * @param[in] high High value edge weight cutoff.
 * @param[in] med medium value edge weight cutoff.
 * @param[out] toProcess Processing queue should toMark becomes well
 *                       connected.
 **********************************************************************/
bool markConnectedVertexDouble(cu8 edgeWeight,
                  vertex<geneData, u8> *toMark, cu8 high, cu8 med,
                                            queue<geneData> &toProcess);


/*******************************************************************//**
 *  Mark vertex as reached as appropriate and return true if edge can be
 * safely removed.  Requires all three connections to be considered
 * reached.
 *
 * @param[in] edgeWeight The weight of the edge connecting the outgoing
 *                       vertex.
 * @param[in,out] toMark The outgoing vertex to mark.
 * @param[in] high High value edge weight cutoff.
 * @param[in] med medium value edge weight cutoff.
 * @param[out] toProcess Processing queue should toMark becomes well
 *                       connected.
 **********************************************************************/
bool markConnectedVertexTriple(cu8 edgeWeight,
                  vertex<geneData, u8> *toMark, cu8 high, cu8 med,
                                            queue<geneData> &toProcess);


/*******************************************************************//**
 *  Appropriately mark all vertexes connected to markFrom, and add ones
 * who are well connected to the current process queue.  Also remove
 * traversed edges on markFrom, as they are only needed once.  Has a
 * requirement of 1 edge of maximum strength for inclusion into the
 * current cluster.
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
void markConnectedVertexesSingle(vertex<geneData, u8> *markFrom,
              cu8 high, cu8 med, graph<geneData, u8> *geneNetwork,
                                          queue<geneData> &toProcessTo);


/*******************************************************************//**
 *  Appropriately mark all vertexes connected to markFrom, and add ones
 * who are well connected to the current process queue.  Also remove
 * traversed edges on markFrom, as they are only needed once.  Has a
 * requirement of top 2 strong edges for inclusion into the current
 * cluster.
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
void markConnectedVertexesDouble(vertex<geneData, u8> *markFrom,
              cu8 high, cu8 med, graph<geneData, u8> *geneNetwork,
                                          queue<geneData> &toProcessTo);


/*******************************************************************//**
 *  Appropriately mark all vertexes connected to markFrom, and add ones
 * who are well connected to the current process queue.  Also remove
 * traversed edges on markFrom, as they are only needed once.  Has a
 * requirement of all 3 edges for inclusion into the current cluster.
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
void markConnectedVertexesTriple(vertex<geneData, u8> *markFrom,
              cu8 high, cu8 med, graph<geneData, u8> *geneNetwork,
                                          queue<geneData> &toProcessTo);


/*******************************************************************//**
 *  Reset connection markers on verticex.
 *
 * @param[out] toReset Remove connection marks from passed vertex.
 **********************************************************************/
inline void untouchVertex(vertex<geneData, u8> *toReset);


/*******************************************************************//**
 *  Remove verticies which are apparent that they can no longer be
 * included in any future cluster.
 *
 * @param[in,out] geneNetwork Graph to search through and prune.
 * @param[in] high High value edge weight cutoff.
 * @param[in] med medium value edge weight cutoff.
 **********************************************************************/
void removeWeakVerticies(graph<geneData, u8> *geneNetwork, cu8 high,
                                                              cu8 med);

////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

bool markConnectedVertexSingle(cu8 edgeWeight,
                  vertex<geneData, u8> *toMark, cu8 high, cu8 med,
                                            queue<geneData> &toProcess){
  geneData *value = &toMark->value;

  if(value->threeSigmaLink)
    return true;

  if(!value->threeSigmaLink && edgeWeight >= high){
                                        value->threeSigmaLink = true; }
  else if(!value->twoSigmaLink && edgeWeight >= med){
                                        value->twoSigmaLink   = true; }
  else if(!value->oneSigmaLink){
                                        value->oneSigmaLink   = true; }

  if(value->threeSigmaLink){
    toProcess.push(toMark->value.nameIndex);
    return true;
  }else{
    return false;
  }
}


bool markConnectedVertexDouble(cu8 edgeWeight,
                  vertex<geneData, u8> *toMark, cu8 high, cu8 med,
                                            queue<geneData> &toProcess){
  geneData *value = &toMark->value;

  if(value->threeSigmaLink & value->twoSigmaLink)
    return true;

  if(!value->threeSigmaLink && edgeWeight >= high){
                                        value->threeSigmaLink = true; }
  else if(!value->twoSigmaLink && edgeWeight >= med){
                                        value->twoSigmaLink   = true; }
  else if(!value->oneSigmaLink){
                                        value->oneSigmaLink   = true; }

  if(value->threeSigmaLink & value->twoSigmaLink){
    toProcess.push(toMark->value.nameIndex);
    return true;
  }else{
    return false;
  }
}


bool markConnectedVertexTriple(cu8 edgeWeight,
                  vertex<geneData, u8> *toMark, cu8 high, cu8 med,
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


//TODO merge these so that they take a function pointer --> collapse
//into a single function.
void markConnectedVertexesSingle(vertex<geneData, u8> *markFrom,
              cu8 high, cu8 med, graph<geneData, u8> *geneNetwork,
                                          queue<geneData> &toProcessTo){
  for(size_t i = 0; i < markFrom->getNumLeftEdges(); i++)
    if(markConnectedVertexSingle(markFrom->getLeftEdges()[i]->weight,
            markFrom->getLeftEdges()[i]->other(markFrom), high, med, toProcessTo)){
      geneNetwork->removeEdge(markFrom->getLeftEdges()[i]);
      i = 0;
    }

  for(size_t i = 0; i < markFrom->getNumRightEdges(); i++)
    if(markConnectedVertexSingle(markFrom->getRightEdges()[i]->weight,
            markFrom->getRightEdges()[i]->other(markFrom), high, med, toProcessTo)){
      geneNetwork->removeEdge(markFrom->getRightEdges()[i]);
      i = 0;
    }
}


void markConnectedVertexesDouble(vertex<geneData, u8> *markFrom,
              cu8 high, cu8 med, graph<geneData, u8> *geneNetwork,
                                          queue<geneData> &toProcessTo){
  for(size_t i = 0; i < markFrom->getNumLeftEdges(); i++)
    if(markConnectedVertexDouble(markFrom->getLeftEdges()[i]->weight,
            markFrom->getLeftEdges()[i]->other(markFrom), high, med, toProcessTo)){
      geneNetwork->removeEdge(markFrom->getLeftEdges()[i]);
      i = 0;
    }

  for(size_t i = 0; i < markFrom->getNumRightEdges(); i++)
    if(markConnectedVertexDouble(markFrom->getRightEdges()[i]->weight,
      markFrom->getRightEdges()[i]->other(markFrom), high, med, toProcessTo)){
        geneNetwork->removeEdge(markFrom->getRightEdges()[i]);
        i = 0;
    }
}


void markConnectedVertexesTriple(vertex<geneData, u8> *markFrom,
              cu8 high, cu8 med, graph<geneData, u8> *geneNetwork,
                                          queue<geneData> &toProcessTo){
  for(size_t i = 0; i < markFrom->getNumLeftEdges(); i++)
    if(markConnectedVertexTriple(markFrom->getLeftEdges()[i]->weight,
            markFrom->getLeftEdges()[i]->other(markFrom), high, med,
                                                          toProcessTo)){
      geneNetwork->removeEdge(markFrom->getLeftEdges()[i]);
      i = 0;
    }

  for(size_t i = 0; i < markFrom->getNumRightEdges(); i++)
    if(markConnectedVertexTriple(markFrom->getRightEdges()[i]->weight,
          markFrom->getRightEdges()[i]->other(markFrom), high, med, toProcessTo)){
      geneNetwork->removeEdge(markFrom->getRightEdges()[i]);
      i = 0;
    }
}


inline void untouchVertex(vertex<geneData, u8> *toReset){
  toReset->value.threeSigmaLink = false;
  toReset->value.twoSigmaLink = false;
  toReset->value.oneSigmaLink = false;
}


void removeWeakVerticies(graph<geneData, u8> *geneNetwork, cu8 high,
                                                              cu8 med){
  //first, we need to remove all nodes which do not have 3 available
  //links.  We make the assumption that all edges below tripleLink3 are
  //removed.
  bool disconnectedVerticiesFound;
  do{
    disconnectedVerticiesFound = false;
    for(size_t i = 0; i < geneNetwork->getNumVertexes(); i++){
      const vertex<geneData, u8> *target =
                                          geneNetwork->getVertexes()[i];
      if(2 > (target->getNumLeftEdges() + target->getNumRightEdges())){
        geneNetwork->removeVertex(target);
        disconnectedVerticiesFound = true;
        continue;
      }

      bool highFound, medFound;
      highFound = medFound = false;
      for(size_t j = 0; j < target->getNumLeftEdges()
                                    && (!highFound || !medFound); j++){
        if(target->getLeftEdges()[j]->weight >= high && !highFound)
          highFound = true;
        else if(target->getLeftEdges()[j]->weight >= med)
          medFound = true;
      }
      for(size_t j = 0; j < target->getNumRightEdges()
                                    && (!highFound || !medFound); j++){
        if(target->getRightEdges()[j]->weight >= high && !highFound)
          highFound = true;
        else if(target->getRightEdges()[j]->weight >= med)
          medFound = true;
      }
      if(!highFound || !medFound){
        geneNetwork->removeVertex(target);
        disconnectedVerticiesFound = true;
      }
    }
  }while(disconnectedVerticiesFound);
}


queue<size_t> tripleLinkIteration(graph<geneData, u8> *geneNetwork,
                                        cu8 threeSigma, cu8 twoSigma){
  queue<size_t> toReturn;
  edge<geneData, u8> *initialEdge;
  vertex<geneData, u8> *firstVertex, *secondVertex;
  vertex<geneData, u8> *connectedVertex;
  queue<geneData> toProcessPrimer, toProcessMain;
  size_t targetEdgeIndex = 0;

  double maxFoundValue = geneNetwork->getEdges()[0]->weight;

  //Find strongest edge, use this to grow the tree
  for(size_t i = 1; i < geneNetwork->getNumEdges(); i++){
    u8 weight = geneNetwork->getEdges()[i]->weight;
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
  geneNetwork->removeEdge(initialEdge);


  //Primer connections (single link phase) for triple link
  markConnectedVertexesSingle(firstVertex, threeSigma, twoSigma,
                                          geneNetwork, toProcessPrimer);
  toReturn.push(firstVertex->value.nameIndex);
  geneNetwork->removeVertex(firstVertex);

  markConnectedVertexesSingle(secondVertex, threeSigma, twoSigma,
                                          geneNetwork, toProcessPrimer);
  toReturn.push(secondVertex->value.nameIndex);
  geneNetwork->removeVertex(secondVertex);

  //Double Link phase

  while(!toProcessPrimer.empty()){
    connectedVertex = geneNetwork->getVertexForValue(
                                              toProcessPrimer.front());
    toProcessPrimer.pop();

    if(NULL == connectedVertex) continue;

    toReturn.push(connectedVertex->value.nameIndex);
    markConnectedVertexesDouble(connectedVertex, threeSigma, twoSigma,
                                            geneNetwork, toProcessMain);
    geneNetwork->removeVertex(connectedVertex);
  }

  //Main triple link phase tree expantion loop
  while(!toProcessMain.empty()){
    connectedVertex = geneNetwork->getVertexForValue(
                                                toProcessMain.front());
    toProcessMain.pop();

    if(NULL == connectedVertex) continue;

    toReturn.push(connectedVertex->value.nameIndex);
    markConnectedVertexesTriple(connectedVertex, threeSigma, twoSigma,
                                            geneNetwork, toProcessMain);
    geneNetwork->removeVertex(connectedVertex);
  }

  return toReturn;
}


queue< queue<size_t> > tripleLink(graph<geneData, u8> *geneNetwork,
                                        const struct config &settings){
  queue< queue<size_t> > toReturn;

  removeWeakVerticies(geneNetwork, settings.threeSigmaAdj,
                                                  settings.twoSigmaAdj);

  while(geneNetwork->getNumEdges() > 0){
    toReturn.push(tripleLinkIteration(geneNetwork,
                        settings.threeSigmaAdj, settings.twoSigmaAdj));
    removeWeakVerticies(geneNetwork, settings.threeSigmaAdj,
                                                  settings.twoSigmaAdj);
  }

  return toReturn;
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
