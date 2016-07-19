#ifndef TRIPLELINK_HPP
#define TRIPLELINK_HPP

#include <queue>
#include <vector>

#include "auxillaryUtilities.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"
#include "graph.t.hpp"

using std::queue;
using std::vector;


/*
Return true if edge can be safely removed
*/
bool markConnectedVertexPrimer(cf64 edgeWeight, vertex<geneData, double> *toMark, 
                      cf64 high, cf64 med, queue<geneData> &toProcess){
  geneData *value = &toMark->value;
  
  if(value->threeSigmaLink & value->twoSigmaLink)
    return true;
  
  if(!value->threeSigmaLink && edgeWeight >= high)    value->threeSigmaLink = true;
  else if(!value->twoSigmaLink && edgeWeight >= med)  value->twoSigmaLink   = true;
  
  if(value->threeSigmaLink & value->twoSigmaLink){
    toProcess.push(toMark->value.nameIndex);
    return true;
  }else{
    return false;
  }
}


void markConnectedVertexesPrimer(vertex<geneData, double> *markFrom, cf64 high, 
                        cf64 med, graph<geneData, double> *geneNetwork, 
                                            queue<geneData> &toProcess){
  for(size_t i = 0; i < markFrom->getNumEdges(); i++)
    if(markConnectedVertexPrimer(markFrom->getEdges()[i]->weight,
                    markFrom->getEdges()[i]->other(markFrom), high, med, 
                                                            toProcess))
      geneNetwork->removeEdge(markFrom->getEdges()[i]);
}


/*
Return true if edge can be safely removed
*/
bool markConnectedVertex(cf64 edgeWeight, vertex<geneData, double> *toMark, 
                      cf64 high, cf64 med, queue<geneData> &toProcess){
  geneData *value = &toMark->value;
  
  if(value->threeSigmaLink & value->twoSigmaLink & value->oneSigmaLink)
    return true;
  
  if(!value->threeSigmaLink && edgeWeight >= high){    value->threeSigmaLink = true; }
  else if(!value->twoSigmaLink && edgeWeight >= med){  value->twoSigmaLink   = true; }
  else if(!value->oneSigmaLink){                       value->oneSigmaLink   = true; }
  
  if(value->threeSigmaLink & value->twoSigmaLink & value->oneSigmaLink){
    toProcess.push(toMark->value.nameIndex);
    return true;
  }else{
    return false;
  }
}


void markConnectedVertexes(vertex<geneData, double> *markFrom, cf64 high, 
                        cf64 med, graph<geneData, double> *geneNetwork, 
                                            queue<geneData> &toProcess){
  for(size_t i = 0; i < markFrom->getNumEdges(); i++)
    if(markConnectedVertex(markFrom->getEdges()[i]->weight,
            markFrom->getEdges()[i]->other(markFrom), high, med, toProcess))
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
      const vertex<geneData, f64> *target = geneNetwork->getVertexes()[i];
      if(2 > target->getNumEdges()){
        geneNetwork->removeVertex(target);
        disconnectedVerticiesFound = true;
        break;
      }
      
      bool highFound, medFound;
      highFound = medFound = false;
      for(size_t j = 0; j < target->getNumEdges() && (!highFound || !medFound); j++){
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
                                                  cf64 high, cf64 med){
  queue<size_t> toReturn;
  edge<geneData, double> *initialEdge;
  vertex<geneData, double> *firstVertex, *secondVertex, *connectedVertex;
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
  markConnectedVertexesPrimer(firstVertex, high, med, geneNetwork, toProcess);
  geneNetwork->removeVertex(firstVertex);
  
  markConnectedVertexesPrimer(secondVertex, high, med, geneNetwork, toProcess);
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

    markConnectedVertexes(connectedVertex, high, med, geneNetwork, toProcess);
    geneNetwork->removeVertex(connectedVertex);
  }

  return toReturn;
}


queue< queue<size_t> > tripleLink(graph<geneData, double> *geneNetwork,
                                                  cf64 high, cf64 med){
  queue< queue<size_t> > toReturn;

  //remove verticies which lack connection strength for inclusion
  removeWeakVerticies(geneNetwork, high, med);

  while(geneNetwork->getNumVertexes() > 0 && geneNetwork->getNumEdges() > 0){
    toReturn.push(tripleLinkIteration(geneNetwork, high, med));
    removeWeakVerticies(geneNetwork, high, med);
  }

  return toReturn;
}

#endif