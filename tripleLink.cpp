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


//static bool killSignal = false;


void markConnectedVertex(const double &edgeWeight, 
            vertex<geneData, double> *toMark, const double &tripleLink1, 
                  const double &tripleLink2, const double &tripleLink3){
  if(edgeWeight >= tripleLink1){
    if(!toMark->value.threeSigmaLink){
      toMark->value.threeSigmaLink = true;
    }else if(!toMark->value.twoSigmaLink){
      toMark->value.twoSigmaLink = true;
    }else{
      toMark->value.oneSigmaLink = true;
    }
  }else if(edgeWeight >= tripleLink2){
    if(!toMark->value.twoSigmaLink){
      toMark->value.twoSigmaLink = true;
    }else{
      toMark->value.oneSigmaLink = true;
    }
  }else if(edgeWeight >= tripleLink3){
    toMark->value.oneSigmaLink = true;
  }
}


graph<geneData, double>* tripleLinkIteration(graph<geneData, double> *geneNetwork,
                    const double tripleLink1, const double tripleLink2, 
                const double tripleLink3, const size_t targetEdgeIndex){
  graph<geneData, double> *toReturn;
  toReturn = new graph<geneData, double>();
  edge<geneData, double> *initialEdge;
  vertex<geneData, double> *firstVertex, *secondVertex;
  
  queue< vertex<geneData, double>* > toProcess;
  
  initialEdge = &geneNetwork->edgeArray[targetEdgeIndex];
  firstVertex = initialEdge->left;
  secondVertex = initialEdge->right;
  
  toReturn->addVertex(*firstVertex);
  toReturn->addVertex(*secondVertex);
  
  for(size_t i = 0; i < firstVertex->numEdges; i++){
    vertex<geneData, double> *connectedVertex;
      
    connectedVertex = firstVertex->edges[i]->other(firstVertex);
    
    markConnectedVertex(firstVertex->edges[i]->weight, connectedVertex,
                                tripleLink1, tripleLink2, tripleLink3);
  }
  geneNetwork->removeVertex(firstVertex);
  
  for(size_t i = 0; i < secondVertex->numEdges; i++){
    vertex<geneData, double> *connectedVertex;
      
    connectedVertex = secondVertex->edges[i]->other(secondVertex);
    
    markConnectedVertex(secondVertex->edges[i]->weight, connectedVertex,
                                tripleLink1, tripleLink2, tripleLink3);
    
    if(connectedVertex->value.threeSigmaLink && 
      connectedVertex->value.twoSigmaLink){
      toProcess.push(connectedVertex);
    }
  }
  geneNetwork->removeVertex(secondVertex);
  
  while(!toProcess.empty()){
    vertex<geneData, double> *connectedVertex;
    
    connectedVertex = toProcess.front();
    toProcess.pop();
    toReturn->addVertex(*connectedVertex);
    
    for(size_t i = 0; i < connectedVertex->numEdges; i++){
      vertex<geneData, double> *otherVertex;
      
      otherVertex = connectedVertex->edges[i]->other(connectedVertex);
      
      if(otherVertex->value.threeSigmaLink && 
         otherVertex->value.twoSigmaLink &&
         otherVertex->value.oneSigmaLink){
        continue;
      }
      
      markConnectedVertex(connectedVertex->edges[i]->weight, otherVertex,
                                tripleLink1, tripleLink2, tripleLink3);
    
      if(otherVertex->value.threeSigmaLink && 
         otherVertex->value.twoSigmaLink &&
         otherVertex->value.oneSigmaLink){
        toProcess.push(connectedVertex);
      }
    }
    
    geneNetwork->removeVertex(connectedVertex);
  }
  
  
  return toReturn;
}


vector< graph<geneData, double>* > tripleLink(graph<geneData, double> *geneNetwork, 
    const double tripleLink1, const double tripleLink2, const double tripleLink3){
  vector< graph<geneData, double>* > toReturn;
  vector< vertex<geneData, double>* > removePass;
  
  
  //because removal of edges will cause suffling anyways, it is cheaper
  //to perform a linear search each time for the maximum edge.  The
  //maximum edge must also have a sigma value >= 3.
  
  
  //remove edges that won't be used for triple link --> with sigma 
  //values < 1
  removeLowEdges(geneNetwork, 1);
  
  //remove verticies which lack the connection strength to ever be 
  //included in the triple link algorithm
  removeWeakVerticies(geneNetwork);
  
  while(geneNetwork->numVertexes > 0 && geneNetwork->numEdges > 0){
    bool foundStrongEdge = false;
    size_t targetEdgeIndex;
    double maxFoundValue = 3;
    
    for(size_t i = 0; i < geneNetwork->numEdges; i++){
      if(geneNetwork->edgeArray[i].weight >= maxFoundValue){
        maxFoundValue = geneNetwork->edgeArray[i].weight;
        targetEdgeIndex = i;
        foundStrongEdge = true;
      }
    }
    if(!foundStrongEdge) break;
    
    toReturn.push_back(tripleLinkIteration(geneNetwork, targetEdgeIndex, 
                                tripleLink1, tripleLink2, tripleLink3));
    
    //remove edges that won't be used for triple link --> with sigma 
    //values < 1
    removeLowEdges(geneNetwork, 1);
  
    //remove verticies which lack the connection strength to ever be 
    //included in the triple link algorithm
    removeWeakVerticies(geneNetwork);
    
  }
  
  return toReturn;
}

#endif