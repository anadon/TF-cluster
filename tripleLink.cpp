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


graph<struct geneData, double>* tripleLinkIteration(graph<struct geneData, double> *geneNetwork,
                                                size_t targetEdgeIndex){
  graph<struct geneData, double> *toReturn;
  toReturn = new graph<struct geneData, double>();
  edge<struct geneData, double> *initialEdge;
  vertex<struct geneData, double> *firstVertex, *secondVertex;
  
  queue< vertex<struct geneData, double>* > toProcess;
  
  initialEdge = geneNetwork->edgeArray[targetEdgeIndex];
  firstVertex = initialEdge->left;
  secondVertex = initialEdge->right;
  
  toReturn->addVertex(*firstVertex);
  toReturn->addVertex(*secondVertex);
  
  for(size_t i = 0; i < firstVertex->numEdges; i++){
    vertex<struct geneData, double> *connectedVertex;
      
    connectedVertex = firstVertex->edges[i]->other(firstVertex);
      
    if(firstVertex->edges[i]->weight >= 3){
      if(!connectedVertex->value.threeSigmaLink){
        connectedVertex->value.threeSigmaLink = true;
      }else if(!connectedVertex->value.twoSigmaLink){
        connectedVertex->value.twoSigmaLink = true;
      }else{
        connectedVertex->value.oneSigmaLink = true;
      }
    }else if(firstVertex->edges[i]->weight >= 2){
      if(!connectedVertex->value.twoSigmaLink){
        connectedVertex->value.twoSigmaLink = true;
      }else{
        connectedVertex->value.oneSigmaLink = true;
      }
    }else{
      connectedVertex->value.oneSigmaLink = true;
    }
  }
  geneNetwork->removeVertex(firstVertex);
  
  for(size_t i = 0; i < secondVertex->numEdges; i++){
    vertex<struct geneData, double> *connectedVertex;
      
    connectedVertex = secondVertex->edges[i]->other(secondVertex);
      
    if(secondVertex->edges[i]->weight >= 3){
      if(!connectedVertex->value.threeSigmaLink){
        connectedVertex->value.threeSigmaLink = true;
      }else if(!connectedVertex->value.twoSigmaLink){
        connectedVertex->value.twoSigmaLink = true;
      }else{
        connectedVertex->value.oneSigmaLink = true;
      }
    }else if(secondVertex->edges[i]->weight >= 2){
      if(!connectedVertex->value.twoSigmaLink){
        connectedVertex->value.twoSigmaLink = true;
      }else{
        connectedVertex->value.oneSigmaLink = true;
      }
    }else{
      connectedVertex->value.oneSigmaLink = true;
    }
    
    if(connectedVertex->value.threeSigmaLink && 
      connectedVertex->value.twoSigmaLink){
      toProcess.push(connectedVertex);
    }
  }
  geneNetwork->removeVertex(secondVertex);
  
  while(!toProcess.empty()){
    vertex<struct geneData, double> *connectedVertex;
    
    connectedVertex = toProcess.front();
    toProcess.pop();
    toReturn->addVertex(*connectedVertex);
    
    for(size_t i = 0; i < connectedVertex->numEdges; i++){
      vertex<struct geneData, double> *otherVertex;
      
      otherVertex = connectedVertex->edges[i]->other(connectedVertex);
      
      if(otherVertex->value.threeSigmaLink && 
         otherVertex->value.twoSigmaLink &&
         otherVertex->value.oneSigmaLink){
        continue;
      }
      
      if(connectedVertex->edges[i]->weight >= 3){
        if(!otherVertex->value.threeSigmaLink){
          otherVertex->value.threeSigmaLink = true;
        }else if(!otherVertex->value.twoSigmaLink){
          otherVertex->value.twoSigmaLink = true;
        }else{
          otherVertex->value.oneSigmaLink = true;
        }
      }else if(connectedVertex->edges[i]->weight >= 2){
        if(!otherVertex->value.twoSigmaLink){
          otherVertex->value.twoSigmaLink = true;
        }else{
          otherVertex->value.oneSigmaLink = true;
        }
      }else{
        otherVertex->value.oneSigmaLink = true;
      }
    
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


vector< graph<struct geneData, double>* > tripleLink(graph<struct geneData, double> *geneNetwork, 
    double tripleLink1, double tripleLink2, double tripleLink3){
  vector< graph<struct geneData, double>* > toReturn;
  vector< vertex<struct geneData, double>* > removePass;
  
  
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
      if(geneNetwork->edgeArray[i]->weight >= maxFoundValue){
        maxFoundValue = geneNetwork->edgeArray[i]->weight;
        targetEdgeIndex = i;
        foundStrongEdge = true;
      }
    }
    if(!foundStrongEdge) break;
    
    toReturn.push_back(tripleLinkIteration(geneNetwork, targetEdgeIndex));
    
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