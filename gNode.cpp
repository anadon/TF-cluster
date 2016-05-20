#include <string.h>

#include "gNode.hpp"
#include "gEdge.hpp"
#include "nodeNetwork.hpp"


gNode::gNode(nodeNetwork *belongTo, size_t ID){
  nodeID = ID;
  net = belongTo;
  numEdges = 0;
  clusterGroup = 0;
  edges = (gEdge**) NULL;
}


gNode::~gNode(){
  for(size_t i = 0; i < numEdges; i++)
    net->removeEdge(edges[i]->edgeID);
}


size_t gNode::addEdge(gEdge *toRegister){
  size_t toReturn, newCount;

  toReturn = numEdges;
  newCount = numEdges + 1;

  edges = (gEdge**) realloc(edges, newCount * sizeof(*edges));
  edges[numEdges] = toRegister;
  numEdges = newCount;

  return toReturn;
}


void gNode::removeEdge(gEdge *toRemove){
  gEdge *tmp;
  size_t targetEdgeIndex;
  numEdges--;
  
  //check if this is the right or left of the edge, error if edge
  //doesn't connect this node/vertex
  if(toRemove == edges[toRemove->leftEdgeIndex]){
    targetEdgeIndex = toRemove->leftEdgeIndex;
  }else if(toRemove == edges[toRemove->rightEdgeIndex]){
    targetEdgeIndex = toRemove->rightEdgeIndex;
  }else{
    fprintf(stderr, "gNode::removeEdge: passed edges to remove is not"
        " connected to this vertex\n");
    fflush(stderr);
    exit(1);
  }
  
  //swap edge this is removing to the end
  tmp = edges[numEdges];
  edges[numEdges] = edges[targetEdgeIndex];
  edges[targetEdgeIndex] = tmp;
  
  //update the swapped edge's location in this structure so it still
  //knows where it is in this vertex/node
  if(this == edges[targetEdgeIndex]->left)
    edges[targetEdgeIndex]->leftEdgeIndex = targetEdgeIndex;
  else
    edges[targetEdgeIndex]->rightEdgeIndex = targetEdgeIndex;

  //deallocate the last edge pointer (but don't actually delete -- 
  //that's the network's job).
  edges = (gEdge**) realloc(edges, numEdges * sizeof(*edges));
}


//TODO make multithreading safe
//Quick merge, keep upper values
//this is not multithreading safe, but can be made safe.
void gNode::prune(u8 keepTopN){
  quickMergeEdges(edges, numEdges);
  
  for(size_t i = numEdges-1; i >= keepTopN; i--)
    net->removeEdge(edges[i]->edgeID);
}


void gNode::quickMergeEdges(gEdge **toSort, const size_t size){
  size_t numRising;
  size_t i;
  
  numRising = 0;
  
  for(i = 0; i < size-1; i++)
    if(toSort[i]->coeff < toSort[i+1]->coeff) numRising++;
  
  if(numRising > (size >> 1)){
    //reverse so that more are in order
    gEdge *tmp;
    for(i = 0; i < size/2; i++){
      tmp = toSort[i];
      toSort[i] = toSort[(size-1) - i];
      toSort[(size-1) - i] = tmp;
    }
  }
  
  vector<size_t> indiciesOfInterest;
  indiciesOfInterest.push_back(0);
  
  i=0;
  while(i < size-1)
    if(toSort[i]->coeff < toSort[i+1]->coeff)
      indiciesOfInterest.push_back(++i);
  indiciesOfInterest.push_back(size);
  
  while(indiciesOfInterest.size() > 2){
    vector<size_t> newIndiciesOfInterest;
    for(i = 0; i < indiciesOfInterest.size()-2; i+=2){
      mergeHelper(toSort, indiciesOfInterest[i], 
                  indiciesOfInterest[i+1],
                  indiciesOfInterest[i+2]);
      newIndiciesOfInterest.push_back(indiciesOfInterest[i]);
    }
    indiciesOfInterest = newIndiciesOfInterest;
    indiciesOfInterest.push_back(size);
  }
  
}


void gNode::mergeHelper(gEdge **toSort, const size_t leftIndex, 
                        const size_t rightIndex, const size_t endIndex){
  size_t leftParser, rightParser, mergedParser;
  gEdge **sortSpace;
  
  sortSpace = (gEdge**) malloc(sizeof(*sortSpace) * (endIndex - leftIndex));
  
  leftParser = leftIndex;
  rightParser = rightIndex;
  mergedParser = 0;
  while(leftParser < rightIndex && rightParser < endIndex)
    sortSpace[mergedParser++] = 
        toSort[leftParser] > toSort[rightParser] ? 
        toSort[leftParser++] : toSort[rightParser++];
  
  memcpy(&toSort[leftIndex], sortSpace, sizeof(*toSort) * (endIndex - leftIndex));
  
}
