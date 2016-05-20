#include <unistd.h>
#include <stdlib.h>

#include "nodeNetwork.hpp"
#include "gEdge.hpp"
#include "gNode.hpp"


nodeNetwork::nodeNetwork(){
  numNodes = numEdges = 0;
  nodeArray = (gNode**) NULL;
}


nodeNetwork::~nodeNetwork(){
  for(size_t i = 0; i < numNodes; i++)
    delete nodeArray[i];
  free(edgeArray);
  free(nodeArray);
}


int nodeNetwork::loadFromFile(FILE *input){

  char *line = NULL;
  vector<string> tokens;

  getline(&line, NULL, stdin);
  tokens = tokenizeString(line);

  for(size_t i = 0; i < tokens.size(); i++) addNode(tokens[i]);

  while(feof(stdin)){

    getline(&line, NULL, stdin);
    tokens = tokenizeString(line);

    for(size_t i = 0; i < tokens.size(); i+=2)
      addEdge(nodeArray[i], nodeArray[geneNameToNodeID.find(tokens[i])->second], atof(tokens[i+1].c_str()));
  }
  
  return 0;
}


void nodeNetwork::prune(u8 keepTopN){
  for(size_t i = 0; i < numNodes; i++) nodeArray[i]->prune(keepTopN);
}


void nodeNetwork::removeEdge(size_t edgeIndex){
  delete edgeArray[edgeIndex];
  edgeArray[edgeIndex] = edgeArray[--numEdges];
  edgeArray[edgeIndex]->edgeID = edgeIndex;

  edgeArray = (gEdge**) realloc(edgeArray, sizeof(*edgeArray) * numEdges);
}


void nodeNetwork::removeNode(size_t nodeIndex){
  delete nodeArray[nodeIndex];
  nodeArray[nodeIndex] = nodeArray[--numNodes];
  nodeArray[nodeIndex]->nodeID = nodeIndex;

  nodeArray = (gNode**) realloc(nodeArray, sizeof(*nodeArray) * numNodes);
}


void nodeNetwork::addNode(string geneID){
  size_t newSize = numNodes + 1;
  nodeArray = (gNode**) realloc(nodeArray, sizeof(*nodeArray) * newSize);
  geneNameToNodeID.emplace(geneID, numEdges);
  nodeArray[numNodes] = new gNode(this, numNodes);

  numNodes = newSize;
}


void nodeNetwork::addEdge(gNode *left, gNode *right, double coeff){
  unordered_map<string, size_t>::const_iterator potentialFind;

  potentialFind = geneNameToNodeID.find(right->geneID);
  if(potentialFind == geneNameToNodeID.end()) return;
  //TODO if found, exit early

  size_t newSize = numEdges + 1;
  edgeArray = (gEdge**) realloc(nodeArray, sizeof(*nodeArray) * newSize);
  geneNameToNodeID.emplace(left->geneID, numEdges);
  edgeArray[numNodes] = new gEdge(left, right, coeff);
}


inline const gEdge** nodeNetwork::getEdges(){
  return edgeArray;
}

inline const size_t nodeNetwork::getNumEdges(){
  return numEdges;
}