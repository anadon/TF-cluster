/***********************************************************************
This is a C++ reimplementation of TF-Cluster.  The goal of this project 
to improve the time and space requirements over the perl implementation.

Input: 1 file name or input from stdin.
***********************************************************************/
#include <stdio.h>

#include "auxillaryUtilities.hpp"
#include "nodeNetwork.hpp"



//nodeNetwork tripleLinkNetwork(pairTable toParse){
  
//}

//vector<vector<string>> tripleLinkResult(nodeNetwork toParse){
  
//}


double calculateSigma(const gEdge ** links, const size_t numEdges){
  double toReturn = 0;
  
  
  
  return toReturn;
}


int main(int argc, char **argv){
  
  int error = verifyInput(argc, argv);
  if(error) return error;
  
  nodeNetwork corrData;
  corrData.loadFromFile(stdin);
  fclose(stdin);
  
  corrData.prune();
  double sigma = calculateSigma(corrData.getEdges(), corrData.getNumEdges());
  //corrData.prepareTripleLink();
  //corrData.TripleLink();
  
  
  return 0;
}