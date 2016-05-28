/***********************************************************************
This is a C++ reimplementation of TF-Cluster.  The goal of this project
to improve the time and space requirements over the perl implementation.

Input: 1 file name or input from stdin.
***********************************************************************/
#include <stdio.h>
#include <vector>

#include "auxillaryUtilities.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"
#include "graph.t.hpp"
#include "tripleLink.hpp"


int main(int argc, char **argv){
  graph<geneData, double> *corrData;
  double sigma;
  struct config settings;
  int error;
  vector< graph<geneData, double>* > result;

  error = verifyInput(argc, argv);
  if(error) return error;

  settings = loadConfig();

  corrData = new graph<geneData, double>();
  loadFromFile(corrData, settings.geneListFile, settings.expressionFile);
  
  //printEdges(corrData);

  pruneGraph(corrData, settings.topPick);
  sigma = calculateSigma(corrData->getEdges(), corrData->getNumEdges());
  
  //convert the raw coefficients to their equivelant sigma values
  convertCoeffToSigmaValue(corrData, sigma);
  

  
  result = tripleLink(corrData, settings.tripleLink1, 
                            settings.tripleLink2, settings.tripleLink3);
  
  printClusters(result);
  
  delete corrData;
  for(size_t i = 0; i < result.size(); i++){
    delete result[i];
  }

  return 0;
}