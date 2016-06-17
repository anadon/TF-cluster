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
  graph<geneData, f64> *corrData;
  struct config settings;
  int error;
  vector< graph<geneData, f64>* > result;
  struct upperDiagonalMatrix protoGraph;

  cout << "Loading configuration...";
  settings = loadConfig();
  cout << "Loaded" << endl;
  cout << "verifying input...";
  error = verifyInput(argc, argv, settings.geneListFile);
  if(error){
    cout << "invalid!" << endl;
    return error;
  }
  cout << "valid" << endl;

  cout << "Making primary graph structure..."; fflush(stdout);
  corrData = new graph<geneData, f64>();
  cout << "complete!" << endl;
  
  cout << "Loading correlation matrix..."; fflush(stdout);
  protoGraph = loadMatrix(settings.tripleLink3);
  cout << "done!" << endl;
  
  cout << "Sorting and selecting edges..."; fflush(stdout);
  protoGraph = sortWeights(protoGraph, settings.keepTopN);
  cout << "finished!" << endl;

  
  //result = tripleLink(corrData, settings.tripleLink1, 
  //                          settings.tripleLink2, settings.tripleLink3);
  
  //printClusters(result);
  
  delete corrData;
  for(size_t i = 0; i < result.size(); i++){
    delete result[i];
  }

  return 0;
}