/***********************************************************************
This is a C++ reimplementation of TF-Cluster.  The goal of this project
to improve the time and space requirements over the perl implementation.

Input: 1 file name or input from stdin.
***********************************************************************/
#include <stdio.h>
#include <vector>

#include "correlation-matrix.hpp"//Make this more formally external

#include "auxillaryUtilities.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"
#include "graph.t.hpp"
#include "tripleLink.hpp"


void printCorrelationMatrix(const struct correlationMatrix &protoGraph){
  //vector< pair<size_t, double>* > matrix;
  //unordered_map<string, size_t> labelLookup;
  //vector<string> labels;
  //vector<size_t> colSize;
  
  for(size_t i = 0; i < protoGraph.matrix.size(); i++){
    printf("%s\t", protoGraph.labels[i].c_str());
    for(size_t j = 0; j < protoGraph.colSize[j]; j++){
      printf("%lf\t", protoGraph.matrix[i][j].second);
    }
    printf("\n");
  }
  
}


int main(int argc, char **argv){
  graph<geneData, f64> *corrData;
  struct config settings;
  int error;
  vector< graph<geneData, f64>* > result;
  struct correlationMatrix protoGraph;

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
  
  cout << "Loading correlation matrix..."; fflush(stdout);
  protoGraph = generateFullMatrixFromFile(settings.expressionFile, settings.tripleLink3);
  cout << "done!" << endl;

  cout << "Making primary graph structure..."; fflush(stdout);
  corrData = new graph<geneData, f64>();
  cout << "complete!" << endl;
  
  cout << "Sorting and selecting edges..."; fflush(stdout);
  protoGraph = sortWeights(protoGraph, settings.keepTopN);
  cout << "finished!" << endl;
  
  //FREE protoGraph
  
  //printCorrelationMatrix(protoGraph);
  
  cout << "Making graph..."; fflush(stdout);
  corrData = constructGraph(protoGraph);
  cout << "complete!" << endl;
  
  cout << "Performing triple link..."; fflush(stdout);
  result = tripleLink(corrData, settings.tripleLink1, 
                            settings.tripleLink2, settings.tripleLink3);
  cout << "done!" << endl;
  
  //printClusters(result);
  
  delete corrData;
  //for(size_t i = 0; i < result.size(); i++){
  //  delete result[i];
  //}

  return 0;
}