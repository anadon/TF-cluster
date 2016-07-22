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


using std::cerr;


void printCorrelationMatrix(const struct correlationMatrix &protoGraph){
  //vector< pair<size_t, double>* > matrix;
  //unordered_map<string, size_t> labelLookup;
  //vector<string> labels;
  //vector<size_t> colSize;
  printf("\t");
  for(size_t i = 0; i < protoGraph.matrix.size();i++)
    printf("%s\t", protoGraph.labels[i].c_str());
  
  
  for(size_t i = 0; i < protoGraph.matrix.size(); i++){
    printf("%s\t", protoGraph.labels[i].c_str());
    for(size_t j = 0; j < protoGraph.colSize[j]; j++){
      printf("%lf\t", protoGraph.matrix[i][j].second);
    }
    printf("\n\n");
  }
  
}


bool isProtoGraphValid(const struct correlationMatrix &protoGraph){
  for(size_t i = 0; i < protoGraph.matrix.size(); i++){
    for(size_t j = 0; j < protoGraph.colSize[i]; j++){
      if(protoGraph.matrix.size() <= protoGraph.matrix[i][j].first) return false;
    }
  }
  return true;
}


int main(int argc, char **argv){
  graph<geneData, f64> *corrData;
  struct config settings;
  int error;
  queue< queue<size_t> > result;
  struct correlationMatrix protoGraph;

  cerr << "Loading configuration...";
  settings = loadConfig();
  cerr << "Loaded" << endl;
  cerr << "verifying input...";
  error = verifyInput(argc, argv, settings.geneListFile);
  if(error){
    cerr << "invalid!" << endl;
    return error;
  }
  cerr << "valid" << endl;
  
  cerr << "Loading correlation matrix..."; fflush(stderr);
  protoGraph = generateFullMatrixFromFile(settings.expressionFile, settings.tripleLink3);
  cerr << "done" << endl;
  
  cerr << "Printing correlation matrix...";
  printCorrelationMatrix(protoGraph);
  cerr << "done" << endl;
  return 0;
  
  /*if(!isProtoGraphValid(protoGraph)){
    cout << "given matrix is invalid!" << endl;
    raise(SIGABRT);
  }*/
  
  cerr << "Sorting and selecting edges..."; fflush(stderr);
  protoGraph = sortWeights(protoGraph, settings.keepTopN);
  cerr << "done" << endl;
  
  if(!isProtoGraphValid(protoGraph)){
    cerr << "sorted matrix is invalid!" << endl;
    raise(SIGABRT);
  }

  //cout << "Making primary graph structure..."; fflush(stdout);
  //corrData = new graph<geneData, f64>();
  //cout << "done" << endl;
  
  cerr << "Making graph..."; fflush(stderr);
  corrData = constructGraph(protoGraph);
  cerr << "done" << endl;
  
  for(size_t i = 0; i < protoGraph.matrix.size(); i++)
    free(protoGraph.matrix[i]);
  protoGraph.matrix.clear();
  protoGraph.labelLookup.clear();
  protoGraph.colSize.clear();
  
  cerr << "Pruning graph..."; fflush(stderr);
  pruneGraph(corrData, settings.keepTopN);
  cerr << "done" << endl;
  
  cerr << "Performing triple link..."; fflush(stderr);
  result = tripleLink(corrData, settings.tripleLink1, 
                            settings.tripleLink2);
  cerr << "done!" << endl;
  
  printClusters(result, protoGraph.labels);
  
  protoGraph.labels.clear();
  
  delete corrData;
  //for(size_t i = 0; i < result.size(); i++)
  //  delete result[i];
  

  return 0;
}