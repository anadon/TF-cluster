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


/*void printCorrelationMatrix(const struct correlationMatrix &protoGraph){
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
}*/


//print graph to make sense of it's contents.
void printGraph(graph<geneData, f64> *corrData, const vector<string> &labels){
  for(size_t i = 0; i < corrData->getNumVertexes(); i++){
    fprintf(stderr, "%s\n[", labels[corrData->getVertexes()[i]->value.nameIndex].c_str());
    for(size_t j = 0; j < corrData->getVertexes()[i]->getNumEdges(); j++){
      fprintf(stderr, "{%s , %lf}\t", 
labels[corrData->getVertexes()[i]->getEdges()[j]->other(corrData->getVertexes()[i])->value.nameIndex].c_str(), 
                        corrData->getVertexes()[i]->getEdges()[j]->weight); 
    }
    fprintf(stderr, "]\n\n"
"========================================================================"
"\n\n");
  }
}


int main(int argc, char **argv){
  graph<geneData, f64> *corrData;
  struct config settings;
  int error;
  queue< queue<size_t> > result;
  struct UDCorrelationMatrix protoGraph;

  //TODO: this can be safer
  cerr << "Loading configuration...";
  settings = loadConfig(argv[1]);
  cerr << "Loaded" << endl;
  cerr << "verifying input...";
  error = verifyInput(argc, argv, settings.geneListFile);
  if(error){
    cerr << "invalid!" << endl;
    return error;
  }
  cerr << "valid" << endl;
  
  cerr << "Loading correlation matrix..."; fflush(stderr);
  protoGraph = generateUDMatrixFromFile(settings.expressionFile.c_str());
  convertCoeffToSigmaValue(protoGraph);
  cerr << "done" << endl;
  
  if(settings.keepTopN >= protoGraph.labels.size()){
    cerr << "Too few genes to perform an analysis." << endl;
    return 0;
  }
  
  /*cerr << "Printing correlation matrix...";
  printCorrelationMatrix(protoGraph);
  cerr << "done" << endl;*/
  
  cerr << "Making graph..."; fflush(stderr);
  corrData = constructGraph(protoGraph, settings.threeSigma, settings.oneSigma, settings.keepTopN);
  cerr << "done" << endl;
  
  cerr << "Pruning graph...";
  cerr << "done" << endl;
  
  //printGraph(corrData, protoGraph.labels);
  
  cerr << "Performing triple link..."; fflush(stderr);
  result = tripleLink(corrData, settings.threeSigma, 
                            settings.twoSigma);
  cerr << "done!" << endl;
  
  printClusters(result, protoGraph.labels);
  
  protoGraph.labels.clear();
  
  delete corrData;
  //for(size_t i = 0; i < result.size(); i++)
    //result[i].clear();
    //delete result[i];
  //result.clear();

  return 0;
}