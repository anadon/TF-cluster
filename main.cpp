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
  struct config settings;
  int error;
  vector< graph<geneData, double>* > result;
  struct loadFromFileReturn fileData;

  cout << "Loading configuration...";
  settings = loadConfig();
  cout << "Loaded" << endl;
  cout << "verifying input...";
  error = verifyInput(argc, argv, settings.geneListFile, settings.expressionFile);
  if(error){
    cout << "Invalid!" << endl;
    return error;
  }
  cout << "valid" << endl;

  cout << "Making primary graph structure...";
  corrData = new graph<geneData, double>();
  cout << "complete" << endl;
  
  cout << "Loading files and calculating matrices...";
  fileData = loadFromFile(settings.geneListFile, settings.expressionFile);
  cout << "calculated" << endl;
  
  cout << "Converting matrix entries to sigma values...";
  convertCoeffToSigmaValue(fileData);
  cout << "complete" << endl;
  
  cout << "Adding strongest edges to primary graph structure...";
  addTopEdges(corrData, fileData);
  cout << "complete" << endl;
  
  cout << "Done!" << endl;
  delete corrData;
  return 0;
  
  //printEdges(corrData);

  pruneGraph(corrData, settings.topPick);
  
  //convert the raw coefficients to their equivelant sigma values
  

  
  result = tripleLink(corrData, settings.tripleLink1, 
                            settings.tripleLink2, settings.tripleLink3);
  
  printClusters(result);
  
  delete corrData;
  for(size_t i = 0; i < result.size(); i++){
    delete result[i];
  }

  return 0;
}