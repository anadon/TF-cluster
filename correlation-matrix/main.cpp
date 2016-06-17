/***********************************************************************
This is a C++ reimplementation of TF-Cluster.  The goal of this project
to improve the time and space requirements over the perl implementation.

Input: 1 file name or input from stdin.
***********************************************************************/
#include <stdio.h>
#include <vector>
#include <iostream>

#include "auxillaryUtilities.hpp"
#include "statistics.hpp"


using std::cerr;
using std::endl;


int main(int argc, char **argv){
  struct config settings;
  int error;
  struct loadFromFileReturn fileData;
  struct correlationMatrix corrMatrix;

  cerr << "Loading configuration...";
  settings = loadConfig();
  
  cerr << "Loaded" << endl << "verifying input...";
  error = verifyInput(argc, argv, settings.geneListFile, settings.expressionFile);
  if(error){
    cerr << "Invalid!" << endl;
    return error;
  }
  cerr << "valid" << endl;
  
  cerr << "Loading files..."; fflush(stdout);
  fileData = loadFromFile(settings.geneListFile, settings.expressionFile);
  cerr << "Done!" << endl;
  
  cerr << "Calculating matrices..."; fflush(stderr);
  corrMatrix = calculateCorrelationMatrix(fileData);
  cerr << "Complete!" << endl;
  
  for(size_t i = 0; i < fileData.numRows; i++)
    free(fileData.matrix[i]);
  free(fileData.matrix);
  
  cerr << "Converting matrix entries to sigma values..."; fflush(stdout);
  convertCoeffToSigmaValue(corrMatrix.UDMatrix, corrMatrix.numElements);
  cerr << "complete" << endl;
  
  cerr << "Printing edges...";  fflush(stdout);
  printEdges(fileData, corrMatrix);
  cerr << "Complete" << endl;
  
  free(corrMatrix.UDMatrix);
  
  return 0;
}