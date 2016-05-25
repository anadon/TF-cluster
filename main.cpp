/***********************************************************************
This is a C++ reimplementation of TF-Cluster.  The goal of this project
to improve the time and space requirements over the perl implementation.

Input: 1 file name or input from stdin.
***********************************************************************/
#include <stdio.h>

#include "auxillaryUtilities.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"
#include "graph.t.hpp"


int main(int argc, char **argv){

  int error = verifyInput(argc, argv);
  if(error) return error;

  graph<struct geneData, double> corrData;
  corrData = loadFromFile(stdin);
  fclose(stdin);

  pruneGraph(corrData, 50);
  double sigma = calculateSigma(corrData.getEdges(), corrData.getNumEdges());
  //corrData.prepareTripleLink();
  //corrData.TripleLink();


  return 0;
}