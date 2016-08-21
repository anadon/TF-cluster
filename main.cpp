/*******************************************************************//**
         FILE:  main.cpp

         BUGS:  Correlation values are larger than perl version
        NOTES:  ---
       AUTHOR:  Josh Marshall <jrmarsha@mtu.edu>
      COMPANY:  Michigan technological University
      VERSION:  See git log
      CREATED:  See git log
     REVISION:  See git log
     LISCENSE:  GPLv3
***********************************************************************/
/***********************************************************************
This is a C++ reimplementation of TF-Cluster.  The goal of this project
to improve the time and space requirements over the perl implementation.

Input: 1 file name or input from stdin.
***********************************************************************/

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <vector>

#include "correlation-matrix.hpp"//Make this more formally external

#include "auxillaryUtilities.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"
#include "graph.t.hpp"
#include "tripleLink.hpp"

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::cerr;

////////////////////////////////////////////////////////////////////////
//PRIVATE FUNCTION DEFINITIONS//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


/*
void printCorrelationMatrix(const CMF &protoGraph){
  //vector< pair<size_t, double>* > matrix;
  //unordered_map<string, size_t> labelLookup;
  //vector<string> labels;
  //vector<size_t> colSize;
  printf("\t");
  for(size_t i = 0; i < protoGraph.GeneLabels.size();i++)
    printf("%s\t", protoGraph.GeneLabels[i].c_str());


  for(size_t i = 0; i < protoGraph.TFLabels.size(); i++){
    printf("%s\t", protoGraph.TFLabels[i].c_str());
    for(size_t j = 0; j < protoGraph.GeneLabels.size(); j++){
      printf("%lf\t", protoGraph.fullMatrix[i][j]);
    }
    printf("\n");
  }

}//*/

/*
bool isProtoGraphValid(const struct correlationMatrix &protoGraph){
  for(size_t i = 0; i < protoGraph.matrix.size(); i++){
    for(size_t j = 0; j < protoGraph.colSize[i]; j++){
      if(protoGraph.matrix.size() <= protoGraph.matrix[i][j].first)
                                                           return false;
    }
  }
  return true;
}*/



/*******************************************************************//**
 * print graph to make sense of it's contents to stderr.  The graph is
 * printed with the vertex name on a line, followed by a pair of square
 * brackets ('[ ' ' ]')  which contain curly bracket pairs which hold
 * the name of a connected vertex followed by the edge weight.  Each of
 * these entries is followed by a blank line.
 *
 * @param[in] toPrint The graph structure to print
 * @param[in] labels
 **********************************************************************/
void printGraph(graph<geneData, f64> *toPrint,
                                          const vector<string> &labels){
  for(size_t i = 0; i < toPrint->getNumVertexes(); i++){
    fprintf(stderr, "%s\n[ ",
          labels[toPrint->getVertexes()[i]->value.nameIndex].c_str());
    for(size_t j = 0; j < toPrint->getVertexes()[i]->getNumEdges();
                                                                  j++){
      fprintf(stderr, "{%s , %lf}\t",
labels[toPrint->getVertexes()[i]->getEdges()[j]->other(toPrint->getVertexes()[i])->value.nameIndex].c_str(),
                    toPrint->getVertexes()[i]->getEdges()[j]->weight);
    }
    fprintf(stderr, " ]\n\n"
"========================================================================"
"\n\n");
  }
}

////////////////////////////////////////////////////////////////////////
//PUBLIC FUNCTION DEFINITIONS///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 * Program entry point.
 *
 * @param[in] argc number of c-strings in argv
 * @param[in] argv arguments passed to program
 **********************************************************************/
int main(int argc, char **argv){
  graph<geneData, f64> *corrData;
  struct config settings;
  int error;
  queue< queue<size_t> > result;
  CMF protoGraph;

  //TODO: this can be safer
  cerr << "Loading configuration...";
  settings = loadConfig(argv[1]);
  cerr << "Loaded" << endl;
  
  cerr << "verifying input" << endl;
  error = verifyInput(argc, argv, settings);
  if(error){
    cerr << "invalid!" << endl;
    return error;
  }

  cerr << "Loading correlation matrix" << endl;
  protoGraph = generateMatrixFromFile(settings.expressionFile.c_str(), 
                                            settings.geneList.c_str());

  if(settings.keepTopN >= protoGraph.GeneLabels.size()){
    cerr << "Too few genes to perform an analysis." << endl;
    return 0;
  }
  
  printCorrelationMatrix(protoGraph);

  corrData = constructGraph(protoGraph, settings.oneSigma, 
                                                    settings.keepTopN);

  cerr << "Performing triple link" << endl;
  result = tripleLink(corrData, settings.threeSigma, settings.twoSigma);

  printClusters(result, protoGraph.TFLabels);

  delete corrData;

  return 0;
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////