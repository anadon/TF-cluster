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

#include <argp.h>
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
//GLOBAL VARIABLE DEFINITIONS///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

cs8 *argp_program_version = "0.9.0 RC1";


cs8 *argp_program_bug_address = "<jrmarsha@mtu.edu>";


argp_parser_t parser;


cs8 *doc = "tf-cluster is a tool used to predict transcriptional "
"networks using expression data and a list of expected transcription "
"factors.  Please refer to the original paper for more details at this "
"time.";


static struct argp_option options[] = {
  {"tf-list", 't', "FILE", 0, "File descriptor for list of transcription factors", 0},
  {"expression-data", 'e', "FILE", 0, "File descriptor for experiemntal data of gene expression", 0},
  {"keep", 'k', "INT", 0, "Number of highest matching genes matches between 1 and 255 for inclusion into the coexpression matrix", 0},
  {"triple-link-1", '1', "FLOAT", 0, "Highest link strength required for triple link.  Must be a positive number of standard deviations.", 0},
  {"triple-link-2", '2', "FLOAT", 0, "Middle link strength required for triple link.  Must be a positive number of standard deviations.", 0},
  {"triple-link-3", '3', "FLOAT", 0, "Lowest link strength required for triple link.  Must be a positive number of standard deviations.", 0},
  {"correlation", 'c', "STRING", 0, "Name of the statistical correlation to use in generating the correlation matrix.  Currently supports Pearson's correlation (pearson) and Spearman Rank (spearman).", 0},
  { 0 , 0, 0, 0, 0, 0}
};


static error_t parse_opt(int key, char *arg, struct argp_state *state){
  config *args;
  long test;
  args = (config*) state->input;
  switch(key){
    case 't':
      cerr << "tf-list" << " set to " << arg << endl;
      args->tflist = arg;
      break;
    case 'e':
      args->exprData = arg;
      cerr << "expression-data" << " set to " << arg << endl;
      break;
    case 'k':
      test = atoi(arg);
      if(test < 1 || test > 255){
        cerr << "keep value out of bounds [1, 255]." << endl;
        exit(EINVAL);
      }
      args->keepTopN = (u8) test;
      cerr << "keep" << " set to " << arg << endl;
      break;
    case '1':
      args->threeSigma = atof(arg);
      cerr << "triple-link-1" << " set to " << arg << endl;
      if(0 >= args->threeSigma)
        exit(EINVAL);
      break;
    case '2':
      args->twoSigma = atof(arg);
      cerr << "triple-link-2" << " set to " << arg << endl;
      if(0 >= args->twoSigma)
        exit(EINVAL);
      break;
    case '3':
      args->oneSigma = atof(arg);
      cerr << "triple-link-3" << " set to " << arg << endl;
      if(0 >= args->oneSigma)
        exit(EINVAL);
      break;
    case 'c':
      args->corrMethod = arg;
      cerr << "corrMethod" << " set to " << arg << endl;
      if(strcmp("pearson", arg)) break;
      if(strcmp("spearman", arg)) break;
      else{
        cerr << "Correlation method \"" << arg << "\" is not supported"
             << endl;
        exit(EINVAL);
      }
    default:
      return ARGP_ERR_UNKNOWN;
  }
  
  return 0;
}


static struct argp interpreter = {options, parse_opt, 0, doc, 0, 0, 0};


////////////////////////////////////////////////////////////////////////
//PRIVATE FUNCTION DEFINITIONS//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//*
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


/*
void graphvizRepresentation(graph<geneData, f64> *corrData, vector<string> labels){
  cout << "digraph rep{" << endl;
  cout << endl;
  cout << "bgcolor=\"transparent\"" << endl;
  cout << "fontcolor = black" << endl;
  cout << endl;
  for(size_t i = 0; i < corrData->getNumVertexes(); i++){
    cout << labels[corrData->getVertexes()[i]->value.nameIndex] << endl;
  }
  
  for(size_t i = 0; i < corrData->getNumVertexes(); i++){
    vertex<geneData, f64> *targetV;
    targetV = corrData->getVertexes()[i];
    for(size_t j = 0; j < targetV->getNumEdges(); j++){
      edge<geneData, f64> *targetE;
      targetE = corrData->getVertexes()[i]->getEdges()[j];
      if(targetV == targetE->left){
        cout << labels[targetE->left->value.nameIndex] << " -> ";
        cout << labels[targetE->right->value.nameIndex] << " ;" << endl;
      }
    }
  }
  
  cout << "}";
}//*/


//*
void printEdgeWeights(graph<geneData, u8> *corrData){
  for(size_t i = 0; i < corrData->getNumEdges(); i++)
    cout << (int) corrData->getEdges()[i]->weight << endl;
}//*/


//*
void printProtoGraph(const CMF &toPrint){
  for(size_t i = 0; i < toPrint.numRows(); i++){
    for(size_t j = 0; j < toPrint.numCols(); j++){
      fprintf(stdout, "%s\t%s\t%lf\n", toPrint.TFLabels[i].c_str(), 
                      toPrint.GeneLabels[j].c_str(), toPrint.fullMatrix[i][j]);
    }
  }
}//*/
  


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
//*
void printGraph(graph<geneData, f64> *toPrint,
                                          const vector<string> &labels){
  for(size_t i = 0; i < toPrint->getNumVertexes(); i++){
    for(size_t j = 0; j < toPrint->getVertexes()[i]->getNumEdges();
                                                                  j++){
      fprintf(stderr, "%s\t%s\t%lf\n",
          labels[toPrint->getVertexes()[i]->value.nameIndex].c_str(),
labels[toPrint->getVertexes()[i]->getEdges()[j]->other(toPrint->getVertexes()[i])->value.nameIndex].c_str(),
                    toPrint->getVertexes()[i]->getEdges()[j]->weight);
    }
    fprintf(stderr, " ]\n\n"
"========================================================================"
"\n\n");
  }
}//*/

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
  graph<geneData, u8> *corrData;
  struct config settings;
  int error;
  queue< queue<size_t> > result;
  CMF protoGraph;

  //parse input
  settings = config{0, 0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 100};
  argp_parse(&interpreter, argc, argv, 0, 0, &settings);
  if(0 != (error = verifyInput(settings))){
    return error;
  }

  cerr << "Loading correlation matrix" << endl;
  protoGraph = generateMatrixFromFile(settings.exprData, 
                                          settings.tflist, "spearman");
  if(NULL == protoGraph.fullMatrix){
    cerr << "There was a fatal error in generating the correlation "
            "matrix" << endl;
  }

  if(settings.keepTopN >= protoGraph.GeneLabels.size()){
    cerr << "Too few genes to perform an analysis." << endl;
    return 0;
  }
  //printProtoGraph(protoGraph);
  //exit(0);
  
  
  corrData = constructGraph(protoGraph, settings);

  //printEdgeWeights(corrData);
  cerr << "Performing triple link" << endl;
  result = tripleLink(corrData, settings);

  printClusters(result, protoGraph.TFLabels);

  delete corrData;

  return 0;
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
