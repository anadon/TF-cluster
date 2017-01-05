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


#include "auxillaryUtilities.hpp"
#include "correlation-matrix.hpp"
#include "diagnostics.hpp"
#include "edge.t.hpp"
#include "vertex.t.hpp"
#include "graph.t.hpp"
#include "tripleLink.hpp"
#include "upper-diagonal-square-matrix.t.hpp"

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
      args->tflist = arg;
      break;
    case 'e':
      args->exprData = arg;
      break;
    case 'k':
      test = atoi(arg);
      if(test < 1 || test > 255){
        cerr << "keep value out of bounds [1, 255]." << endl;
        exit(EINVAL);
      }
      args->keepTopN = (u8) test;
      break;
    case '1':
      args->threeSigma = atof(arg);
      if(0 >= args->threeSigma)
        exit(EINVAL);
      break;
    case '2':
      args->twoSigma = atof(arg);
      if(0 >= args->twoSigma)
        exit(EINVAL);
      break;
    case '3':
      args->oneSigma = atof(arg);
      if(0 >= args->oneSigma)
        exit(EINVAL);
      break;
    case 'c':
      args->corrMethod = arg;
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
  queue< queue<size_t> > result;
  CMF protoGraph;
  UpperDiagonalSquareMatrix<u8> *sccm;

  //parse input
  settings = config{0, 0, 0, 0.0, 0.0, 0.0, 0, 0, 0, 100};
  argp_parse(&interpreter, argc, argv, 0, 0, &settings);


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
  
  sccm = constructCoincidenceMatrix(protoGraph, settings);
  corrData = constructGraph(sccm, protoGraph, settings);
  delete sccm;

  result = tripleLink(corrData, settings);

  delete corrData;

  printClusters(result, protoGraph.TFLabels);

  return 0;
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
