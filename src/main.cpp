/*Copyright 2016-2017 Josh Marshall************************************/

/***********************************************************************
    This file is part of TF-Cluster.

    TF-Cluster is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TF-Cluster.  If not, see <http://www.gnu.org/licenses/>.
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
#include <graph.tpp>
#include <stdio.h>
#include <vector>


#include <auxillaryUtilities.hpp>
#include <tf-cluster-files.hpp>
#include <diagnostics.hpp>
#include <tripleLink.hpp>
#include <upper-diagonal-square-matrix.tpp>

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::cerr;
using std::endl;

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
      args->tflistFile = arg;
      break;
    case 'e':
      args->exprDataFile = arg;
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


config getDefaultSettings(){
  return config{0, "csv", 0, "tf-cluster-legacy", 0, 0.0, 0.0, 0.0, 0, 0, 0, 100};
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
  graph<geneData, u8> *corrData;
  struct config settings;
  queue< queue<size_t> > result;
  double **protoGraph;
  UpperDiagonalSquareMatrix<u8> *sccm;
  unordered_map<string, char> TFCheck;

  //parse input
  settings = getDefaultSettings();
  argp_parse(&interpreter, argc, argv, 0, 0, &settings);


  TFCheck = readInTFList(settings.tflistFile);
  ED fileData = loadExpressionData(settings.exprDataFile, TFCheck);
  protoGraph = generateMatrixFromFile(fileData, TFCheck, "spearman");



  if(NULL == protoGraph){
    cerr << "There was a fatal error in generating the correlation "
            "matrix" << endl;
  }

  if(settings.keepTopN >= fileData.numGenes){
    cerr << "Too few genes to perform an analysis." << endl;
    return 0;
  }


  sccm = constructCoincidenceMatrix((const double**)protoGraph, fileData, settings);

  //printCoincidenceMatrix(sccm, settings.keepTopN, protoGraph.TFLabels);

  corrData = constructGraph(sccm, fileData, settings);
  delete sccm;

  result = tripleLink(corrData, settings);

  delete corrData;

  printClusters(result, fileData.genes);

  return 0;
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
