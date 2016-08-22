/*******************************************************************//**
         FILE:  auxillaryUtilities.cpp

  DESCRIPTION:  Miscelaneous functions used in TF-cluser

         BUGS:  Correlation values are larger than perl version
        NOTES:  ---
       AUTHOR:  Josh Marshall <jrmarsha@mtu.edu>
      COMPANY:  Michigan technological University
      VERSION:  See git log
      CREATED:  See git log
     REVISION:  See git log
     LISCENSE:  GPLv3
***********************************************************************/

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <queue>
#include <string>
#include <thread>
#include <utility>


#include "auxillaryUtilities.hpp"
#include "edge.t.hpp"
#include "graph.t.hpp"
#include "statistics.h"
#include "vertex.t.hpp"

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::pair;
using std::queue;
using std::string;
using std::thread;
using std::vector;


////////////////////////////////////////////////////////////////////////
//PRIVATE STRUCTS
////////////////////////////////////////////////////////////////////////

struct quickMergeDoubleSizeTPair_recurse_struct{
  pair<f64, size_t> *toSort;
  pair<f64, size_t> *workSpace;
  size_t size;
};


struct constructGraphHelperStruct{
  size_t numerator;
  size_t denominator;
  size_t numRows;
  size_t numCols;
  cf64 **fullMatrix;
  pair<f64, size_t> **intermediateGraph;
  u8 maxNumEdges;
};

////////////////////////////////////////////////////////////////////////
//PRIVATE FUNCTION DECLARATIONS/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  A helper function to sortDoubleSizeTPairLowToHigh().
 **********************************************************************/
void sortDoubleSizeTPairLowToHighHelper(pair<f64, size_t> *toSort,
                csize_t leftIndex, csize_t rightIndex, csize_t endIndex,
                                          pair<f64, size_t> *sortSpace);


/*******************************************************************//**
 *  A helper function to sortDoubleSizeTPairHighToLow().
 **********************************************************************/
void sortDoubleSizeTPairHighToLowHelper(pair<f64, size_t> *toSort,
                csize_t leftIndex, csize_t rightIndex, csize_t endIndex,
                                          pair<f64, size_t> *sortSpace);


/*******************************************************************//**
 *  A helper function to constructGraph() for edge selection.
 **********************************************************************/
void *constructGraph(void *arg);

////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//TODO: this can be made more complete
int verifyInput(int argc, char **argv, struct config settings){

  if( argc != 2 ) return EINVAL;

  int errorCode = 0;

  return errorCode;
}


struct config loadConfig(const char *configFilePath){
  struct config settings;
  ifstream configFile;
  string line;
  bool errorInFile = false;

  settings.keepTopN = 100;
  settings.expressionFile = "";
  settings.geneList = "";//labeled as geneList, but really contains TFs 
  settings.kickSize = 0;
  settings.threeSigma = -1;
  settings.twoSigma = -1;
  settings.oneSigma = -1;

  configFile.open(configFilePath);

  if(!configFile.is_open()){
    cerr << "ERROR: cannot open required file \"SCCM_pipe.cfg\""
         << endl;
    exit(1);
  }

  while(getline(configFile, line)){
    if(0 == line.size()) continue;
    if('#' == line[0]) continue;
    if(string::npos != line.find("expression")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting expression to " << readValue << endl;
      settings.expressionFile = readValue;
      continue;
    }
    if(string::npos != line.find("geneList")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting geneList to " << readValue << endl;
      settings.geneList = readValue;
      continue;
    }
    if(string::npos != line.find("topPick")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting topPick to " << readValue << endl;
      int readInt = atoi(readValue.c_str());
      if(0 > readInt || 255 < readInt){
        cerr << "Keep Top N value out of bounds (0,255]" << endl;
        exit(1);
      }
      settings.keepTopN = (u8) readInt;
      continue;
    }
    if(string::npos != line.find("kickSize")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting kickSize to " << readValue << endl;
      settings.kickSize = atoi(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink1")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting threeSigma to " << readValue << endl;
      settings.threeSigma = atof(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink2")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#") + 1;
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting twoSigma to " << readValue << endl;
      settings.twoSigma = atof(readValue.c_str());
      continue;
    }
    if(string::npos != line.find("tripleLink3")){
      size_t startIndex = line.find("=") + 1;
      size_t endIndex = line.find("#");
      if(string::npos == endIndex) endIndex = line.size();
      string readValue = line.substr(startIndex, endIndex - startIndex);
      cerr << "setting oneSigma to " << readValue << endl;
      settings.oneSigma = atof(readValue.c_str());
      continue;
    }
    cerr << "Skipped line \"" << line << "\"" << endl;
  }
  configFile.close();

  if(settings.expressionFile != ""){
    cerr << "expression is in \"" << settings.expressionFile << "\""
         << endl;
  }else{
    cerr << "ERROR: expression file not set" << endl;
    errorInFile = true;
  }

  cerr << "topPick is " << (short) settings.keepTopN << endl;

  cerr << "kickSize is " << settings.kickSize << endl;

  if(settings.threeSigma < 0){
    cerr << "ERROR: threeSigma is less than 0" << endl;
    errorInFile = true;
  }else if(settings.threeSigma < settings.twoSigma){
    cerr << "ERROR: threeSigma is smaller than twoSigma" << endl;
    errorInFile = true;
  }else{
    cerr << "threeSigma is " << settings.threeSigma << endl;
  }
  if(settings.twoSigma < 0){
    cerr << "ERROR: twoSigma is less than 0" << endl;
    errorInFile = true;
  }else if(settings.twoSigma < settings.oneSigma){
    cerr << "ERROR: twoSigma is smaller than oneSigma" << endl;
    errorInFile = true;
  }else{
    cerr << "twoSigma is " << settings.twoSigma << endl;
  }
  if(settings.oneSigma < 0){
    cerr << "ERROR: oneSigma is smaller than 0" << endl;
    errorInFile = true;
  }else{
    cerr << "oneSigma is " << settings.oneSigma << endl;
  }

  if(errorInFile) exit(1);


  return settings;
}


void printClusters(queue< queue<size_t> > clusters,
                                            const vector<string> &TFs){
  for(size_t i = 0; !clusters.empty(); i++){
    cout << "cluster: " << (i+1) << endl;
    while(!clusters.front().empty()){
      cout << TFs[clusters.front().front()] << endl;
      clusters.front().pop();
    }
    clusters.pop();
  }
}


void *constructGraphHelper(void *arg){
  struct constructGraphHelperStruct *args = 
                              (struct constructGraphHelperStruct*) arg;
  csize_t numerator = args->numerator;
  csize_t denominator = args->denominator;
  csize_t numRows = args->numRows;
  csize_t numCols = args->numCols;
  cf64 **fullMatrix = args->fullMatrix;
  pair<f64, size_t> **intermediateGraph = args->intermediateGraph;
  cu8 maxNumEdges = args->maxNumEdges;
  
  void *tmpPtr;
  
  for(size_t i = (numerator * numRows) / denominator; 
                    i < ((numerator + 1) * numRows) / denominator; i++){
    pair<f64, size_t> *toSort;
    size_t shrinkSize;
    
    tmpPtr = malloc(sizeof(*toSort) * (numCols));
    toSort = (pair<f64, size_t>*) tmpPtr;

    for(size_t j = 0; j < numCols; j++)
      toSort[j] = pair<f64, size_t>(fullMatrix[i][j], j);


    if(numCols < maxNumEdges){
      shrinkSize = numCols;
    }else{
      shrinkSize = maxNumEdges;
      sortDoubleSizeTPairHighToLow(toSort, numCols);
      tmpPtr = realloc(toSort, sizeof(*toSort) * shrinkSize);
      toSort = (pair<f64, size_t>*) tmpPtr;
    }

    intermediateGraph[i] = toSort;
  }
  
  return NULL;
}


size_t XYToW(csize_t x, csize_t y, size_t n){
  n--;
  return (n*(n-1)/2) - (n-y)*((n-y)-1)/2 + x - y - 1;
}


void printCoincidenceMatrix(cf64 *matrix, csize_t width, cu8 maxMatch, 
                                              const vector<string> TFs){
  f64 **mtr;
  mtr = (f64**) malloc(sizeof(*mtr) * width);
  for(size_t i = 0; i < width; i++)
    mtr[i] = (f64*) malloc(sizeof(**mtr) * width);
  
  for(size_t i = 0; i < width; i++){\
    mtr[i][i] = maxMatch;
    for(size_t j = i+1; j < width; j++){
      mtr[i][j] = mtr[j][i] = matrix[XYToW(i, j, width)];
    }
  }
  
  for(size_t i = 0; i < width; i++){
    cout << TFs[i] << "\t";
    for(size_t j = 0; j < width; j++){
      cout << mtr[i][j] << "\t";
    }
    cout << endl;
    free(mtr[i]);
  }
  free(mtr);
  cout << endl << endl;
}


graph<geneData, f64>* constructGraph(const CMF &protoGraph, 
                                        cf64 oneSigma, cu8 maxNumEdges){
  graph<geneData, f64>* tr;
  void *tmpPtr;
  pair<f64, size_t> **intermediateGraph;
  pthread_t *workers;
  struct constructGraphHelperStruct *instructions;
  int *toIgnore;
  f64 *coincidenceMatrix;
  f64 sigma;
  
  csize_t n = protoGraph.numRows;
  csize_t actualNumEdges = maxNumEdges < protoGraph.numCols ? 
                                      maxNumEdges : protoGraph.numCols;
  csize_t numCPUs = thread::hardware_concurrency() < n-1 ?
                    thread::hardware_concurrency() : n-1 ;
  csize_t UDMSize = (n * (n-1) / 2) -1;

  cerr << "allocating preliminary memory" << endl;
  
  tmpPtr = malloc(sizeof(*intermediateGraph) * n);
  intermediateGraph = (pair<f64, size_t>**) tmpPtr;
  tmpPtr = malloc(sizeof(*workers) * numCPUs);
  workers = (pthread_t*) tmpPtr;
  tmpPtr = malloc(sizeof(*instructions) * numCPUs);
  instructions = (struct constructGraphHelperStruct*) tmpPtr;

  for(size_t i = 0; i < numCPUs; i++){
    instructions[i] = {
      i, 
      numCPUs, 
      protoGraph.numRows,
      protoGraph.numCols,
      (cf64**) protoGraph.fullMatrix,
      intermediateGraph,
      maxNumEdges
    };
  }

  cerr << "sorting edges" << endl;

  if(numCPUs > 1){
    for(size_t i = 0; i < numCPUs; i++)
      pthread_create(&workers[i], NULL, constructGraphHelper,
                                                      &instructions[i]);
    for(size_t i = 0; i < numCPUs; i++)
      pthread_join(workers[i], (void**) &toIgnore);
  }else{
    constructGraphHelper((void*) instructions);
  }

  free(workers);
  free(instructions);
  
  //Don't need the very large UDMatrix in protoGraph; free it.
  for(size_t i = 0; i < protoGraph.numRows; i++)
    free(protoGraph.fullMatrix[i]);
  free(protoGraph.fullMatrix);

  cerr << "constructing correlation matrix" << endl;

  vector<bool> checks[protoGraph.numRows];
  
  for(size_t i = 0; i < protoGraph.numRows; i++)
    checks[i] = vector<bool>(protoGraph.numCols, false);

  for(size_t i = 0; i < protoGraph.numRows; i++)
    for(size_t j = 0; j < actualNumEdges; j++)
      checks[i][intermediateGraph[i][j].second] = true;

  tmpPtr = calloc(sizeof(*coincidenceMatrix), UDMSize);
  coincidenceMatrix = (f64*) tmpPtr;

  for(size_t i = 0; i < protoGraph.numRows; i++)
    for(size_t j = i+1; j < protoGraph.numRows; j++)
      for(size_t k = 0; k < actualNumEdges; k++)
        if(checks[i][intermediateGraph[j][k].second])
          coincidenceMatrix[(n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1]++;

  for(size_t i = 0; i < protoGraph.numRows; i++)
    checks[i].clear();

  cerr << "calculating sigma" << endl;
  sigma = 0;

  inplaceCenterMean(coincidenceMatrix, UDMSize);

  for(size_t i = 0; i < UDMSize; i++)
    sigma += (coincidenceMatrix[i] * coincidenceMatrix[i]);
  
  sigma = sqrt(sigma / (UDMSize-1));
  for(size_t i = 0; i < UDMSize; i++)
    coincidenceMatrix[i] /= sigma;

  //now prepare the graph for all the data it is about to recieve, else
  //after the fact memory allocations can take minutes.
  cerr << "constructing graph" << endl;
  tr = new graph<geneData, f64>();

  tr->hintNumVertexes(protoGraph.numRows);
  tr->hintNumEdges(protoGraph.numRows * actualNumEdges);

  for(size_t i = 0; i < n; i++)
    tr->addVertex(geneData(i))->hintNumEdges(actualNumEdges);

  for(size_t i = 0; i < protoGraph.numRows; i++){
    vertex<geneData, double> *left = tr->getVertexForValue(geneData(i));
    for(size_t j = i+1; j < protoGraph.numRows; j++){
      f64 weight = coincidenceMatrix[XYToW(j, i, n)];
      if(weight >= oneSigma){
        vertex<geneData, double> *right = tr->getVertexForValue(geneData(j));
        tr->addEdge(left, right, weight);
      }
    }
  }

  for(size_t i = 0; i < protoGraph.numRows; i++)
    tr->getVertexForValue(geneData(i))->shrinkToFit();

  return tr;
}


void sortDoubleSizeTPairHighToLow(pair<f64, size_t> *toSort,
                                                          csize_t size){
  size_t numRising;
  size_t i;
  void *tmpPtr;
  size_t *indiciesOfInterest, *newIndiciesOfInterest;

  if(1 >= size) return;

  numRising = 0;

  for(i = 0; i < size-1; i++){
    if(toSort[i].first < toSort[i+1].first) numRising++;
  }

  if(numRising > (size >> 1)){
    //reverse so that more are in order
    pair<f64, size_t> tmp;
    for(i = 0; i < size/2; i++){
      tmp = toSort[i];
      toSort[i] = toSort[(size-1) - i];
      toSort[(size-1) - i] = tmp;
    }
  }

  tmpPtr = malloc(sizeof(*indiciesOfInterest) * size);
  indiciesOfInterest = (size_t*) tmpPtr;
  indiciesOfInterest[0] = 0;
  size_t IOISize = 1;

  for(i = 0; i < size-1; i++){
    if(toSort[i].first < toSort[i+1].first)
      indiciesOfInterest[IOISize++] = (i+1);
  }
  indiciesOfInterest[IOISize++] = size;

  pair<f64, size_t> *sortSpace;
  tmpPtr = malloc(sizeof(*sortSpace) * size);
  sortSpace = (pair<f64, size_t>*) tmpPtr;

  tmpPtr = malloc(sizeof(*newIndiciesOfInterest) * size);
  newIndiciesOfInterest = (size_t*) tmpPtr;

  while(IOISize > 2){
    size_t NIOISize = 0;
    for(i = 0; i < IOISize-2; i+=2){

      sortDoubleSizeTPairHighToLowHelper(toSort, indiciesOfInterest[i],
                      indiciesOfInterest[i+1], indiciesOfInterest[i+2],
                                                            sortSpace);

      newIndiciesOfInterest[NIOISize++] = indiciesOfInterest[i];
    }
    if(!(IOISize & 1)){
      newIndiciesOfInterest[NIOISize++] = indiciesOfInterest[IOISize-2];
    }
    newIndiciesOfInterest[NIOISize++] = size;
    memcpy(indiciesOfInterest, newIndiciesOfInterest,
                                NIOISize * sizeof(*indiciesOfInterest));
    IOISize = NIOISize;
  }

  free(indiciesOfInterest);
  free(newIndiciesOfInterest);
  free(sortSpace);
  
  /*
  for(size_t i = 0; i < size-1; i++)
    if(toSort[i] < toSort[i+1])
      raise(SIGABRT); //*/
}


void sortDoubleSizeTPairHighToLowHelper(pair<f64, size_t> *toSort,
                csize_t leftIndex, csize_t rightIndex, csize_t endIndex,
                                          pair<f64, size_t> *sortSpace){
  size_t leftParser, rightParser, mergedParser;

  leftParser = leftIndex;
  rightParser = rightIndex;
  mergedParser = 0;
  while(leftParser < rightIndex && rightParser < endIndex)
    sortSpace[mergedParser++] =
        toSort[leftParser].first > toSort[rightParser].first ?
        toSort[leftParser++] : toSort[rightParser++];

  while(leftParser < rightIndex)
    sortSpace[mergedParser++] = toSort[leftParser++];
  while(rightParser < endIndex)
    sortSpace[mergedParser++] = toSort[rightParser++];
  memcpy(&toSort[leftIndex], sortSpace,
                              sizeof(*toSort) * (endIndex - leftIndex));
  
  /*
  for(size_t i = leftIndex; i < endIndex-1; i++)
    if(toSort[i] < toSort[i+1])
      raise(SIGABRT);//*/
}


void sortDoubleSizeTPairLowToHigh(pair<f64, size_t> *toSort,
                                                          csize_t size){
  size_t numFalling;
  size_t i;
  void *tmpPtr;

  if(1 >= size) return;

  numFalling = 0;

  for(i = 0; i < size-1; i++){
    if(toSort[i].first > toSort[i+1].first) numFalling++;
  }

  if(numFalling > (size >> 1)){
    //reverse so that more are in order
    pair<f64, size_t> tmp;
    for(i = 0; i < size/2; i++){
      tmp = toSort[i];
      toSort[i] = toSort[(size-1) - i];
      toSort[(size-1) - i] = tmp;
    }
  }

  size_t* indiciesOfInterest = new size_t[size];
  indiciesOfInterest[0] = 0;
  size_t IOISize = 1;

  for(i = 0; i < size-1; i++){
    if(toSort[i].first > toSort[i+1].first)
      indiciesOfInterest[IOISize++] = (i+1);
  }
  indiciesOfInterest[IOISize++] = size;

  pair<f64, size_t> *sortSpace;
  tmpPtr = malloc(sizeof(*sortSpace) * size);
  sortSpace = (pair<f64, size_t>*) tmpPtr;

  size_t *newIndiciesOfInterest = new size_t[size];
  while(IOISize > 2){
    size_t NIOISize = 0;
    for(i = 0; i < IOISize-2; i+=2){

      sortDoubleSizeTPairLowToHighHelper(toSort, indiciesOfInterest[i],
                      indiciesOfInterest[i+1], indiciesOfInterest[i+2],
                                                            sortSpace);

      newIndiciesOfInterest[NIOISize++] = indiciesOfInterest[i];
    }
    if(!(IOISize & 1)){
      newIndiciesOfInterest[NIOISize++] = indiciesOfInterest[IOISize-2];
    }
    newIndiciesOfInterest[NIOISize++] = size;
    memcpy(indiciesOfInterest, newIndiciesOfInterest,
                                NIOISize * sizeof(*indiciesOfInterest));
    IOISize = NIOISize;
  }

  delete indiciesOfInterest;
  delete newIndiciesOfInterest;
  free(sortSpace);

}


void sortDoubleSizeTPairLowToHighHelper(pair<f64, size_t> *toSort,
                csize_t leftIndex, csize_t rightIndex, csize_t endIndex,
                                          pair<f64, size_t> *sortSpace){
  size_t leftParser, rightParser, mergedParser;

  leftParser = leftIndex;
  rightParser = rightIndex;
  mergedParser = 0;
  while(leftParser < rightIndex && rightParser < endIndex)
    sortSpace[mergedParser++] =
        toSort[leftParser].first < toSort[rightParser].first ?
        toSort[leftParser++] : toSort[rightParser++];

  while(leftParser < rightIndex)
    sortSpace[mergedParser++] = toSort[leftParser++];
  while(rightParser < endIndex)
    sortSpace[mergedParser++] = toSort[rightParser++];
  memcpy(&toSort[leftIndex], sortSpace,
                              sizeof(*toSort) * (endIndex - leftIndex));
}


void inPlaceAbsoluteValue(f64 *array, csize_t size){
  for(size_t i = 0; i < size; i++)
    if(0 > array[i])
      array[i] = (-1.0) * array[i];
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////