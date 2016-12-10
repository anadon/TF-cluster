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
 * verify parsed input as complete and valid
 **********************************************************************/
int verifyInput(const struct config settings);


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


//TODO: add doc
size_t XYToW(csize_t &x, csize_t &y, size_t n);


//TODO: add doc

pair<u8, size_t>* countingSortHighToLow(pair<u8, size_t> *toSort, 
                                                            csize_t n);

////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

int verifyInput(const struct config settings){
  if(NULL == settings.tflist){
    cerr << "file path to transcription factor names not given" << endl;
    return EINVAL;
  }
  if(NULL == settings.exprData){
    cerr << "file path to expression data not given" << endl;
    return EINVAL;
  }
  if(NULL == settings.corrMethod){
    cerr << "correlation calculation method not given" << endl;
    return EINVAL;
  }
  if(0.0  == settings.threeSigma){
    cerr << "triple link 1 not set" << endl;
    return EINVAL;
  }
  if(0.0  == settings.twoSigma){
    cerr << "triple link 2 not set" << endl;
    return EINVAL;
  }
  if(0.0  == settings.oneSigma){
    cerr << "triple link 3 not set" << endl;
    return EINVAL;
  }
  if(settings.threeSigma <= settings.twoSigma){
    cerr << "triple link 1 is less than or equal to triple link 2" << endl;
    return EINVAL;
  }
  if(settings.twoSigma <= settings.oneSigma){
    cerr << "triple link 2 is less than or equal triple link 1" << endl;
    return EINVAL;
  }
  if(0    == settings.keepTopN){
    cerr << "number of top links to keep is not set" << endl;
    return EINVAL;
  }
  
  return 0;
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


/***********************************************************************
 * Returns for each TF an array of pairs containing an index to the gene
 * name and its correlation coefficient sorted from highest correlation
 * to lowest correlation.
 * ********************************************************************/
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


/***********************************************************************
 * Returns for each TF an array of pairs containing an index to the gene
 * name and its correlation coefficient sorted from highest correlation
 * to lowest correlation.
 * ********************************************************************/
void *sortCoindicenceMatrixHelper(void *arg){
  struct sortCoindicenceMatrixHelperStruct *args =
                        (struct sortCoindicenceMatrixHelperStruct*) arg;
  csize_t numerator = args->numerator;
  csize_t denominator = args->denominator;
  cu8 *coincidenceMatrix = args->coindicenceMatrix;
  csize_t n = args->n;
  pair<u8, size_t> **sortedCoincidenceMatrix;
  cu8 keepTopN = args->keepTopN;
  void *tmpPtr;
  pair<u8, size_t> *sortColumn;
  
  sortedCoincidenceMatrix = args->sortedCoincidenceMatrix;
  
  
  tmpPtr = malloc(sizeof(**sortedCoincidenceMatrix) * (n-1));
  sortColumn = (pair<u8, size_t>*) tmpPtr;


  for(size_t itr = (numerator * n) / denominator;
                      itr < ((numerator + 1) * n) / denominator; itr++){
    
    
      for(size_t j = 0; j < itr; j++){
        pair<u8, size_t> toAdd;
        toAdd = pair<u8, size_t>(coincidenceMatrix[XYToW(j, itr, n)], j);
        sortColumn[j] = toAdd;
      }
      for(size_t j = itr+1; j < n; j++){
        pair<u8, size_t> toAdd;
        toAdd = pair<u8, size_t>(coincidenceMatrix[XYToW(itr, j, n)], j);
        sortColumn[j-1] = toAdd;
      }
    
    pair<u8, size_t> *sorted;
    sorted = countingSortHighToLow(sortColumn, n-1);
    tmpPtr = realloc(sorted, sizeof(*sorted) * keepTopN);
    sorted = (pair<u8, size_t>*) tmpPtr;
    sortedCoincidenceMatrix[itr] = sorted;
  }
  
  free(sortColumn);

  return NULL;
}


/* turn X, Y coordinates into a 1-D index value for a upper-diagonal
 * matrix */
size_t XYToW(csize_t &x, csize_t &y, size_t n){
  n--;
  return (n*(n-1)/2) - (n-y)*((n-y)-1)/2 + x - y - 1;
  //or? (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
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


graph<geneData, u8>* constructGraph(const CMF &protoGraph,
                                                struct config &settings){
  graph<geneData, u8>* tr;
  void *tmpPtr;
  pair<f64, size_t> **intermediateGraph;
  pthread_t *workers;
  struct constructGraphHelperStruct *instructions;
  int *toIgnore;
  //f64 *coincidenceMatrix; //This can't represent the data in a close 
  //enough way to the originalities nuances.  Keep note here for legacy.
  u8 *coincidenceMatrix;
  pair<u8, size_t> **sortedCoincidenceMatrix;
  f64 sigma;
  size_t clen, sum;
  struct sortCoindicenceMatrixHelperStruct *sortInstructions;
  

  /*This is setting up for a general multithreading job dispatch*/
  csize_t n = protoGraph.numRows();
  cu8 actualNumEdges = (u8) settings.keepTopN < protoGraph.numCols() ?
                              settings.keepTopN : protoGraph.numCols();
  csize_t numCPUs = thread::hardware_concurrency() < n-1 ?
                    thread::hardware_concurrency() : n-1 ;
  //csize_t numCPUs = 1;

  csize_t UDMSize = (n * (n-1) / 2) -1;

  cerr << "allocating preliminary memory" << endl;

  tmpPtr = malloc(sizeof(*intermediateGraph) * n);
  intermediateGraph = (pair<f64, size_t>**) tmpPtr;
  tmpPtr = malloc(sizeof(*workers) * numCPUs);
  workers = (pthread_t*) tmpPtr;
  tmpPtr = malloc(sizeof(*instructions) * numCPUs);
  instructions = (struct constructGraphHelperStruct*) tmpPtr;

  /*Setting up directions for each thread*/
  for(size_t i = 0; i < numCPUs; i++){
    instructions[i] = {
      i,
      numCPUs,
      protoGraph.numRows(),
      protoGraph.numCols(),
      (cf64**) protoGraph.fullMatrix,
      intermediateGraph,
      actualNumEdges
    };
  }

  cerr << "sorting edges" << endl;
  /*See constructGraphHelper() for the results*/
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
  for(size_t i = 0; i < protoGraph.numRows(); i++)
    free(protoGraph.fullMatrix[i]);
  free(protoGraph.fullMatrix);

  cerr << "constructing coincidence matrix" << endl;

  /*The coincidence matrix logic is best detailed in the paper*/
  //TODO: detail coincidence matrix here

  vector<bool> checks[protoGraph.numRows()];

  for(size_t i = 0; i < protoGraph.numRows(); i++)
    checks[i] = vector<bool>(protoGraph.numCols(), false);

  for(size_t i = 0; i < protoGraph.numRows(); i++)
    for(size_t j = 0; j < actualNumEdges; j++)
      checks[i][intermediateGraph[i][j].second] = true;

  tmpPtr = calloc(sizeof(*coincidenceMatrix), UDMSize);
  coincidenceMatrix = (u8*) tmpPtr;

  for(size_t i = 0; i < protoGraph.numRows(); i++)
    for(size_t j = i+1; j < protoGraph.numRows(); j++)
      for(size_t k = 0; k < actualNumEdges; k++)
        if(checks[i][intermediateGraph[j][k].second])
          coincidenceMatrix[XYToW(j, i, n)]++;
  
  //printCoincidenceMatrix(coincidenceMatrix, n, actualNumEdges, 
  //                                                protoGraph.TFLabels);

  for(size_t i = 0; i < protoGraph.numRows(); i++)
    checks[i].clear();
  
  //Now the coincidence matrix needs to be sorted again in order to only
  //add the top keepN entries into the graph for consideration.
  
  cerr << "Sorting coincidence matrix" << endl;
  
  tmpPtr = malloc(sizeof(*sortInstructions) * numCPUs);
  sortInstructions = (struct sortCoindicenceMatrixHelperStruct*) tmpPtr;
  tmpPtr = malloc(sizeof(*sortedCoincidenceMatrix) * n);
  sortedCoincidenceMatrix = (pair<u8, size_t>**) tmpPtr;
  
  for(size_t i = 0; i <numCPUs; i++){
    sortInstructions[i] = {
      i,
      numCPUs,
      coincidenceMatrix,
      n,
      sortedCoincidenceMatrix,
      actualNumEdges
    };
  }
  
  if(numCPUs > 1){
    tmpPtr = malloc(sizeof(*workers) * numCPUs);
    workers = (pthread_t*) tmpPtr;
  
    for(size_t i = 0; i < numCPUs; i++)
      pthread_create(&workers[i], NULL, sortCoindicenceMatrixHelper, 
                                                  &sortInstructions[i]);
    for(size_t i = 0; i < numCPUs; i++)
      pthread_join(workers[i], (void**) &toIgnore);
      
    free(workers);
  }else{
    sortCoindicenceMatrixHelper((void*) sortInstructions);
  }
  
  free(sortInstructions);
  
  cerr << "Calculating statistics" << endl;

  sigma = 0;
  sum = 0;
  clen = n * actualNumEdges * 2;
  
  
  for(size_t i = 0; i < n; i++)
    for(size_t j = 0; j < actualNumEdges-1; j++)
      sum += sortedCoincidenceMatrix[i][j].first;
  sum += (n * settings.keepTopN);
  
  cf64 avg = sum / ((f64) clen);
  
  for(size_t i = 0; i < n; i++){
    for(size_t j = 0; j < actualNumEdges-1; j++){
      cf64 tmp = sortedCoincidenceMatrix[i][j].first - avg;
      sigma += (tmp * tmp);
    }
  }
  sigma += (n * (settings.keepTopN * settings.keepTopN));
  sigma = sqrt(sigma / ((f64)clen -1));

  cerr << "avg:\t" << avg << endl;
  cerr << "clen:\t" << clen << endl;
  cerr << "std:\t" << sigma << endl;
  
  //now change the sigma levels accordingly so that the storing matrix
  //for data can stay as 1-byte storage, not 8-byte.
  
  settings.threeSigma = (settings.threeSigma * sigma) + avg;
  settings.threeSigmaAdj = (u8) settings.threeSigma;
  settings.twoSigma = (settings.twoSigma * sigma) + avg;
  settings.twoSigmaAdj = (u8) settings.twoSigma;
  settings.oneSigma = (settings.oneSigma * sigma) + avg;
  settings.oneSigmaAdj = (u8) settings.oneSigma;
  
  cerr << "Adjusted sigmas are " << (int) settings.threeSigmaAdj << ", "
       << (int) settings.twoSigmaAdj << ", " << (int) settings.oneSigmaAdj << endl;
  
  //*/

  //now prepare the graph for all the data it is about to recieve, else
  //after the fact memory allocations can take minutes.
  cerr << "constructing graph" << endl;
  tr = new graph<geneData, u8>();

  tr->hintNumVertexes(protoGraph.numRows());
  tr->hintNumEdges(protoGraph.numRows() * actualNumEdges);

  for(size_t i = 0; i < n; i++)
    tr->addVertex(geneData(i))->hintNumEdges(actualNumEdges);

  for(size_t i = 0; i < protoGraph.numRows(); i++){
    vertex<geneData, u8> *left = tr->getVertexForValue(geneData(i));
    for(size_t j = 0; j < actualNumEdges; j++){
      u8 weight = sortedCoincidenceMatrix[i][j].first;
      if(weight >= settings.oneSigmaAdj){
        vertex<geneData, u8> *right;
        right = tr->getVertexForValue(geneData(sortedCoincidenceMatrix[i][j].second));
        if(false == left->areConnected(right))
          tr->addEdge(left, right, weight);
      }
    }
  }


  
  //for(size_t i = 0; i < tr->getNumEdges(); i++)
  //  if(tr->getEdges()[i]->weight != 0.0)
  //    clen++;
  /*
  clen = tr->getNumEdges();
  if(protoGraph.numRows() * 
  
  for(size_t i = 0; i < tr->getNumEdges(); i++)
    sum += abs(tr->getEdges()[i]->weight);
  
  cf64 avg = sum / clen;

  
  for(size_t i = 0; i < tr->getNumEdges(); i++)
    tr->getEdges()[i]->weight -= avg;
  
  for(size_t i = 0; i < tr->getNumEdges(); i++){
    cf64 tmp = tr->getEdges()[i]->weight;
    sigma += (tmp * tmp);
  }

  sigma = sqrt(sigma / (clen -1));
  
  for(size_t i = 0; i < tr->getNumEdges(); i++)
    tr->getEdges()[i]->weight /= sigma;
    
  cerr << "clen: " << clen << endl;
  cerr << "avg: " << avg << endl;
  cerr << "std: " << sigma << endl;
  
  
  for(size_t i = 0; i < tr->getNumEdges(); i++){
    edge<geneData, f64> *test = tr->getEdges()[i];
    if(oneSigma > test->weight){
      tr->removeEdge(test);
      i = 0;
    }
  }

  for(size_t i = 0; i < protoGraph.numRows(); i++)
    tr->getVertexForValue(geneData(i))->shrinkToFit();
  //*/

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


pair<u8, size_t>* countingSortHighToLow(pair<u8, size_t> *toSort, 
                                                            csize_t n){
  size_t counts[256];
  pair<u8, size_t> *sortSpace;
  void *tmpPtr;
  
  memset(counts, 0, sizeof(counts));
  
  tmpPtr = malloc(sizeof(*sortSpace) * n);
  sortSpace = (pair<u8, size_t>*) tmpPtr;
  
  for(size_t i = 0; i < n; i++)
    counts[toSort[i].first]++;
  for(int i = 255-1; i >= 0; i--)
    counts[i] += counts[i+1];
  
  for(size_t i = n-1; i != ((size_t)0)-1; i--){
    counts[toSort[i].first]--;
    sortSpace[counts[toSort[i].first]] = toSort[i];
  }
  
  return sortSpace;
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
