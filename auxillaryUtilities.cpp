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

////////////////////////////////////////////////////////////////////////
//INCLUDES//////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstring>
#include <fstream>
#include <graph.tpp>
#include <iostream>
#include <queue>
#include <string>
#include <thread>
#include <utility>

#include "auxillaryUtilities.hpp"
#include "diagnostics.hpp"
#include "statistics.h"

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


struct multithreadLoad{
  size_t numerator;
  size_t denominator;
  void *specifics;
};


struct constructGraphHelperStruct{
  size_t numRows;
  size_t numCols;
  cf64 **fullMatrix;
  pair<f64, size_t> **intermediateGraph;
};


struct constructSCCMHelperStruct{
  u8 numEdges;
  pair<f64, size_t> **intermediateGraph;
  unordered_map<size_t, bool> *hashChecks;
  pthread_mutex_t *rowLocks;
  UpperDiagonalSquareMatrix<u8> *coincidenceMatrix;
};


/*******************************************************************//**
 *  Struct to ease parameter passing from addTopEdges() to
 * addTopEdgesHelper()
 **********************************************************************/
struct addTopEdgesHelperStruct{
  struct correlationMatrix *protoGraph;
  u8 keepTopN;
  size_t startIndex;
  size_t endIndex;
};


//TODO
struct sortCoindicenceMatrixHelperStruct{
  UpperDiagonalSquareMatrix<u8> *coindicenceMatrix;
  size_t n;
  pair<u8, size_t>** sortedCoincidenceMatrix;
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


/*******************************************************************//**
 * Sorts toSort based on the values in toSort[x].first() in linear time.
 **********************************************************************/
pair<u8, size_t>* countingSortHighToLow(pair<u8, size_t> *toSort, 
                                                            csize_t n);


/***********************************************************************
 * TODO
 * ********************************************************************/
void *constructPreSCCMHelper(void *arg);


/***********************************************************************
 * TODO
 * ********************************************************************/
void *constructSCCMHelper(void *arg);


/***********************************************************************
 * Returns for each TF an array of pairs containing an index to the gene
 * name and its correlation coefficient sorted from highest correlation
 * to lowest correlation.
 * ********************************************************************/
void *sortCoindicenceMatrixHelper(void *arg);


/***********************************************************************
 * 
 * ********************************************************************/
void autoThreadLauncher(void* (*func)(void*), void *sharedArgs);


////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


void autoThreadLauncher(void* (*func)(void*), void *sharedArgs){
  void *tmpPtr;
  
  csize_t numCPUs = thread::hardware_concurrency();
  //csize_t numCPUs = 1;
  
  
  struct multithreadLoad *instructions;
  tmpPtr = malloc(sizeof(*instructions) * numCPUs);
  instructions = (struct multithreadLoad*) tmpPtr;
  
  for(size_t i = 0; i < numCPUs; i++){
    instructions[i] = {i, numCPUs, sharedArgs};
  }
  
  if(numCPUs < 2){
    func((void*) instructions);
  }else{
    int *toIgnore;
    pthread_t *workers;
    
    tmpPtr = malloc(sizeof(*workers) * numCPUs);
    workers = (pthread_t*) tmpPtr;
  
    for(size_t i = 0; i < numCPUs; i++)
      pthread_create(&workers[i], NULL, func, (void*) &instructions[i]);
      
    for(size_t i = 0; i < numCPUs; i++)
      pthread_join(workers[i], (void**) &toIgnore);
    
    free(workers);
  }
  
  free(instructions);
}


void *constructPreSCCMHelper(void *arg){
  struct multithreadLoad *argPrime = (struct multithreadLoad*) arg;
  csize_t numerator = argPrime->numerator;
  csize_t denominator = argPrime->denominator;
  
  struct constructGraphHelperStruct *args =
              (struct constructGraphHelperStruct*) argPrime->specifics;
  csize_t numRows = args->numRows;
  csize_t numCols = args->numCols;
  cf64 **fullMatrix = args->fullMatrix;
  pair<f64, size_t> **intermediateGraph = args->intermediateGraph;


  for(size_t i = (numerator * numRows) / denominator;
                    i < ((numerator + 1) * numRows) / denominator; i++){

    for(size_t j = 0; j < numCols; j++)
      intermediateGraph[i][j] = pair<f64, size_t>(fullMatrix[i][j], j);

    sortDoubleSizeTPairHighToLow(intermediateGraph[i], numCols);
  }

  return NULL;
}


void *constructSCCMHelper(void *arg){
  struct multithreadLoad *argPrime = (struct multithreadLoad*) arg;
  csize_t numerator = argPrime->numerator;
  csize_t denominator = argPrime->denominator;
  
  
  struct constructSCCMHelperStruct *args =
                (struct constructSCCMHelperStruct*) argPrime->specifics;
  cu8 numEdges = args->numEdges;
  pair<f64, size_t> **intermediateGraph = args->intermediateGraph;
  unordered_map<size_t, bool> *hashChecks = args->hashChecks;
  pthread_mutex_t *rowLocks = args->rowLocks;
  UpperDiagonalSquareMatrix<u8> *coincidenceMatrix = 
                                                args->coincidenceMatrix;
  
  
  csize_t itrStart = (numerator * 
              (coincidenceMatrix->numberOfElements()-1) ) / denominator;
  csize_t itrEnd = ((numerator + 1) * 
              (coincidenceMatrix->numberOfElements()-1) ) / denominator;
  
  const pair<size_t, size_t> startCoord = 
                                    coincidenceMatrix->WToXY(itrStart);
  const pair<size_t, size_t> endCoord = 
                                      coincidenceMatrix->WToXY(itrEnd);
  csize_t xStart = startCoord.first;
  csize_t yStart = startCoord.second;
  csize_t xEnd = endCoord.first;
  csize_t yEnd = endCoord.second;
  csize_t n = coincidenceMatrix->getSideLength();

  size_t j = xStart;
  for(size_t i = yStart; i <= yEnd; i++){
    for(; (j < n && i < yEnd) || j < xEnd; j++){
      for(size_t k = 0; k < numEdges; k++){
        csize_t target = intermediateGraph[j][k].second;
        if(hashChecks[i].count(target)){
          pthread_mutex_lock(&rowLocks[i]);
          u8 *ptr = coincidenceMatrix->getReferenceForIndex(i, j);
          (*ptr)++;
          pthread_mutex_unlock(&rowLocks[i]);
        }
      }
    }
    j = i+2;
  }

  return NULL;
}


void *sortCoindicenceMatrixHelper(void *arg){  
  struct multithreadLoad *argPrime = (struct multithreadLoad*) arg;
  csize_t numerator = argPrime->numerator;
  csize_t denominator = argPrime->denominator;
  
  
  struct sortCoindicenceMatrixHelperStruct *args =
        (struct sortCoindicenceMatrixHelperStruct*) argPrime->specifics;
  
  
  UpperDiagonalSquareMatrix<u8> *coincidenceMatrix = 
                                                args->coindicenceMatrix;
  csize_t n = args->n;
  pair<u8, size_t> **sortedCoincidenceMatrix = 
                                          args->sortedCoincidenceMatrix;
  
  void *tmpPtr;
  pair<u8, size_t> *sortColumn;

  
  tmpPtr = malloc(sizeof(**sortedCoincidenceMatrix) * (n-1));
  sortColumn = (pair<u8, size_t>*) tmpPtr;


  for(size_t itr = (numerator * n) / denominator;
                      itr < ((numerator + 1) * n) / denominator; itr++){
    
    
    for(size_t j = 0; j < itr; j++){
      sortColumn[j] = pair<u8, size_t>(coincidenceMatrix->getValueAtIndex(j, itr), j);
    }
    for(size_t j = itr+1; j < n; j++){
      sortColumn[j-1] = pair<u8, size_t>(coincidenceMatrix->getValueAtIndex(itr, j), j);
    }
    
    pair<u8, size_t> *sorted;
    sorted = countingSortHighToLow(sortColumn, n-1);
    sortedCoincidenceMatrix[itr] = sorted;
  }
  
  free(sortColumn);

  return NULL;
}


UpperDiagonalSquareMatrix<u8>* constructCoincidenceMatrix(const CMF &protoGraph,
                                                struct config &settings){

  void *tmpPtr;
  pair<f64, size_t> **intermediateGraph;
  UpperDiagonalSquareMatrix<u8> *coincidenceMatrix;
  

  /*This is setting up for a general multithreading job dispatch*/
  csize_t n = protoGraph.numRows();
  
  cu8 actualNumEdges = settings.keepTopN;

  //Allocating preliminary memory
  tmpPtr = malloc(sizeof(*intermediateGraph) * n);
  intermediateGraph = (pair<f64, size_t>**) tmpPtr;
  for(size_t i = 0; i < n; i++){
    tmpPtr = malloc(sizeof(**intermediateGraph) * protoGraph.numCols());
    intermediateGraph[i] = (pair<f64, size_t>*) tmpPtr;
    memset(intermediateGraph[i], 0, sizeof(**intermediateGraph) * protoGraph.numCols());
  }

  struct constructGraphHelperStruct preSCCMInstr;
  preSCCMInstr = {
      protoGraph.numRows(), 
      protoGraph.numCols(),
      (cf64**) protoGraph.fullMatrix, 
      intermediateGraph
    };

  autoThreadLauncher(constructPreSCCMHelper, (void*) &preSCCMInstr);
  

  //Don't need the very large UDMatrix in protoGraph; free it.
  for(size_t i = 0; i < protoGraph.numRows(); i++)
    free(protoGraph.fullMatrix[i]);
  free(protoGraph.fullMatrix);

  for(size_t i = 0; i < n; i++){
    size_t allocSize = sizeof(**intermediateGraph) * actualNumEdges;
    tmpPtr = realloc(intermediateGraph[i], allocSize);
    intermediateGraph[i] = (pair<f64, size_t>*) tmpPtr;
  }

  //Allocating coincidence matrix
  /*The coincidence matrix logic is best detailed in the paper*/  
  
  unordered_map<size_t, bool> *hashChecks;
  hashChecks = new unordered_map<size_t, bool>[n];

  for(size_t i = 0; i < n; i++){
    hashChecks[i].max_load_factor(0.5);
    hashChecks[i].reserve(actualNumEdges);
    for(size_t j = 0; j < actualNumEdges; j++){
      size_t target = intermediateGraph[i][j].second;
      pair<size_t, bool> toInsert(target, true);
      hashChecks[i].insert(toInsert);
    }
  }

  coincidenceMatrix = new UpperDiagonalSquareMatrix<u8>(n);
  coincidenceMatrix->zeroData();


  //Constructing coincidence matrix
  
  pthread_mutex_t *rowLocks;
  tmpPtr = malloc(sizeof(*rowLocks) * n);
  rowLocks = (pthread_mutex_t*) tmpPtr;
  for(size_t i = 0; i < n; i++){
    pthread_mutex_init(&rowLocks[i], NULL);
  };
  
  struct constructSCCMHelperStruct SCCMInstr;
  SCCMInstr = {
    actualNumEdges,
    intermediateGraph,
    hashChecks,
    rowLocks,
    coincidenceMatrix
  };
  
  autoThreadLauncher(constructSCCMHelper, (void*) &SCCMInstr);
      
  delete[] hashChecks;
  for(size_t i = 0; i < n; i++){
    pthread_mutex_destroy(&rowLocks[i]);
  };
  free(rowLocks);
  
  for(size_t i = 0; i < n; i++)
    free(intermediateGraph[i]);
  free(intermediateGraph);
  
  return coincidenceMatrix;
}


graph<geneData, u8>* constructGraph(UpperDiagonalSquareMatrix<u8> *SCCM,
                        const CMF &protoGraph, struct config &settings){
  graph<geneData, u8>* tr;
  void *tmpPtr;
  pthread_t *workers;
  int *toIgnore;
  pair<u8, size_t> **sortedCoincidenceMatrix;
  f64 sigma;
  size_t clen, sum;
  struct sortCoindicenceMatrixHelperStruct sortInstructions;
  
  csize_t n = protoGraph.numRows();
  
  cu8 actualNumEdges = (u8) settings.keepTopN;
  csize_t numCPUs = thread::hardware_concurrency();
  
  
  //Now the coincidence matrix needs to be sorted again in order to only
  //add the top keepN entries into the graph for consideration.
  
  //Sorting coincidence matrix
  tmpPtr = malloc(sizeof(*sortedCoincidenceMatrix) * n);
  sortedCoincidenceMatrix = (pair<u8, size_t>**) tmpPtr;
  
  sortInstructions = {
      SCCM, 
      n, 
      sortedCoincidenceMatrix
    };
  
  autoThreadLauncher(sortCoindicenceMatrixHelper, 
                                            (void*) &sortInstructions);
 
  
  for(size_t i = 0; i < n; i++){
    tmpPtr = realloc(sortedCoincidenceMatrix[i], 
                    actualNumEdges * sizeof(**sortedCoincidenceMatrix));
    sortedCoincidenceMatrix[i] = (pair<u8, size_t>*) tmpPtr;
  }
  
  
  //Calculating statistics
  sigma = 0;
  sum = 0;
  clen = n * actualNumEdges;
  
  
  for(size_t i = 0; i < n; i++)
    for(size_t j = 0; j < actualNumEdges; j++)
      sum += sortedCoincidenceMatrix[i][j].first;
  
  cf64 avg = sum / ((f64) clen);
  
  for(size_t i = 0; i < n; i++){
    for(size_t j = 0; j < actualNumEdges; j++){
      cf64 tmp = sortedCoincidenceMatrix[i][j].first - avg;
      sigma += (tmp * tmp);
    }
  }
  sigma = sqrt(sigma / ((f64)clen -1));

  cerr << "avg:\t" << avg << endl;
  cerr << "clen:\t" << clen << endl;
  cerr << "std:\t" << sigma << endl;
  
  //now change the sigma levels accordingly so that the storing matrix
  //for data can stay as 1-byte storage, not 8-byte.
  
  settings.threeSigma = (settings.threeSigma * sigma) + avg;
  settings.threeSigmaAdj = (u8) ceil(settings.threeSigma);
  settings.twoSigma = (settings.twoSigma * sigma) + avg;
  settings.twoSigmaAdj = (u8) ceil(settings.twoSigma);
  settings.oneSigma = (settings.oneSigma * sigma) + avg;
  settings.oneSigmaAdj = (u8) ceil(settings.oneSigma);
  
  //cerr << "Adjusted sigmas are " << (int) settings.threeSigmaAdj << ", "
  //     << (int) settings.twoSigmaAdj << ", " 
  //     << (int) settings.oneSigmaAdj << endl;


  //now prepare the graph for all the data it is about to recieve, else
  //after the fact memory allocations can take minutes.
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

  for(size_t i = 0; i < n; i++){
    free(sortedCoincidenceMatrix[i]);
  }
  free(sortedCoincidenceMatrix);
  
  tr->shrinkToFit();

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

  delete[] indiciesOfInterest;
  delete[] newIndiciesOfInterest;
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
