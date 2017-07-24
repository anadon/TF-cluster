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

#ifndef QUICKMERGE_HPP
#define QUICKMERGE_HPP


#include <utility>
#include <cstdlib>

using std::pair;
using std::size_t;


////////////////////////////////////////////////////////////////////////
//FUNCTION DECLARATIONS/////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


void _sortDoubleSizeTPairHighToLow(pair<f64, size_t> *toSort,
                                                          csize_t size);


void _sortDoubleSizeTPairHighToLowHelper(pair<f64, size_t> *toSort,
                csize_t leftIndex, csize_t rightIndex, csize_t endIndex,
                                          pair<f64, size_t> *sortSpace);


void _sortDoubleSizeTPairLowToHigh(pair<f64, size_t> *toSort,
                                                          csize_t size);


void _sortDoubleSizeTPairLowToHighHelper(pair<f64, size_t> *toSort,
                csize_t leftIndex, csize_t rightIndex, csize_t endIndex,
                                          pair<f64, size_t> *sortSpace);



////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void _sortDoubleSizeTPairHighToLow(pair<f64, size_t> *toSort,
                                                          csize_t size){
  size_t numRising;
  size_t i;
  void *tmpPtr;
  size_t *indiciesOfInterest, *newIndiciesOfInterest;
  pair<f64, size_t> *sortSpace;

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

  tmpPtr = malloc(sizeof(*sortSpace) * size);
  sortSpace = (pair<f64, size_t>*) tmpPtr;

  tmpPtr = malloc(sizeof(*newIndiciesOfInterest) * size);
  newIndiciesOfInterest = (size_t*) tmpPtr;

  while(IOISize > 2){
    size_t NIOISize = 0;
    for(i = 0; i < IOISize-2; i+=2){

      _sortDoubleSizeTPairHighToLowHelper(toSort, indiciesOfInterest[i],
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


void _sortDoubleSizeTPairHighToLowHelper(pair<f64, size_t> *toSort,
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


void _sortDoubleSizeTPairLowToHigh(pair<f64, size_t> *toSort,
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
  
  size_t *indiciesOfInterest;
  tmpPtr = malloc(sizeof(*indiciesOfInterest) * size);
  indiciesOfInterest = (size_t*) tmpPtr;
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

  size_t *newIndiciesOfInterest;
  tmpPtr = malloc(sizeof(*newIndiciesOfInterest) * size);
  newIndiciesOfInterest = (size_t*) tmpPtr;
  
  while(IOISize > 2){
    size_t NIOISize = 0;
    for(i = 0; i < IOISize-2; i+=2){

      _sortDoubleSizeTPairLowToHighHelper(toSort, indiciesOfInterest[i],
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


void _sortDoubleSizeTPairLowToHighHelper(pair<f64, size_t> *toSort,
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

#endif
