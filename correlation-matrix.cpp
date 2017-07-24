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

#include <cctype>
#include <cerrno>
#include <cmath>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <math.h>
#include <pthread.h>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>

#include "correlation-matrix.hpp"
#include "statistics.h"
#include "quickmerge.hpp"

////////////////////////////////////////////////////////////////////////
//DEFINES///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#define STRING1 "Error opening correlation file %s for reading\n"
#define STRING2 "WARNING: redundant data entry for %s, omitting from " \
                "analysis -- clean up your data!\n"
#define STRING3 "WARNING: missing gene data for Transcription Factor " \
                "%s, omitting from analysis -- clean up your data!\n"
 

////////////////////////////////////////////////////////////////////////
//NAMESPACE USING///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

using std::abs;
using std::ifstream;
using std::iostream;
using std::make_pair;
using std::ofstream;
using std::pair;
using std::make_pair;
using std::string;
using std::thread;
using std::unordered_map;
using std::vector;


////////////////////////////////////////////////////////////////////////
//STRUCTS///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

struct corrHelpStruct{
  size_t numerator;
  size_t denominator;
  double **results;
  const f64 **sumsOfSquares;
  const ED *input;
};

typedef struct corrHelpStruct CHS;


struct tauCorrHelpStruct{
  size_t numerator;
  size_t denominator;
  size_t vecLeng;
  double **results;
  const vector<size_t> *TFCorrData;
  const vector<GER> *geneCorrData;
};

typedef struct tauCorrHelpStruct TCHS;


struct rankHelpStruct{
  size_t numerator;
  size_t denominator;
  size_t vecLeng;
  vector<size_t> *TFCorrData;
  vector<GER> *geneCorrData;
};

typedef struct rankHelpStruct RHS;


struct rankCorrHelpStruct{
  size_t numerator;
  size_t denominator;
  size_t vecLeng;
  const vector<size_t> *TFCorrData;
  const vector<GER> *geneCorrData;
  double **results;
};

typedef struct rankCorrHelpStruct RCHS;


struct CAPHStruct{
    vector<size_t>* TFCorrData;
    vector<GER>* geneCorrData;
    size_t length;
    f64 **sumsOfSquares;
    size_t denominator;
    size_t numerator;
};

typedef struct CAPHStruct CAPHS;


struct corrRecord{
  size_t nameIndex;
  f64 *corrData;
  f64 sumOfSquares;
};

////////////////////////////////////////////////////////////////////////
//PRIVATE FUNCTION DECLARATIONS/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

/*******************************************************************//**
 *  Helper function to calculatePearsonCorrelationMatrix().
 **********************************************************************/
void *correlationHelper(void *protoArgs);


/*******************************************************************//**
 *  Header function which dispatches multiple threads to center row data
 * in fileData about their mean and pre-compute a sum of squares for
 * each row data.
 **********************************************************************/
f64* centerAndPrecompute(const ED &fileData);


/*******************************************************************//**
 *  Helper function to centerAndPrecompute().
 **********************************************************************/
void *centerAndPrecomputeHelper(void *arg);


//TODO:: add doc
void *rankHelper(void *protoArgs);

void *rankCorrHelper(void *protoArgs);


/*******************************************************************//**
 * Helper to calculate a coefficient matrix  for Kendall's tau
 * coefficient.
 **********************************************************************/
void *tauCorrelationHelper(void *protoArgs);


////////////////////////////////////////////////////////////////////////
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


ED loadExpressionData(cs8 *expressionFile, 
                            const unordered_map<string, char> &TFCheck){

  vector<string> colHeaders, rowHeaders;
  vector<char*> fileContents;
  string tmp, headerGeneName;
  char *word;
  //ifstream exprCoeffReader;
  FILE *exprCoeffReader;
  ED tr = {vector< size_t >(), vector< GER >(), 0, vector<string>(), 
                                                      vector<string>()};  
  void *tmpPtr;
  unordered_map<string, char> uniqueEntryCheck;
  string geneLabel;
  GER tmpRecord;
  
  tmpPtr = calloc(sizeof(*word), 1024);
  word = (char*) tmpPtr;

  exprCoeffReader = fopen(expressionFile, "r");
  if(!exprCoeffReader){
    fprintf(stderr, STRING1, expressionFile);
    exit(1);
  }
  
  int readStatus = 0;
  while (EOF != (readStatus = fscanf(exprCoeffReader, "%s1022", word))){
    tmpPtr = realloc(word, strlen(word)+1);
    word = (char*) tmpPtr;
    fileContents.push_back(word);
    tmpPtr = calloc(sizeof(*word), 1024);
    word = (char*) tmpPtr;
  }
  fclose(exprCoeffReader);
  free(word);
  
  tr.corrVecLeng = 0;
  errno = 0;
  while(0.0 != strtod(fileContents[1+tr.corrVecLeng], NULL)){
    tr.corrVecLeng++;
  }
  
  
  if(fileContents.size() % (tr.corrVecLeng+1) != 0){
    fprintf(stderr, "Experimental data in file has a dataset length "
                                                    "inconsistancy\n");
    //exit(EINVAL);
  }
  

  //loop
  size_t entryIndex = 0;
  while(entryIndex < fileContents.size()){
    
    geneLabel = string(fileContents[entryIndex++]);
    
    if(geneLabel == string("")) continue;

    if(uniqueEntryCheck.count(geneLabel)){
      fprintf(stderr, STRING2, fileContents[entryIndex-1]);
      entryIndex += tr.corrVecLeng;
      continue;
    }

    uniqueEntryCheck.insert(pair<string, char>(geneLabel, '\0'));
    tmpRecord.nameIndex = tr.genes.size();
    tr.genes.push_back(geneLabel);

    
    tmpPtr = malloc(sizeof(f64) * tr.corrVecLeng);
    tmpRecord.exprData = (f64*) tmpPtr;
    
    for(size_t i = 0; i < tr.corrVecLeng; i++)
      tmpRecord.exprData[i] = strtod(fileContents[entryIndex++], NULL);

    tr.geneCorrData.push_back(tmpRecord);

    if(TFCheck.count(geneLabel)){
      tr.TFCorrData.push_back(tr.geneCorrData.size()-1);
      tr.TFs.push_back(geneLabel);
    }
    
  }
  
  for(size_t i = 0; i < fileContents.size(); i++){
    free(fileContents[i]);
  }

  for(auto i = TFCheck.begin(); i != TFCheck.end(); i++)
    if(!uniqueEntryCheck.count(i->first))
      fprintf(stderr, STRING3, i->first.c_str());

  return tr;
}


/***********************************************************************
From expression data, construct a upper-diagonal section of a
correlation matrix, omitting the x=y entries.
***********************************************************************/
CMF calculatePearsonCorrelationMatrix(const ED &input, 
                                                  cf64 **sumsOfSquares){

  CMF tr;
  pthread_t *workers;
  struct corrHelpStruct *instructions;
  int *toIgnore;
  void *tmpPtr;
  
  tr.TFLabels = input.TFs;
  tr.GeneLabels = input.genes;

  tmpPtr = malloc(sizeof(*tr.fullMatrix) * input.TFCorrData.size());
  tr.fullMatrix = (f64**) tmpPtr; //[TF index][gene index]
  for(size_t i = 0; i < input.TFCorrData.size(); i++){
    tmpPtr = malloc(sizeof(**tr.fullMatrix)* input.geneCorrData.size());
    tr.fullMatrix[i] = (f64*) tmpPtr;
  }

  //* TODO NOTE: this is a very GPU friendly set of operations -- OpenCL
  csize_t numCPUs = thread::hardware_concurrency() < tr.numRows() ?
                    thread::hardware_concurrency() : tr.numRows() ;
  //*/

  /*
  csize_t numCPUs = 1;  //*/

  tmpPtr = malloc(sizeof(*workers) * numCPUs);
  workers = (pthread_t*) tmpPtr;
  tmpPtr = malloc(sizeof(*instructions) * numCPUs);
  instructions = (struct corrHelpStruct*) tmpPtr;


  for(size_t i = 0; i < numCPUs; i++){
    instructions[i] = {i,
                       numCPUs,
                       tr.fullMatrix,
                       sumsOfSquares,
                       &input
                      };
  }

  if(numCPUs > 1){
    //*
    for(size_t i = 0; i < numCPUs; i++)
      pthread_create(&workers[i], NULL, correlationHelper,
                                                      &instructions[i]);

    for(size_t i = 0; i < numCPUs; i++)
      pthread_join(workers[i], (void**) &toIgnore);
    //*/
  }else{
    correlationHelper((void*) instructions);
  }

  free(workers);
  free(instructions);

  return tr;
}


/***********************************************************************
From expression data, construct a upper-diagonal section of a
correlation matrix, omitting the x=y entries.
***********************************************************************/
void calculateRankCorrelationMatrix(ED &input){

  pthread_t *workers;
  RHS *RHSinstructions;
  //RCHS *RCHSinstructions;
  int *toIgnore;
  void *tmpPtr;

  //* TODO NOTE: this is a very GPU friendly set of operations -- OpenCL
  csize_t numCPUs = thread::hardware_concurrency() < input.TFCorrData.size() ?
                    thread::hardware_concurrency() : input.TFCorrData.size() ;
  //*/

  /*
  csize_t numCPUs = 1;  //*/
  
  //Create the ranks////////////////////////////////////////////////////

  tmpPtr = malloc(sizeof(*workers) * numCPUs);
  workers = (pthread_t*) tmpPtr;
  tmpPtr = malloc(sizeof(*RHSinstructions) * numCPUs);
  RHSinstructions = (RHS*) tmpPtr;

  for(size_t i = 0; i < numCPUs; i++){
    RHSinstructions[i].numerator = i;
    RHSinstructions[i].denominator = numCPUs;
    RHSinstructions[i].vecLeng = input.corrVecLeng;
    RHSinstructions[i].TFCorrData = &input.TFCorrData;
    RHSinstructions[i].geneCorrData = &input.geneCorrData;
  }

  if(numCPUs > 1){
    for(size_t i = 0; i < numCPUs; i++)
      pthread_create(&workers[i], NULL, rankHelper, &RHSinstructions[i]);

    for(size_t i = 0; i < numCPUs; i++)
      pthread_join(workers[i], (void**) &toIgnore);
  }else{
    rankHelper((void*) RHSinstructions);
  }

  free(workers);
  free(RHSinstructions);
  
  //Calculate rank correlation//////////////////////////////////////////

  /*tmpPtr = malloc(sizeof(*workers) * numCPUs);
  workers = (pthread_t*) tmpPtr;
  tmpPtr = malloc(sizeof(*RCHSinstructions) * numCPUs);
  RCHSinstructions = (RCHS*) tmpPtr;

  for(size_t i = 0; i < numCPUs; i++){
    RCHSinstructions[i].numerator = i;
    RCHSinstructions[i].denominator = numCPUs;
    RCHSinstructions[i].vecLeng = input.corrVecLeng;
    RCHSinstructions[i].TFCorrData = &input.TFCorrData;
    RCHSinstructions[i].geneCorrData = &input.geneCorrData;
    RCHSinstructions[i].results = tr.fullMatrix;
  }

  if(numCPUs > 1){
    for(size_t i = 0; i < numCPUs; i++)
      pthread_create(&workers[i], NULL, rankCorrHelper, 
                                                  &RCHSinstructions[i]);
    for(size_t i = 0; i < numCPUs; i++)
      pthread_join(workers[i], (void**) &toIgnore);
  }else{
    correlationHelper((void*) RCHSinstructions);
  }

  free(workers);
  free(RCHSinstructions);*/

}


void *correlationHelper(void *protoArgs){
  const struct corrHelpStruct *args = (struct corrHelpStruct*)
                                                              protoArgs;

  csize_t numerator = args->numerator;
  csize_t denominator = args->denominator;
  csize_t corrVecLeng = args->input->corrVecLeng;
  f64 **results = args->results;
  const vector<size_t> *TFCorrData = &(args->input->TFCorrData);
  const vector<GER> *geneCorrData = &(args->input->geneCorrData);
  cf64 **sumsOfSquares = args->sumsOfSquares;
  csize_t numTFs = TFCorrData->size();
  csize_t numGenes = geneCorrData->size();


  for(size_t y = (numTFs * numerator) / denominator; 
                      y < (numTFs * (numerator+1)) / denominator; y++){
    for(size_t x = 0; x < numGenes; x++){
      cf64 abCrossSum = getSumOfMultipliedArrays(
                        (*geneCorrData)[(*TFCorrData)[y]].exprData, 
                        (*geneCorrData)[x].exprData, corrVecLeng);
      results[y][x] = getCenteredCorrelationBasic(sumsOfSquares[0][y],
                                      sumsOfSquares[1][x], abCrossSum);
    }
  }

  return NULL;
}


void *tauCorrelationHelper(void *protoArgs){
  const struct tauCorrHelpStruct *args = (struct tauCorrHelpStruct*)
                                                              protoArgs;

  csize_t numerator = args->numerator;
  csize_t denominator = args->denominator;
  csize_t corrVecLeng = args->vecLeng;
  f64 **results = args->results;
  const vector<size_t> *TFCorrData = args->TFCorrData;
  const vector<GER> *geneCorrData = args->geneCorrData;
  csize_t numTFs = TFCorrData->size();
  csize_t numGenes = geneCorrData->size();


  for(size_t y = (numTFs * numerator) / denominator; 
                      y < (numTFs * (numerator+1)) / denominator; y++){
    for(size_t x = 0; x < numGenes; x++){
      ssize_t coordinateDisccordinatePairTally = 0;
      for(size_t i = 0; i < corrVecLeng; i++){
        double left, right;
        left = (*geneCorrData)[x].exprData[i];
        right = (*geneCorrData)[(*TFCorrData)[y]].exprData[i];
        if(left == right)
          coordinateDisccordinatePairTally++;
        else
          coordinateDisccordinatePairTally--;
      }
      coordinateDisccordinatePairTally *= 2;
      results[y][x] = ((f64) coordinateDisccordinatePairTally) / 
                      (corrVecLeng*(corrVecLeng-1));
    }
  }

  return NULL;
}


void *rankHelper(void *protoArgs){
  const RHS *args = (RHS*) protoArgs;
  void *tmpPtr;
  pair<f64, size_t> *toSort;
  
  csize_t numerator = args->numerator;
  csize_t denominator = args->denominator;
  csize_t corrVecLeng = args->vecLeng;
  vector<GER> *geneCorrData = args->geneCorrData;
  
  csize_t numGenes = geneCorrData->size();

  tmpPtr = malloc(sizeof(*toSort) * corrVecLeng);
  toSort = (pair<f64, size_t>*) tmpPtr;

  for(size_t i = (numGenes * numerator) / denominator; 
                    i < (numGenes * (numerator+1)) / denominator; i++){
    for(size_t j = 0; j < corrVecLeng; j++){
      toSort[j] = pair<f64, size_t>((*geneCorrData)[i].exprData[j],
                                                                    j);
    }
    _sortDoubleSizeTPairLowToHigh(toSort, corrVecLeng);
    for(size_t j = 0; j < corrVecLeng; j++){
      (*geneCorrData)[i].exprData[toSort[j].second] = j;
    }
  }
  
  free(toSort);

  return NULL;
}


f64** centerAndPrecompute(ED &fileData){
  f64 **tr; //sumsOfSquares
  pthread_t *workers;
  struct CAPHStruct *instructions;
  int *toIgnore;
  void *tmpPtr;

  tr = (f64**) malloc(sizeof(*tr) * 2);
  tr[0] = (f64*) malloc(sizeof(**tr) * fileData.TFCorrData.size());
  tr[1] = (f64*) malloc(sizeof(**tr) * fileData.geneCorrData.size());

  //TODO NOTE: this is a very GPU friendly set of operations -- OpenCL
  size_t numCPUs;
  csize_t a = thread::hardware_concurrency();
  csize_t b = fileData.TFCorrData.size();
  csize_t c = fileData.geneCorrData.size();
  if(a < b && a < c)  numCPUs = a;
  else if(b < c)      numCPUs = b;
  else                numCPUs = c;

  /*
  numCPUs = 1;//*/

  tmpPtr = malloc(sizeof(*workers) * numCPUs);
  workers = (pthread_t*) tmpPtr;
  tmpPtr = malloc(sizeof(*instructions) * numCPUs);
  instructions = (struct CAPHStruct*) tmpPtr;


  for(size_t i = 0; i < numCPUs; i++){
    instructions[i].TFCorrData = &(fileData.TFCorrData);
    instructions[i].geneCorrData = &(fileData.geneCorrData);
    instructions[i].length = fileData.corrVecLeng;
    instructions[i].sumsOfSquares = tr;
    instructions[i].denominator = numCPUs;
    instructions[i].numerator = i;
  }


  if(numCPUs > 1){
    for(size_t i = 0; i < numCPUs; i++)
      pthread_create(&workers[i], NULL, centerAndPrecomputeHelper,
                                                      &instructions[i]);

    for(size_t i = 0; i < numCPUs; i++)
      pthread_join(workers[i], (void**) &toIgnore);
  }else{
    centerAndPrecomputeHelper((void*) instructions);
  }
  
  for(size_t i = 0; i < fileData.TFCorrData.size(); i++){
    tr[0][i] = tr[1][fileData.TFCorrData[i]];
  }

  free(workers);
  free(instructions);

  return tr;
}


void *centerAndPrecomputeHelper(void *arg){
  struct CAPHStruct *args = (struct CAPHStruct*) arg;
  vector<GER> *geneCorrData = args->geneCorrData;
  csize_t length = args->length;
  f64 **sumsOfSquares = args->sumsOfSquares;
  csize_t denominator = args->denominator;
  csize_t numerator = args->numerator;

  for(size_t i = (geneCorrData->size() * numerator) / denominator;
        i < (geneCorrData->size() * (numerator+1)) / denominator; i++){
    inplaceCenterMean((*geneCorrData)[i].exprData, length);
    sumsOfSquares[1][i] = getSumOfSquares(
                              (*geneCorrData)[i].exprData, length);
  }

  return NULL;
}


unordered_map<string, char> TFList(const char *geneList){
  string geneName;
  FILE *geneListReader;
  unordered_map<string, char> tr;
  char *word;
  void *tmpPtr;

  errno = 0;
  geneListReader = fopen(geneList, "r");
  if(!geneListReader || errno){
    fprintf(stderr, STRING1, geneList);
    exit(1);
  }

  tmpPtr = malloc(sizeof(*word) * 1024);
  word = (char*) tmpPtr;
  while(!feof(geneListReader)){
    memset(word, 0, sizeof(*word) * 1024);
    fscanf(geneListReader, "%s1022", word);
    geneName = string(word);
    if(geneName == string("")) continue;
    if(tr.count(geneName)){
      fprintf(stderr, STRING2, geneName.c_str());
    }else{
      tr.insert(pair<string, char>(geneName, '\0'));
    }
  }
  free(word);

  fclose(geneListReader);

  return tr;
}


CMF calculateKendallsTauCorrelationCorrelationMatrix(const ED &input){
  
  CMF tr;
  pthread_t *workers;
  struct tauCorrHelpStruct *instructions;
  int *toIgnore;
  void *tmpPtr;
  
  tr.TFLabels = input.TFs;
  tr.GeneLabels = input.genes;

  tmpPtr = malloc(sizeof(*tr.fullMatrix) * input.TFCorrData.size());
  tr.fullMatrix = (f64**) tmpPtr; //[TF index][gene index]
  for(size_t i = 0; i < input.TFCorrData.size(); i++){
    tmpPtr = malloc(sizeof(**tr.fullMatrix) * input.geneCorrData.size());
    tr.fullMatrix[i] = (f64*) tmpPtr;
  }

  //* TODO NOTE: this is a very GPU friendly set of operations -- OpenCL
  csize_t numCPUs = thread::hardware_concurrency() < tr.numRows() ?
                    thread::hardware_concurrency() : tr.numRows() ;
  //*/

  /*
  csize_t numCPUs = 1;  //*/

  tmpPtr = malloc(sizeof(*workers) * numCPUs);
  workers = (pthread_t*) tmpPtr;
  tmpPtr = malloc(sizeof(*instructions) * numCPUs);
  instructions = (struct tauCorrHelpStruct*) tmpPtr;

  for(size_t i = 0; i < numCPUs; i++){
    instructions[i].numerator = i;
    instructions[i].denominator = numCPUs;
    instructions[i].vecLeng = input.corrVecLeng;
    instructions[i].results = tr.fullMatrix;
    instructions[i].TFCorrData = &input.TFCorrData;
    instructions[i].geneCorrData = &input.geneCorrData;
  }

  if(numCPUs > 1){
    for(size_t i = 0; i < numCPUs; i++)
      pthread_create(&workers[i], NULL, tauCorrelationHelper,
                                                      &instructions[i]);

    for(size_t i = 0; i < numCPUs; i++)
      pthread_join(workers[i], (void**) &toIgnore);
  }else{
    tauCorrelationHelper((void*) instructions);
  }

  free(workers);
  free(instructions);

  return tr;
}


extern CMF generateMatrixFromFile(cs8 *expressionFile, cs8 *geneList, 
                                                  cs8 *correlationType){
  CMF tr;
  ED fileData;
  f64 **sumsOfSquares;
  unordered_map<string, char> TFCheck;
  s8 *corrTypeLowerCase;

  corrTypeLowerCase = (s8*) malloc(strlen(correlationType)+1);
  for(size_t i = 0; i < strlen(correlationType); i++)
    corrTypeLowerCase[i] = tolower(correlationType[i]);
  corrTypeLowerCase[strlen(correlationType)]=0;
  
  TFCheck = TFList(geneList);
  fileData = loadExpressionData(expressionFile, TFCheck);
  
  if(!strcmp(corrTypeLowerCase, "pearson")){
    sumsOfSquares =  centerAndPrecompute(fileData);
    tr = calculatePearsonCorrelationMatrix(fileData, 
                                                (cf64**) sumsOfSquares);
    free(sumsOfSquares[0]);
    free(sumsOfSquares[1]);
    free(sumsOfSquares);
  }else if(!strcmp(corrTypeLowerCase, "spearman")){
    calculateRankCorrelationMatrix(fileData);
    //The difference between spearman and pearson
    sumsOfSquares =  centerAndPrecompute(fileData);
    tr = calculatePearsonCorrelationMatrix(fileData, 
                                                (cf64**) sumsOfSquares);
    free(sumsOfSquares[0]);
    free(sumsOfSquares[1]);
    free(sumsOfSquares);
  }else if(!strcmp(corrTypeLowerCase, "kendall")){
    calculateRankCorrelationMatrix(fileData);
    tr = calculateKendallsTauCorrelationCorrelationMatrix(fileData);
  }else{
    errno = EINVAL;
    tr.fullMatrix = NULL;
  }
  free(corrTypeLowerCase);
  
  for(size_t i = 0; i < fileData.geneCorrData.size(); i++)
                          free(fileData.geneCorrData[i].exprData);

  return tr;
}

////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
