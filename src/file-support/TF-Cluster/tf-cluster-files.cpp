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
#include <tf-cluster-files.hpp>
#include <statistics.h>
#include <quickmerge.hpp>

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
//FUNCTION DEFINITIONS//////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
ED loadExpressionData(const char *expressionFile,
                                           unordered_map<string, char> TFCheck){

  vector<string> colHeaders, rowHeaders;
  vector<char*> fileContents;
  string tmp, headerGeneName;
  char *word;
  //ifstream exprCoeffReader;
  FILE *exprCoeffReader;
  ED tr = {0, 0, NULL, 0, NULL, vector<string>()};
  void *tmpPtr;
  unordered_map<string, char> uniqueEntryCheck;
  string geneLabel;

  tmpPtr = calloc(sizeof(*word), 1024);
  word = (char*) tmpPtr;

  exprCoeffReader = fopen(expressionFile, "r");
  if(!exprCoeffReader){
    fprintf(stderr, STRING1, expressionFile);
    exit(1);
  }

  //Parse the initial data, and kind of tokenize it.
  int readStatus = 0;
  while (EOF != (readStatus = fscanf(exprCoeffReader, "%s1022", word))){
    tmpPtr = realloc(word, strlen(word)+1);
    word = (char*) tmpPtr;
    if(!strcmp(word, "")){
      free(word);
      continue;
    }
    if(word[0] == '#'){
      free(word);
      continue;
    }
    fileContents.push_back(word);
    tmpPtr = calloc(sizeof(*word), 1024);
    word = (char*) tmpPtr;
  }
  fclose(exprCoeffReader);
  free(word);

  //See how many experiments/samples we can expect per variable (gene)
  tr.corrVecLeng = 0;
  errno = 0;
  while(0.0 != strtod(fileContents[1+tr.corrVecLeng], NULL)){
    tr.corrVecLeng++;
  }

  //Check that the diamensions can work (but this isn't fool proof).
  //TODO: make fool proof
  if(fileContents.size() % (tr.corrVecLeng+1) != 0){
    fprintf(stderr, "Experimental data in file has a dataset length "
                                                    "inconsistancy\n");
    //exit(EINVAL);
  }
  tr.numGenes = fileContents.size() / (tr.corrVecLeng+1);

  //allocate memory to store experimental data
  tr.exprData = (f64**) malloc(sizeof(f64*) * tr.numGenes);
  for(size_t i = 0; i < tr.numGenes; i++)
    tr.exprData[i] = (f64*) malloc(sizeof(*tr.exprData) * tr.corrVecLeng);

  //Start parsing in experimental data
  size_t entryIndex = 0;
  tmpPtr = malloc(sizeof(*tr.TFCorrData) * TFCheck.size());
  tr.TFCorrData = (size_t*) tmpPtr;
  tr.TFCorrDataLength = 0;
  for(size_t rowItr = 0; rowItr < tr.numGenes; rowItr++){
    geneLabel = string(fileContents[entryIndex++]);
    free(fileContents[entryIndex-1]);

    //Make sure it isn't redundant (actually a huge issue)
    if(uniqueEntryCheck.count(geneLabel)){
      fprintf(stderr, STRING2, fileContents[entryIndex-1]);
      entryIndex += tr.corrVecLeng;
      tr.numGenes--;
      continue;
    }

    uniqueEntryCheck.insert(pair<string, char>(geneLabel, '\0'));
    tr.genes[rowItr] = strdup(geneLabel.c_str());//TODO: make optimal

    //read in samples
    for(size_t j = 0; j < tr.corrVecLeng; j++){
      tr.exprData[rowItr][j] = strtod(fileContents[entryIndex++], NULL);
      free(fileContents[entryIndex-1]);
    }

    //if this is a Transcription Factor we want to focus on, remember it
    if(TFCheck.count(geneLabel)){
      tr.TFCorrData[tr.TFCorrDataLength++] = rowItr;
    }

  }
  tr.TFCorrData = (std::size_t*)realloc(tr.TFCorrData, sizeof(*tr.TFCorrData) * tr.TFCorrDataLength);

  return tr;
}


unordered_map<string, char> readInTFList(const char *TFListFile){
  string geneName;
  FILE *geneListReader;
  unordered_map<string, char> tr;
  char *word;
  void *tmpPtr;

  errno = 0;
  geneListReader = fopen(TFListFile, "r");
  if(!geneListReader || errno){
    fprintf(stderr, STRING1, TFListFile);
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


////////////////////////////////////////////////////////////////////////
//END///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
